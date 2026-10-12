# Published PT programs along the package pseudospace: notebooks 79 (human vs mouse) and 80 (AKI)

Frozen on 2026-10-12, with the coordinator's decisions of the same day built in (section 12), before any fit on our data.
Nothing has been computed on our data. The program lists were fixed from the published tables alone. Our data entered
only in two ways: gene membership of the notebook 71 and 76 universes (the "measured in our panel" rule), and the
top-100 cap, which ranks genes by each paper's own statistic.

**Frozen files** (sha256 in section 1.7). They are committed under `docs/results/known-programs/` and are the source of
truth. The notebooks read them relative to REPO and assert their hashes.
- `programs_hvm.csv` and `programs_aki.csv`: one row per program, with its members and the per-program gene
  accounting (non-one-to-one, unmapped symbol, not in universe, label/orientation, capped);
- `programs_{hvm,aki}_members.csv`: one row per program × member, with the paper's statistic;
- `programs_{hvm,aki}_dropped.csv`: programs below the size rule.

**The builder,** `build_programs.py --scratch <scratch> --repo <repo>`, stays in the session scratch, because it reads
the scratch lit-full extraction (`<scratch>/litfull/tables`). It holds no absolute path.

## Why

The user's idea: with two specimens per group, a 1,513-set discovery screen cannot support new pathway claims. The
record says the same thing:
- The AKI screens fail calibration: 0.120 in notebook 71 and 0.144 on notebook 76's universe.
- The human vs mouse screen passes narrowly (0.071 against the 0.075 limit).
- Notebook 77's 83 calls replicate in atlases no better than random genes of the same positional strength (p 0.35 and
  0.64). Replication supports genes, not gene sets.

So go hypothesis-driven. Take what atlases and microdissection studies already published (human vs mouse PT
differences, PT zonation programs, injured-PT programs after ischaemia–reperfusion), fix those programs before looking
at our data, and ask a small set of questions of each:
1. Does our cohort reproduce the published claim?
2. Where along the PT is the difference: which segment, uniform or segment-specific, and is that replicated across
   specimen pairs?
3. What does that add to the paper?

**Why position can add something here.** The published species claims are level claims within one segment, from
single-cell or single-nucleus DE with one species against the other two pooled. Our level comparison of human and
mouse is exposed to the two probe panels. Notebook 75's level control fails (18 of 27 literature controls have the
expected sign), and probe count leaks into the level (Spearman 0.46). A per-gene probe efficiency is a constant
factor along the PT, however, so it cancels in every between-segment contrast. The question "does the human–mouse
difference change between S1, S2 and S3?" is therefore the probe-invariant part of a species comparison. No published
table asks it: the lit-full check found no table that compares the S1→S3 profile between species. For AKI the panel is
shared, so level is not probe-exposed. The published injury lists name states, and sometimes segments (Kirita's
"Injured S1/2" and "Injured S3"), but none tests where along the PT the IR − control difference sits.

**What would make us abandon the claim that position adds information.** Suppose fewer than 2 of the 6 species
programs and fewer than 3 of the 9 injury programs reach a resolved localization class (localized, graded or uniform
within resolution; section 5.3). Then, at n = 2 + 2, the positional readout adds no resolvable information to the
published claims, and the notebooks say so.

**Bounds inherited from notebook 72.**
- C14 is not met: 10 of 12 human vs mouse specimen × segment cells, and 11 of 12 in AKI, are noise-dominated. Every
  location claim is therefore at segment level, and within-segment shape is descriptive.
- C8 is not met: the species registration gap is 0.37, so cross-species curves are read in segment units.
- Notebook 77's O6 applies: the coordinate-free (label-based) location is primary and the coordinate peak is a
  sensitivity.

## Labels

| | Notebook 79 | Notebook 80 |
|---|---|---|
| File | `analysis/notebooks/79_pt_known_species_programs.ipynb` | `analysis/notebooks/80_pt_known_injury_programs.ipynb` |
| Logic version | `79.known_species_programs.1` | `80.known_injury_programs.1` |
| Results | `results/pt_known_species_programs/` | `results/pt_known_injury_programs/` |
| Comparison | human (case) vs mouse (reference) | IR (case) vs control (reference) |

Each results folder has `tables/`, `figures/`, `stage_cache/` and `run_record.json`. New code goes in
`pseudospace/program_localization.py`, with `tests/test_program_localization.py` (synthetic only). No existing
module, notebook or test changes.

## 1 · The programs (frozen before any look at our data)

### 1.1 What was read

- The lit-full extraction (`SCRATCH/litfull/tables/<source>__genes.csv`): each paper's own rows, statistic and
  threshold.
- The gene lists of notebook 71's human vs mouse universe (11,055 genes) and of notebook 76's U76 (10,701 genes).

No count, curve, statistic or call of our data was read.

### 1.2 Shared rules (both notebooks)

- **R1 · Symbols, with no alias mapping (decision O5).** A program member is the paper's gene mapped to its mouse
  symbol (`gene_mouse`). Human genes map through the project's one-to-one ortholog map, as the extraction did. No alias
  table is applied.
  - A mouse symbol that is neither in the ortholog map nor on the mouse probe panel is counted as `unmapped_symbol`
    and listed. Examples are older symbols in Kirita 2020 or Zhou 2023, such as Atp5a1 and Qk.
  - Range over the confirmatory programs: 0–26 in 79 and 0–36 in 80.
  - **Ponytail:** the upgrade path is an MGI synonym table.
- **R2 · One-to-one orthologs only (decision O3).** Excluded and counted as `non_one_to_one`:
  - a human gene without a one-to-one mouse ortholog;
  - a mouse gene on the probe panel that is not in the one-to-one map.

  Range over the confirmatory programs: 0–25 in 79 and 3–86 in 80, the most for human AKI.
- **R2b · Universe.** Members must be in the universe the gene models were fitted on:
  - 79: notebook 71's human vs mouse universe;
  - 80: notebook 76's U76 (one-to-one-ortholog mouse panel, AKI detection rule).

  One-to-one genes outside it are counted and listed (`not_in_universe`), and not tested.
- **R3 · Label and orientation genes out.** The 15 genes that name the package's PT segments or orient its curve are
  removed from every program and listed: the PT-S1/S2/S3 label panel (Cyp2e1, Cyp7b1, Gatm, Slc7a13, Slc13a3, Slc22a6,
  Slc5a2, Slc22a7, Slc5a12) and the axis anchors (also Aqp2, Aqp4, Atp6v0d2, Rhcg, Slc14a2, Wnt7b). This is notebook
  72's evaluation exclusion and 77's K10. Their zonation is circular, and their species × segment contrast is forced
  towards zero by the labelling.
- **R4 · Cap.** If more than 100 members remain, the top 100 by the paper's own effect statistic in the published
  direction are kept (ties by name). CellSpectra blocks are already a published top 99.
- **R5 · Size.** A program needs at least 5 members after R1–R4, or it is dropped and recorded.
- **R6 · Roles.**
  - *Confirmatory* programs enter a BH family (section 6).
  - *Secondary* programs get the same analyses with nominal p, labelled secondary. They are other sources for the same
    claim, or segment-resolved sources too narrow or too different from our design to be primary.

### 1.3 Notebook 79, species family S (confirmatory, 6 programs)

**Source.** CellSpectra / SISKA 1.0, Klötzer et al., Nat Genet 2025, Supplementary Table 4
([PMID 40775269](https://pubmed.ncbi.nlm.nih.gov/40775269/)). It lists species DEGs per PT segment: in segment k, the
top 99 genes higher in one species than in the other two pooled (Wilcoxon, healthy and diseased cells pooled).

**Rule.**
- Block "Human PTSk" gives program S-Hk (human-high in Sk). Block "Mouse PTSk" gives S-Mk (mouse-high in Sk).
- A gene listed in both blocks of the same segment is removed, because such a listing implies no human-vs-mouse
  direction (rat low). This removes 4 genes in S3 (Kcnj15, Mast4, Pbx1, Slc16a9) and none in S1 or S2.

| id | published claim (one line) | direction | segment | n listed → eligible → kept | first members |
|---|---|---|---|---|---|
| S-H1 | higher in human than in rodents in PTS1 | human-high | S1 | 99 → 95 → 95 | Cubn, Dlg2, L3mbtl4, Slc2a9, Mapk10, Ptprm |
| S-H2 | higher in human than in rodents in PTS2 | human-high | S2 | 99 → 93 → 93 | Cubn, L3mbtl4, Slc2a9, Dlg2, Nlgn1, Ptprm |
| S-H3 | higher in human than in rodents in PTS3 | human-high | S3 | 95 → 90 → 90 | Cubn, Nlgn1, L3mbtl4, Insr, Farp1, Cntn4 |
| S-M1 | higher in mouse than in human and rat in PTS1 | mouse-high | S1 | 99 → 96 → 96 | Fut9, Dab2, Esrrg, Ptprd, Osbpl6, Pde4d |
| S-M2 | higher in mouse than in human and rat in PTS2 | mouse-high | S2 | 99 → 96 → 96 | Fut9, Dab2, Esrrg, Pde4d, Chpt1, Ghr |
| S-M3 | higher in mouse than in human and rat in PTS3 | mouse-high | S3 | 95 → 94 → 94 | Ghr, Pde4d, Fut9, Chpt1, Wdr17, Slc25a21 |

The blocks overlap across segments: S-H1 ∩ S-H2 = 72 genes and S-M1 ∩ S-M2 = 78. They are separate published claims,
so they stay separate programs. BH is valid under this positive dependence.

**Caveats that travel with these programs.**
- The comparator is the other two species pooled, not mouse alone.
- Disease and healthy cells are pooled.
- Only the top 99 genes per block are published, so absence from a block is not a claim.
- Human data in the atlas are mostly single-nucleus. The human-high blocks are rich in long, intron-rich genes (Dlg2,
  Nlgn1, Ptprm, Cntn4), a known single-nucleus signature. A failure to replicate them is compatible with a technology
  effect in the published contrast; our Visium HD probes are neither nuclear nor whole-cell. This is tested by the
  pre-specified gene-length sensitivity (section 4a), not only stated.

### 1.4 Notebook 79, zonation family Z (confirmatory, 6 programs)

| id | source (PMID) | rule | segment | n listed → eligible → kept | first members |
|---|---|---|---|---|---|
| Z-M1 | Chen 2021, microdissection, Supp Table 3 ([33769951](https://pubmed.ncbi.nlm.nih.gov/33769951/)) | up in PTS1 vs all other tubule segments (edgeR FDR < 0.05, \|log2FC\| ≥ 2) and not listed up for PTS2 or PTS3 | S1 | 223 → 173 → 100 (cap) | Ugt2b34, Slc34a3, Gc, Spp2, Cox6a2, Prodh2 |
| Z-M2 | same | the same rule for PTS2 | S2 | 288 → 252 → 100 (cap) | Dpys, Cyp4b1, Susd2, Azgp1, Slc5a10, Hpd |
| Z-M3 | same | the same rule for PTS3 | S3 | 109 → 68 → 68 | Akr1c18, Hsd11b1, Kif20b, Apobec2, Calml4, Chrdl1 |
| Z-H1 | Hinze 2022, Table S6 ([36085050](https://pubmed.ncbi.nlm.nih.gov/36085050/)) | human PT-S1 marker vs the other PT subclusters (Wilcoxon) | S1 | 50 → 39 → 39 | Pcdh15, Plppr1, Slc7a8, Rnf212b, Man1c1, Ppp2r2b |
| Z-H2 | same | PT-S2 marker vs the other PT subclusters | S2 | 14 → 6 → 6 | Sox5, Dach1, Slc34a1, Bhmt, Prkn, Cgnl1 |
| Z-H3 | same | PT-S3 marker vs the other PT subclusters | S3 | 103 → 61 → 61 | Gpm6a, Slc5a8, Nebl, Slc30a8, Zfp804b, Naaladl2 |

**Why these sources.**
- **Mouse: Chen 2021.** It is the brief's first-named source and the only anatomical (microdissected) ground truth,
  with no clustering. Its PT tables are vs_all, so a segment program takes the lit-full rule: listed for that PT
  segment and for no other PT segment.
- **Human: Hinze 2022, not Lake 2023 (decision O1, accepted).** Hinze's Table S6 is the only human within-PT table. Lake
  2023, which the brief names, publishes vs_all markers only. Under the same differential-listing rule its programs
  have 8 (S1), 3 (S2, fails R5) and 11 (S3) genes, so it becomes secondary.
- **Hinze caveat.** The comparator includes the AKI-associated PT states, and nuclei from AKI kidneys are pooled.

### 1.5 Notebook 80, injury family I (confirmatory, 9 programs)

**Rule for state programs** (Kirita et al., PNAS 2020, Dataset S2,
[PMID 32571916](https://pubmed.ncbi.nlm.nih.gov/32571916/)):
- Design: male mice, bilateral IRI, a time course from 4 h to 6 wk.
- Each injured PT state is compared against all other PT nuclei (padj < 0.05, positive markers).
- Markers of the healthy S1, S2 or S3 subclusters of the same table are removed, so each program carries
  injury-specific genes, not segment identity.

**Why Kirita is primary.** Among the eligible IRI sources it is the only one that has all of these:
- male mice and bilateral IRI;
- a within-PT comparator;
- named healthy S1/S2/S3 subclusters, which the subtraction needs;
- five named injured states in one table;
- a published segment identity for two of them.

**Core IRI response** (Zhou et al., JASN 2023, Dataset 3, [PMID 37639336](https://pubmed.ncbi.nlm.nih.gov/37639336/);
these are Balzer 2022's IRI data):
- The comparison is IRI vs control within PT, in 6 groups: short or long ischaemia × 1, 3 or 14 days.
- A gene is kept if listed in the same direction in at least 3 of the 6 groups and never in the opposite direction
  (decision O6, accepted).
- The rule is time-agnostic because our IRI time point is unknown.

**Human AKI transfer** (Hinze 2022, Table S2): AKI vs control kidneys, whole PT, pseudobulk DESeq2, padj < 0.001 and
|log2FC| > 1, mapped one-to-one to mouse.

| id | published claim (one line) | direction | segment | n listed → eligible → kept | first members |
|---|---|---|---|---|---|
| K-REP | Repairing PT: arises at 2 d; mitotic cell cycle (Top2a) | up in IR | PT | 152 → 118 → 100 (cap) | Top2a, Cdca2, Fyb2, Il7, Zwilch, Prkch |
| K-INJ12 | Injured S1/2: 4–12 h; Myc, Havcr1; shares healthy S1/S2 genes | up in IR | S1–S2 | 258 → 224 → 100 (cap) | Tmem108, Bmp6, Plin2, Gbe1, Col18a1, Nrg1 |
| K-INJ3 | Injured S3: 4–12 h; Myc, Havcr1; shares healthy S3 genes | up in IR | S3 | 472 → 419 → 100 (cap) | Anxa2, Havcr1, Rtn4, Lrrfip1, Cep85, Clic4 |
| K-SEV | Severe injured PT: Krt20, heat shock; cell-cycle arrest | up in IR | PT | 399 → 353 → 100 (cap) | Hspa1a, Hsp90aa1, Ubc, Hsp90ab1, Hspa8, Jun |
| K-FR | Failed-repair PT (FR-PTC): Vcam1, Dcdc2a; ~30% of PT at 14 d | up in IR | PT | 393 → 324 → 100 (cap) | Dock10, Pdgfd, Erbb4, Psd3, Vcam1, Sema5a |
| Z-IRI-UP | higher in IRI PT than control in ≥ 3 of 6 IRI groups | up in IR | PT | 405 → 324 → 100 (cap) | Spp1, Actg1, Mpc1, Tmsb4x, Mt1, S100a9 |
| Z-IRI-DN | lower in IRI PT than control in ≥ 3 of 6 IRI groups | down in IR | PT | 644 → 533 → 100 (cap) | Inmt, Slc34a1, Zbtb20, Cyp4b1, Timp3, Slc22a30 |
| H-AKI-UP | human AKI PT higher than control PT | up in IR | PT | 406 → 279 → 100 (cap) | Fosl1, Lcn2, Sfn, Map1b, Fga, Lrp8 |
| H-AKI-DN | human AKI PT lower than control PT | down in IR | PT | 307 → 155 → 100 (cap) | Afm, Ido2, Slc34a1, Mycn, Slc22a13, Dao |

**Published locations that our data can test directly.** These are stated before the run:
- K-INJ12 should differ most in S1–S2, and K-INJ3 in S3 (Kirita's segment identities).
- The other programs have no published segment, so their location is new information.
- Literature priors on location, for wording only, never for selection:
  - injured Krt20-lineage PT cells shift from the cortico-medullary boundary to the outer cortex between day 7 and
    day 28 (Gerhardt 2021, [34183416](https://pubmed.ncbi.nlm.nih.gov/34183416/));
  - Vcam1+/Ccl2+ late-injured cells sit "at the site of primary injury in the hypersensitive S3" (same paper);
  - notebook 76's priors (outer-stripe S3 injury; its cited PMIDs) apply unchanged.

### 1.6 Secondary programs (same analyses, nominal p, labelled secondary)

**Notebook 79** (zonation; the sources the brief also names):
- Ransick 2019 male clusters ([31689386](https://pubmed.ncbi.nlm.nih.gov/31689386/)): S1, S2 and S3, each 100 (cap).
- Zhou 2023 Dataset 6 within-PT ([37639336](https://pubmed.ncbi.nlm.nih.gov/37639336/)): S1, S2 and S3, each 100
  (cap). The comparator includes injured PT.
- Chen 2025 Extended Data Fig. 1 source data ([40259083](https://pubmed.ncbi.nlm.nih.gov/40259083/)): adult male strata
  W12, W52 and W92, differential listing pooled over those strata. S1, S2 and S3, each 100 (cap).
- Lake 2023 ([37468583](https://pubmed.ncbi.nlm.nih.gov/37468583/)): ST5 curated plus ST8 NSForest (snCv3); S1 8 and
  S3 11. S2 has 3 genes and is dropped.

**Notebook 80** (segment-resolved published injury results, and other sources of the same states):
- Gerhardt 2023, Supp Table 9 ([36735940](https://pubmed.ncbi.nlm.nih.gov/36735940/)):
  - repaired (Ki67-lineage) vs control nuclei within PTS1, PTS2 and PTS3, up and down: 6 programs of 69–100 genes;
  - this is late repair, at 4 wk and 6 mo;
  - they are the only published segment-resolved IRI-vs-control single-nucleus lists.
- Dixon 2022, Supp Data File 2 ([34853151](https://pubmed.ncbi.nlm.nih.gov/34853151/)):
  - Visium PTS1-2 and PTS3-TAL spots, IRI at 4 h, 12 h or 2 d vs sham, up and down: 4 programs of 20–100 genes;
  - female mice; the PTS3-TAL spots mix S3 and TAL.
- Same-state alternates:
  - Gerhardt 2023 FR-PTC (Supp Table 4);
  - Gerhardt 2021 cycling early-injured (cluster 10) and Vcam1+/Ccl2+ late-injured (cluster 11) PT
    ([34183416](https://pubmed.ncbi.nlm.nih.gov/34183416/));
  - Lake 2023 ST11, human adaptive PT up and down (100 each).

### 1.7 Considered and not used

| Candidate | Reason |
|---|---|
| CellSpectra injured-PT species blocks | our specimens are healthy |
| A CellSpectra "segment-reversing" program (human-high in one segment, mouse-high in another) | after removing same-segment double listings, only Setbp1 qualifies (human-high S1, mouse-high S3); it is shown as a gene |
| Fernandes 2026 SLC34A1/2/3 ([42494142](https://pubmed.ncbi.nlm.nih.gov/42494142/)) | 3 genes, below R5; shown as genes |
| Raghubar 2022 Visium cortex human vs mouse DE | whole-cortex resolution, not PT; confounded by age and depth by its authors' account |
| Chen 2023 PT subsegments by sex | supplementary tables not obtainable |
| Balzer 2022 Supp Data 9 | comparator not stated (group vs rest); the same data, against controls, are Zhou's Dataset 3 |
| Ide 2021, Polonsky 2024, Xuanyuan 2025, Zhou 2023 Injured1/2 | unilateral IRI or pooled across non-IRI models, or vs-all marker panels; not closer to our design than the primary lists |
| Lake 2023 ST12/ST13 | pan-cell-type degenerative/adaptive sets, not PT-specific |

**Frozen file hashes (sha256).** These are the files the coordinator copies to `docs/results/known-programs/`.
- `programs_hvm.csv`: `be61b8183c51d798ee33314e6226d5d8dd269188b82eb60ae65c888a355db79e`
- `programs_hvm_members.csv`: `2a3889075ca4130b102661409aee1390f1dd1f936a4171c9c9b7b0a7365299bf`
- `programs_hvm_dropped.csv`: `0671974d9ced8bb10007251579f744df17d07ee12bb72c82c2d0f75c7b1b58eb`
- `programs_aki.csv`: `1d01fe15a61f7333983c680e3c7cded1aa3ee4f5ab4cdcf1e7445f06a9949d40`
- `programs_aki_members.csv`: `e00750ad066a8a76b3c870efb6a4c96e56bc7163aa7830ca1bd0ed6df52c667d`
- `programs_aki_dropped.csv`: `b7e0cbe5f5659b62ca2c9ef2f6b8599c25097491fd4c89f9c2d63527b1fdaaad`

The notebooks assert these hashes on the committed copies (decision O2). Any later change needs a new hash and a
recorded reason.

### 1.8 Why the confirmatory families are small (6, 6 and 9)

- With two specimens per group, a program-level p has at best the resolution of the specimen-level variance.
- BH over a family of m pre-specified programs asks a raw p of about 0.05·i/m, so with m = 6–9 a single strong
  program can pass.
- A family of hundreds would rebuild the discovery screen that cannot be calibrated here. Notebook 76's AKI screen
  calibrates at 0.144.
- Every confirmatory program is one published claim, fixed before the data, so no selection enters the family.
- Each source contributes the programs that carry its published claim, and no more: species per direction × segment;
  zonation per species × segment from one primary source each; injury per named state from one primary source, plus
  the core IRI response and its human transfer.

## 2 · Inputs and guards

### Inputs (read, never refit, unless stated)

**Notebook 79.**
- Units: notebook 70's `package_pseudospace/hvm/pt_units.h5ad` through `package_inputs.hvm_inputs`, with the 13 `mt-`
  genes excluded (71's call). This gives the units, counts, mt-free library, segments, specimens, positions and
  universe.
- Notebook 71: the condition stage cache `gene_model_human_vs_mouse_species` (Δ(x), SE, group and specimen curves on
  the 101-point grid) and the relabeling caches `gene_model_human_vs_mouse_swap_*` (Δ(x)); also
  `gene_statistics_human_vs_mouse.csv` (positional F-test z and q per partition), `gene_sets_human_vs_mouse.json` (the
  1,513 library sets, for calibration only) and `universe_human_vs_mouse.csv`.
- Notebook 75: `gene_statistics_segment_human_vs_mouse.csv` (segment F-test z and q; Δ_S1–3 with HC3 SE, per
  partition) and `deseq2_genes_human_vs_mouse.csv` (the conventional comparator).
- Notebook 72: `transitions.csv` (pooled S1→S2 and S2→S3 positions, 77's rule).
- `pt_literature_deep/probe_counts.csv` (`probe_balance`), for the equal-probe sensitivity, as in 77.
- `data/external/ensembl_r110_human_transcripts/Homo_sapiens.GRCh38.cdna.all.fa.gz` (public; sha256
  `1881e80298f19207159b60bf1eb79b0c5fe34eb38c24c5b4e796a35f016c1e7d`), for the gene-length sensitivity (4a). Only
  the header lines are read.
- The frozen programs: `REPO/docs/results/known-programs/programs_hvm{,_members}.csv`, with their sha256 asserted.

**Notebook 80.**
- Units: notebook 70's AKI export through `package_inputs.aki_inputs`, with notebook 76's two calls (pool, then
  detection and `exclude_genes = MT`). This gives U76.
- Notebook 76: stage caches `gene_model_primary_{species,swap_*}` (Δ(x), specimen curves),
  `specimen_shape_primary_*`, `gene_statistics.csv` (positional F-test z and q on U76), `gene_sets.json` (1,495 library
  sets, for calibration), `segment_units.csv` (control cuts), and `placement_units.csv` with the run folder's
  `02_qc/units_qc.h5ad` counts (the c19 arm).
- Notebook 75: `gene_statistics_segment_aki_vs_control.csv` (guard only) and `deseq2_genes_aki_vs_control.csv`
  (comparator, on 71's AKI universe).
- The frozen programs: `REPO/docs/results/known-programs/programs_aki{,_members}.csv`, with their sha256 asserted.

### The one refit

80's segment F-test on U76. Notebook 75 fitted it on 71's AKI universe (10,198 genes), which lacks 503 U76 genes,
among them Havcr1, Lcn2, Mki67 and Top2a. 80 reruns `discrete_comparators.segment_f_test`, the same code, on U76 for
the condition and both relabelings, with 76's condition α. This is the brief's "fit with the same model code and say
so". 79 refits nothing.

### Guards (not cached)

- **G1.** The export contract and the confirmed pass-1 fingerprints: hvm `efafc9ac2e99`, AKI `bb4a80c1b484`.
- **G2.** The universes equal 71's (11,055 genes) and 76's (10,701 genes); the units are 71's (11,383) and 76's
  (14,048).
- **G3.** The cached condition models reproduce their gene tables: 71's T_spatial and z (75's G3 rule), and 76's
  `gene_statistics.csv`.
- **G4.**
  - 79: 75's segment table covers 71's universe.
  - 80: on the 10,198 genes U76 shares with 71, the refitted T_step, Δ_S and SE equal 75's (rtol 1e-6) wherever 76's
    α equals 71's. 76's G5 records two genes whose α differ by 1.4 × 10⁻⁵; they are reported, not a stop.
- **G5.** The program files match their sha256. Every member is in the universe, and no member is a label/orientation
  gene.
- **G6.** The segment scales reproduce their sources: 72's pooled transitions (79) and 76's control cuts (80).

## 3 · Per-gene statistics

Δ, all contrasts, case and reference follow the same convention: case − reference, so human − mouse and IR − control.
Curves are in natural log, and everything reported is in log2.

### 3.1 Reused, with no refit (except 80's segment F-test on U76)

| Statistic | 79 source | 80 source | Role |
|---|---|---|---|
| z_seg: segment F-test z (species × segment change, specimen-shape denominator) | 75 | refit on U76 (section 2) | between-segment test, coordinate-free (primary) |
| z_pos: positional F-test z (species × position change) | 71 | 76 | the coordinate version (sensitivity) |
| Δ_S1–3 (NB steps fit, equal specimen weights) | 75 | the U76 refit | check of the pseudobulk segment summaries |
| Δ(x), specimen curves, relabeled Δ(x) | 71 caches | 76 caches | program curves and figures |
| DESeq2 Wald stat per scope (whole PT, S1, S2, S3) | 75 | 75 (71's universe) | conventional comparator, descriptive |

### 3.2 New and coordinate-free: specimen × segment pseudobulks

The units are the 71 and 76 analysed units (the common-support trim removes about 2% per specimen at the coordinate
tails). The untrimmed version is a sensitivity (section 7).

- **Scopes:** S1, S2 and S3 (package labels), PT (all), and S1–S2 (S1 ∪ S2, for K-INJ12 and the Dixon PTS1-2 lists).
- **Rate:** r_gjs = log2((C_gjs + 0.5) / L_js), where C is the summed raw count of gene g over specimen j's units in
  scope s, and L is their summed mt-free library. The 0.5 floor is 71's floor rule.
- **Contrasts,** each linear in r and each centred by its median over the universe genes (76's centring, the
  positional analogue of median-of-ratios):
  - the **group difference** e_gs: the mean of the case specimens minus the mean of the reference specimens;
  - the **pair differences** for the 4 case × reference pairs;
  - the **within-group differences** d_case = case₁ − case₂ and d_ref = ref₁ − ref₂;
  - the **segment contrast** (zonation only): for specimen j and segment k, r_gjk minus the mean of r_gj over the other
    two segments.

### 3.3 New: moderated level t per scope and partition

This is limma's model on the four pseudobulks of a scope, in closed form.
- **Condition** (design ~ group): effect e = the centred group difference, s² = (d_case² + d_ref²)/4, 2 residual df.
  This is the variance rule of 71's S3 arm (there on a whole-PT level from curves), which passed its swaps (114 calls;
  0 and 0). Here it is applied per scope to pseudobulks.
- **Relabeling** (design ~ species + relabel; 1 residual df):
  - relabeling 1 has effect (d_case + d_ref)/2 and s² = (d_case − d_ref)²/4;
  - relabeling 2 has effect (d_case − d_ref)/2 and s² = (d_case + d_ref)²/4.

  The case and reference roles follow `partition_groups` (hvm: {HUK1_COR1, Ctrl1A2} and {HUK1_COR1, Ctrl1A4};
  AKI: {IR2A2, Ctrl1A2} and {IR2A2, Ctrl1A4}).
- **Moderation:** `pathway_pipelines.squeeze_var` (covariate log mean raw count per unit, robust FALSE), then
  t = e / √s²_post on df + d₀, with a two-sided p.

### 3.4 New: moderated within-species segment contrast (zonation family only)

- For one species and segment k: the effect is the mean over its 2 specimens of the centred segment contrast, and
  s² = (c₁ − c₂)²/4 on 1 df.
- squeeze_var as in 3.3, then t on 1 + d₀.
- It is computed in each species for each segment.

### 3.5 Program curves (coordinate version)

- The curves follow 76's rule: `aki_along_pt.centre_rows` on Δ(x) (the median gene at each x), then
  `program_profile`, the median over usable members, in log2. Usable means converged and not separated in the
  condition fit.
- **Pair curves:** each case specimen curve minus each reference specimen curve, centred the same way; NaN outside a
  specimen's own support.
- **Relabeled curves:** the relabelings' Δ(x), centred.
- **Segment units.** 79 uses 72's pooled transitions (77's reading 2: HUK1_MED1's transitions cross, so the
  species-specific mappings are undefined). 80 uses 76's control cuts (t₁₂ 0.349, t₂₃ 0.643).
- **Median primary, mean as a sensitivity (decision O4).** The primary program curve, and every segment summary, is
  the median over the members' log-ratio curves (76's rule). Human-only or mouse-only members reach floor-limited log
  ratios of 5–10, which would dominate a mean. The mean is a pre-specified sensitivity (section 7).

## 4 · Analyses (per program)

### 4a · Replication of the published claim

The claim's scope is its published segment: S1, S2, S3, S1–S2, or PT when no segment is published.

**Species and injury programs (families S and I).**
- **Gene level** (descriptive):
  - the share of members whose centred group difference at the claim's scope has the published sign;
  - the share with that sign in all 4 pairs;
  - both next to the universe base rate.
- **Program level** (the test of the claim): `camera_pr` (inter.gene.cor 0.01) on the moderated level t at the
  claim's scope, with the family's programs as the sets.
  - The one-sided p in the claimed direction is p/2 when cameraPR's direction matches the claim, and 1 − p/2
    otherwise.
  - The opposite-direction one-sided p is also kept.
- **Calibration of the level statistic.**
  - cameraPR of the same statistic runs on all library sets (79: 1,513; 80: 1,495) for the condition and both
    relabelings, per scope.
  - Notebook 75's `swap_rule` is applied (calibration p = the relabelings' two-sided p, pooled).
  - If the calibration share at that scope is above 0.075, the program's p is replaced by the empirical p against
    that pooled null (`aki_along_pt.empirical_p` on two-sided p, then halved in the claimed direction as above).
  - The swap-rule table is an output.
- **Specificity:** the program's own relabeled one-sided p is above 0.05 in both relabelings (76's set-claim rule).
- **The verdict uses all program genes (decision O8).** Two pre-specified sensitivities rerun the whole 5.1 rule. A
  program whose verdict (class) changes under either one carries that label beside its class.
  - **Equal probe counts (79):** members and background restricted to genes with `probe_balance` equal (77's
    subset), when at least 5 members remain. Label: "probe-sensitive".
  - **Gene length (79, family S; the CellSpectra single-nucleus caveat):**
    - Gene length is the human genomic span of the gene: max(end) − min(start) + 1 over its transcripts in the
      Ensembl 110 GRCh38 cDNA headers, matched to the export's `human_symbol`, case-insensitive. A loading check found
      every universe gene covered, and parsing takes about 3 s.
    - The level t is rank-normalised within deciles of log₁₀ span over the universe genes, then the same cameraPR
      test runs. Length effects on level then cancel within strata.
    - The members' median span against the universe's is reported.
    - Label: "length-sensitive".
    - The between-segment test is not length-exposed (a per-gene constant cancels), so the stratified version is
      reported for it as descriptive only.
- **Conventional comparator** (descriptive): cameraPR on 75's DESeq2 Wald stat at the claim's scope (condition and
  relabelings).

**Zonation programs (family Z).**
- **Source species:** `camera_pr` (Up) on the within-species t (3.4) for the program's segment, in the species the
  claim comes from.
- **Transfer:** the same test in the other species.
- **Segment contrast medians** per specimen are reported, together with the genome-wide ratio of human to mouse
  contrast for the same segment. That ratio is the global flattening of C12 and 77, so a program's attenuation can be
  read beyond it.

### 4b · Localization along the PT

This is the human − mouse difference for families S and Z, and the IR − control difference for family I.

**Coordinate-free (primary).**
- **Segment summaries:** P_k is the median over members of the centred group difference in scope k (S1, S2, S3), for
  the condition and for each of the 4 pairs. The relabeled medians come from the relabeled effects of 3.3.
- **Signed values:** Q_k = d · P_k, where d is the claimed direction. For zonation programs, d is the sign of the
  condition P at the program's own segment.
- **Between-segment test:** `camera_pr` (Up) of the program on z_seg against all universe genes. The p is the one-sided
  Up p.
  - 79: z_seg passes 75's swap rule (0.055), so the p is nominal.
  - 80: 75's AKI segment F-test fails calibration (0.115). The swap rule is recomputed on U76; if it still fails, the
    p is empirical against the relabeled one-sided p of the library sets (76's rule).
- **Peak segment:** argmax_k Q_k, in the condition and in each pair.
- **Margin:** Q at the peak minus the next-largest Q (condition).
- **Relabeled ratio:** |relabeled Q at the peak| / |condition Q at the peak|, for both relabelings.

**Coordinate version (sensitivity and display).**
- The program curve, 4 pair curves and 2 relabeled curves (3.5).
- `profile_summary`: level, drift, centre of mass, peak in segment units, and the change within each segment.
- The between-position test: `camera_pr` (Up) on z_pos. Its p is calibrated as for z_seg (hvm passes at 0.071; AKI
  fails, so empirical).
- **Axis agrees:** the floor of the curve's peak (segment units) equals the label-based peak segment.
- Within-segment change is reported, labelled "descriptive (C14)".

**Check:** the label-based peak segment against the peak of the median of 75's Δ_S (79), or of the refit's Δ_S (80).

### 4c · Added information (one line per program)

Each line comes from a fixed template filled from sections 4a, 4b and 5, so no line is hand-written after the run:

> "{source}: {claim}. Ours: {replication class} at {scope} (median {P} log2; {n}/4 pairs; q {q}). Along the PT:
> {localization class}{, largest in S_k}{; S_j: no difference resolved}. Axis {agrees | differs}. {flags}."

The flags are: probe-sensitive; human S2 small (632 units); S3 compares human cortical S3 with mouse S3 including the
outer stripe; label limitation (AKI c15/c19); technology (CellSpectra).

Zonation lines add the transfer class: "published as mouse S3 program; conserved / attenuated / absent / reversed in
human; the species difference sits in S_k / is uniform".

### 4d · Figures

**Per family, small multiples.** One row per program, three panels:
1. **Program curve.** Centred Δ(x) median in log2 on segment units: the condition thick, the 4 pairs thin, the 2
   relabelings dashed grey, with segment boundaries marked.
2. **Member heatmap trio.** Rows are members, at most 40, ordered by the reference group's peak bin. Columns are 20
   equal-count bins per specimen on the analysed position (72's panel-a rule; Σcounts/Σlibrary per bin). The three
   panels are the reference-group mean, the case-group mean (both centred per gene over both groups' bins) and their
   difference (centred by the median gene per bin). Each bin's majority label is marked.
3. **Segment bars.** The condition P_k as bars, the 4 pairs as dots and the 2 relabelings as crosses. The published
   segment is outlined. The title gives the replication and localization classes.

**Summary heatmap.** Programs × S1/S2/S3, coloured by Q_k in log2 (diverging, clipped at ±2). Each cell shows a dot
per agreeing pair. The localized segment is boxed, the published segment underlined, and the replication class marked
at the row end. 79 has one panel per family; 80 has one.

**Supplementary.** Secondary programs (the same summary heatmap); the equal-probe sensitivity (79); the c19 and
untrimmed arms (80); the coordinate vs label peak agreement.

## 5 · Decision rules (fixed now)

All thresholds below come from the record: 0.25 log2 and the 0.25 relabeled ratio are 76's O5 values, q 0.05 is the
project's set rule, and "all pairs" is 76's.

### 5.1 Replication (families S and I)

| Class | Rule |
|---|---|
| replicated | BH q ≤ 0.05 (calibrated p, section 6) **and** specific (relabeled one-sided p > 0.05 in both) **and** the condition P at the claim's scope has the claimed sign with \|P\| ≥ 0.25 **and** the claimed sign holds in all 4 pairs |
| direction only | the claimed sign in the condition and in all 4 pairs, but any other part fails |
| not replicated | none of the above |
| opposite | opposite-direction BH q ≤ 0.05 **and** the opposite sign in all 4 pairs |

The class always uses all program genes (O8). The labels "probe-sensitive" and "length-sensitive" (79) are added when
the class changes under the equal-probe or gene-length sensitivity of 4a. "Probe check not possible" is shown when
fewer than 5 equal-probe members remain.

### 5.2 Zonation: source species and transfer (family Z)

| Class | Rule (per species; "both specimens" = that species' two) |
|---|---|
| zonated | BH q ≤ 0.05 for the within-species test **and** the segment contrast ≥ 0.25 log2 in both specimens |
| attenuated | the claimed sign in both specimens, but not zonated |
| absent | the signs differ between the specimens, or both are within ±0.25 |
| reversed | the opposite sign with \|contrast\| ≥ 0.25 in both specimens |

- **The confirmatory claim of family Z is the positive control (decision O7):** "replicated in the source species"
  means *zonated* in the species the claim comes from. Failure flags the labels or the units.
- **The added-information readout is the cross-species transfer:** the class in the other species. Is the mouse S3
  program also S3 in human, and the human S3 program also S3 in mouse? It is "conserved" when zonated. It is read
  together with the program's species-difference localization (5.3): a conserved program should show a uniform or
  absent species difference, and a mouse-only program a difference localized to its own segment.
- The attenuation ratio is reported against the genome-wide ratio (4a).

### 5.3 Localization (all families)

| Class | Rule |
|---|---|
| localized to S_k | between-segment BH q ≤ 0.05 **and** both relabeled ratios < 0.25 **and** the argmax of Q is S_k in the condition and in all 4 pairs **and** the margin is ≥ 0.25 log2 |
| graded | between-segment BH q ≤ 0.05 and specific, but the peak is not unanimous across pairs or the margin is < 0.25. The ordering of Q is reported |
| uniform within resolution | between-segment q > 0.05 **and** the condition range of Q over S1–S3 is < 0.25 **and** Q has the same sign in all three segments and in all 4 pairs **and** \|Q\| ≥ 0.25 somewhere. This means a difference of the same size along the PT |
| no difference | every \|Q_k\| < 0.25 in the condition and q > 0.05 |
| not resolved | anything else, for example a range ≥ 0.25 without a significant between-segment test, or pairs that disagree |

Per-segment flags:
- "present in S_k": Q_k ≥ 0.25 in the condition, with the sign in all 4 pairs.
- "no difference resolved in S_k": |Q_k| < 0.25. Equality is never claimed.

"Agrees with the published segment" applies to programs with a published segment: the label-based peak lies in it
(S1 or S2 for S1–S2 programs).

For the coordinate version, "changes along the axis but not between segments" (z_pos significant, z_seg not) is
reported as within-segment and descriptive only (C14).

### 5.4 What counts as informative

- **Position adds information** for a program when its localization class is localized, graded or uniform within
  resolution. The published lists cannot make any of these three statements: where the difference is largest, that it
  changes between segments, or that it is the same size across segments.
- **Not informative:** "not replicated" or "direction only", together with "not resolved". Our data then neither
  reproduce nor place the published claim.
- **Overall readout:** per family, the counts per class, and how many located programs fall outside the published
  segment.
- The abandonment rule in the Why section applies to the confirmatory families.

## 6 · Multiplicity and calibration

- BH runs within each family × question, over the confirmatory programs only:
  - S: replication (6) and between-segment (6);
  - Z: source species (6; the confirmatory positive control), transfer (6) and between-segment (6). Transfer and
    between-segment are the family's added-information readouts, BH-adjusted the same way;
  - I: replication (9) and between-segment (9).
- The coordinate between-position test (z_pos) gets the same BH within its family and is labelled a sensitivity.
- Secondary programs carry nominal p, labelled secondary, and are never pooled into a family.
- Calibration is checked per statistic and scope with 75's swap rule on the library sets (section 4a). A failing
  statistic's program p becomes the empirical p against the relabelings (76's rule), before BH.

The ceiling of that rule: it treats relabeled sets as exchangeable nulls. They carry real between-specimen
differences (AKI's A2/A4 factor), so it errs conservative.

## 7 · Sensitivities (descriptive; a sign or class that changes is stated beside the primary)

**Both notebooks.**
- The untrimmed units (all PT units, no common-support trim).
- The program mean instead of the median (O4): the segment summaries, the classes of 5.3 and the curve summaries
  recomputed with means.
- 75's (or the refit's) Δ_S peak against the pseudobulk peak.
- The secondary programs.

**79 only.**
- The equal-probe and gene-length sensitivities of 4a. A changed class is labelled (5.1).
- The two human sections separately: the human specimen-pair columns are already in the 4 pairs. HUK1_MED1's crossing
  transitions (C1) are flagged on curves.

**80 only.**
- **The c19 arm.** This is 76's arm (ii) in pseudobulk only: 76's `placement_units.csv` places the 80 c19 units at
  their majority PT label (73 S3, 7 S2), and their counts come from the run folder. No model refit.
- c15 is not added (76's O4).

## 8 · Outputs (`tables/`)

- `programs.csv` and `programs_members.csv`, the frozen copies with members in and out of the universe.
- `pseudobulk_rates.parquet` (genes × specimen × scope).
- `gene_level.csv`: per program × member, the centred group and pair differences per scope, t, z_seg, z_pos, the
  flags (converged, separated, probe balance), and the published statistic.
- `calibration.csv`: the swap rule per statistic × scope.
- `program_tests.csv`: per program × question, the raw p, calibrated p, BH q, relabeled p and the conventional
  comparator.
- `segment_summary.csv`: program × partition or pair × scope, with median, mean, n and sign shares.
- `profiles.npz` and `profile_summaries.csv`: the curves.
- `classes.csv`: replication, transfer and localization classes, peak, margin, axis agreement and flags.
- `added_information.csv`: one line per program.
- `sensitivity_*.csv`.
- `run_record.json`.

Figures follow section 4d (PNG and PDF).

## 9 · Caveats that travel with every result

- The cohort is two male control mice and two cortex sections of one male human donor (79), or two control and two
  IR male mice (80). Specimens are the replicates, and everything is descriptive.
- Human sections are cortex only, so in 79 the S3 comparison is human cortical S3 against mouse S3 including the outer
  stripe. Human PT-S2 is 632 units.
- Probe panels differ between species: level claims in 79 are probe-exposed, and between-segment claims are not.
- The 80 cohort lacks the most injured units: c15 (2,425 units, 72% IR) is dropped, and c19 (80 injured PT units) is
  labelled cTAL. Injured-state programs are expected to be attenuated, most in S3. The IRI time point is unknown, so
  which Kirita states replicate is read only as "compatible with the states Kirita observed at …", never as a
  time-point claim.
- The published lists differ from our measurement: single-cell or single-nucleus, other comparators (one species
  against two, state against all PT), female mice (Dixon), human AKI of mixed causes (Hinze).

## 10 · New code: `pseudospace/program_localization.py` (plan)

Thin, reusing `aki_along_pt` (`centre_rows`, `program_profile`, `profile_summary`, `empirical_p`, `segment_units`),
`pathway_pipelines` (`camera_pr`, `squeeze_var`, `partition_groups`), `discrete_comparators` (`segment_f_test`,
`swap_rule`) and `species_calls.to_segment_units`.

| Function | Does |
|---|---|
| `read_programs(path, universe, exclude, sha256)` | reads a committed program file (`docs/results/known-programs/`), checks the hash, that members lie in the universe and that none is excluded; returns `{id: members}` and the metadata |
| `gene_span(fasta_gz, symbols)` | human genomic span per symbol from the Ensembl cDNA headers (4a) |
| `stratified_rank_normal(statistic, strata)` | `rank_normal_scores` within each stratum (gene-length deciles) |
| `segment_pseudobulk(counts, library, specimen, segment, scopes, floor=0.5)` | log2 rates, genes × (specimen, scope) |
| `centred_contrast(rates, weights, reference)` | a linear contrast over specimens, minus its median over the reference genes |
| `partition_effects(rates, case, reference, partition)` | the effect, s² and df of 3.3 for the condition or a relabeling, in closed form |
| `segment_contrast(rates, specimens, segment, segments)` | the within-species effect, s² and df of 3.4 |
| `moderated_t(effect, s2, df, covariate)` | squeeze_var, then t, df_total and the two-sided p |
| `program_test(statistic, programs, direction)` | cameraPR on the programs, with one-sided p in the claimed and the opposite direction |
| `calibrated(p, null_p)` | the empirical p (wraps `empirical_p`) |
| `segment_summary(contrasts, members, min_members=5)` | the long table of section 8 |
| `replication_class`, `transfer_class`, `localization_class` | the rules of section 5 as pure functions of the summary rows and q |
| `added_information(row)` | the 4c template |
| `bin_rates(counts, library, position, specimen, n_bins=20)` | the heatmap bins |

Tests (`tests/test_program_localization.py`, synthetic; `pytest.importorskip('rpy2')` for cameraPR):
1. Pseudobulk sums and the floor on a hand-made matrix.
2. The closed forms of 3.3 equal an OLS fit of the 4-sample design for the condition and both relabelings.
3. `moderated_t` with d₀ = 0 equals the ordinary t, and the moderation shrinks s² towards the trend.
4. `program_test` on a planted human-high set: a small one-sided p for the "human-high" claim and a large one for
   "mouse-high".
5. `localization_class`:
   - a planted S3-only difference is "localized to S3";
   - a constant difference is "uniform within resolution";
   - a pair with another argmax gives "graded";
   - zero difference gives "no difference";
   - a relabeled ratio ≥ 0.25 blocks "localized".
6. `replication_class` and `transfer_class` at each boundary, including "opposite" and "reversed".
7. `segment_summary` returns NaN below `min_members`.
8. `read_programs` refuses a wrong hash, a member outside the universe and an excluded gene.
9. `gene_span` on a two-gene synthetic gzip header file: the span is the outermost transcript bounds.
   `stratified_rank_normal` is rank-normal within each stratum and unchanged by a stratum-constant shift.

## 11 · Compute

| Notebook | Main costs | Cold run | Cached rerun |
|---|---|---|---|
| 79 | no NB fit; pseudobulks about 1 min; about 40 cameraPR calls on up to 1,513 sets, about 1–2 min; figures about 2 min | under 10 min | under 2 min |
| 80 | the segment F-test on U76 (3 partitions, 24 jobs, about 3–5 min, cached); the c19 pseudobulk about 1 min; the rest as in 79 | about 15 min | under 2 min |

Memory stays below 16 GB.

## 12 · Coordinator decisions (2026-10-12)

| # | Decision | Where |
|---|---|---|
| O1 | Accepted: Hinze 2022 (within-PT) is primary for human zonation; Lake 2023 is secondary | 1.4, 1.6 |
| O2 | The coordinator copies the frozen CSVs to `docs/results/known-programs/` and commits them with the protocol. The notebooks read them from there and check the sha256 (1.7). `build_programs.py` takes `--scratch`/`--repo` (no absolute path) and stays in scratch, because it reads the scratch lit-full extraction. The committed CSVs are the source of truth | header, 1.7, 2, 10 |
| O3 | Non-one-to-one genes are excluded and counted per program (`n_non_one_to_one`, listed) | R2, CSVs |
| O4 | The primary program curve is the median of the members' log-ratio curves; the mean is a pre-specified sensitivity | 3.5, 7 |
| O5 | No alias mapping; unmapped symbols are counted per program (`n_unmapped_symbol`, listed) | R1, CSVs |
| O6 | Accepted: the Zhou core rule (≥ 3 of 6 IRI groups, never opposite) | 1.5 |
| O7 | Family Z is confirmatory as a positive control (replication of published zonation in the same species). Its added-information readout is the cross-species transfer | 5.2, 6 |
| O8 | "Replicated" uses all program genes. The equal-probe restriction is a pre-specified sensitivity, and a program whose verdict changes under it is labelled "probe-sensitive" | 4a, 5.1 |
| CellSpectra caveat | Turned into a pre-specified gene-length sensitivity (human genomic span from the public Ensembl 110 cDNA headers; length-decile stratified rank-normal level t). It is cheap (a header parse) and covers the whole universe. A changed verdict is labelled "length-sensitive" | 1.3, 4a, 5.1 |

Nothing is open. Phase 2 starts on the coordinator's go, after the commit.
