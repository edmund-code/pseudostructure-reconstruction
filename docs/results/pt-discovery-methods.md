# Ways of finding biology along the PT: the running scoreboard

The user's request (2026-10-18): "keep analyzing different ways of finding biological information in our data that we
will eventually write into the paper … just keep trying stuff to decide on the best method."

This page keeps one row per way of finding biology that has been tried. Each row records what it is good for, how it
did on the shared yardsticks, and its verdict. It is updated after each experiment, and the evidence lives in the
linked results notes. Notebook numbers are internal; shareable pages use "reconstructed pseudospace" and no notebook
names.

## The yardsticks (in order of weight)

1. **External replication.** Does the call hold in independent data?
   - AKI: Kirita 2020's IRI time course, 18 mice (notebook 94), and human AKI biopsies from KPMP (notebook 98).
   - Human vs mouse: the cell atlases (Lake/KPMP, Census) and microdissection (Chen 2021).
2. **Internal specificity.** Calls under balanced specimen relabeling against the condition: NCDR and calibration
   share.
   - AKI caveat (notebook 95): the A2/A4 relabeling correlates with the IRI response itself (ρ 0.23–0.66). Relabeled
     "false" calls on AKI are partly shared, injury-like biology, so external replication outranks internal
     calibration for AKI.
3. **Replication between specimens:** pair agreement, and independent-pair shape correlation.
4. **Stability across reconstructions:** notebook 69's coordinate against the package's.
5. **Probe exposure (human vs mouse):** the share of calls that follow probe-panel differences rather than biology.
6. **Recovery of named literature expectations,** fixed before the fit.
7. **Readability:** how many programs the calls reduce to, and whether the reduction is coherent.

## Calling methods (which sets or genes differ)

| Method | Question | Specificity | External | Other | Verdict |
|---|---|---|---|---|---|
| **Positional F-test** (NB T_spatial, moderated F, cameraPR; 71) | Does the shape of the difference along the PT differ? | hvm pass (0.071); AKI fails calibration (0.120; 76: 0.144) | AKI shape calls replicate 6.8–8.0× the matched null; set ranking ρ 0.54 with IRI shape (94). Human AKI: partly (genes 1.4×, sets not; location unmeasurable in biopsies, 98). hvm atlas replication not set-specific (77). | Reconstruction Jaccard 0.84 / 0.72; probe exposure 0.17; hub robustness 0.92 / 0.99 (93) | **Primary for positional differences** |
| Segment F-test (S1/S2/S3 steps; 74, 75) | The same, with steps | hvm pass (76 calls, 71 shared); AKI fails | – | The spline adds 12 sets, the steps 5 | Discrete twin; continuous adds little set-level |
| **Clustering** (package gene groups and clusters, matched ORA; 71) | Which sets sit in up, down or shape gene groups? | passes both | AKI: replicates as **level** (1.5×), not as shape (0.5×) (94); sets do not generalise to human AKI (1.06×, 98) | Reconstruction Jaccard 0.42 / 0.66; probe exposure 0.69 | **Complementary level layer**; not for shape |
| **Level / location screen** (L, M, A; 91) | Whole-PT level change, location beyond level, gradient survival | relabeling bounds per call (level 40, 3 against 3,159 genes) | AKI level sets the most consistent level list: 0.57, 1.7×, ties clustering. Location sets 5.3× at 12 h (94). Best of our lists in human AKI (sets 1.55×, 98). | Filters show the male-biased and fasting-down programs fall in IRI too (94) | **Best AKI level description** |
| Position-split ORA (scoranet-style row N; 93) | Gene-first ORA, split up/down × segment | passes | AKI location 5.0× (94) | Reconstruction Jaccard 0.45 / 0.43; hub robustness 0.43 / 0.72 | Not better than the F-test (93) |
| DESeq2 + ORA, whole PT or per segment (75) | Level differences | passes | – | hvm probe exposure 0.60–0.71 | Not for species; level only |
| cameraPR or GSEA on raw T_spatial, ssGSEA, DE + signed GSEA (66, 75) | various | fail | – | – | Rejected |
| Regulon activity (CollecTRI, ulm; 95) | Which regulators' targets change? | AKI 0.62 relabeled per call, against member sets' 0.34 | at chance (1.17×) | Mouse zonation regulons rarely carry to human | **Worse than member sets**; names are hypotheses only |
| Co-expression conservation (Crow 2022's partner AUROC, level-free; 97) | Does a gene keep its co-expression partners across species? | within-species median 0.67 (bimodal) | rewiring replicates in the atlases (ρ 0.30; 0.34 for reproducible genes) | reversed-gradient genes most rewired; only xenobiotic metabolism rewired at set level | Not shown useful genome-wide; **supporting level-free cross-check** of species gene lists |
| Signalling footprints (PROGENy, ulm; 95) | Which signalling pathways act? | 3 AKI calls | EGFR, MAPK, PI3K replicate (z 5.3–7.8) | KEGG ErbB and MAPK miss them; footprints miss p53, hypoxia, JAK-STAT | **Add-on for signalling only** |
| Published programs re-read (79–82) | Do published findings hold along our PT? | per program | AKI programs replicate 2.0× at 12 h (94, Kirita's own excluded) and best of all lists in human AKI (13/22, 19.9×; 98) | – | Hypothesis-driven checks; keep |

## Organising and describing calls

| Method | What it gives | Verdict |
|---|---|---|
| Overlap programs (75) and redundancy clusters (71) | Groups of calls sharing genes | Current grouping; F-test calls not coherent beyond chance |
| Pathway networks (scoranet's Jaccard + Leiden, with positional additions; 93) | 1.1–3.9× fewer programs, coherent beyond chance | Settings-sensitive; **descriptive layer**, shown with its sweep |
| Gene shape classes against segment steps (90 A) | Continuous favoured for 28% of non-flat genes against 3% for steps | Supports continuous over discrete at gene level |
| Positional information decoders (90 C) | 26 single pathways localise a unit better than its label | Novel methods number; keep |
| Change density along the axis (90 E) | Change concentrates beside, not at, the label transitions | Descriptive; keep |
| Level-free species comparison (90 D, 91 S) | 74 reversed gradients; about a quarter of the gradient shared; no pathway-level species program | Keep; species differences are gene-level |
| Positional programs against a cluster view (89) | Broad opposite gradients, no boundary programs | Keep for the continuous claim |
| Phase decomposition with an external time course (Kirita as a ruler; 99) | Acute/sustained present, repair partial, late absent: our IR ≈ 12 h–2 d; acute injury S3-weighted (H1) beyond the controls' spread | **New and useful**: dates a single-time-point dataset and locates each phase |
| Tissue location (depth, neighbours; 96) | Nephron position organises injury; depth adds nothing; injury graded per tubule | Supports the nephron coordinate as sufficient |

## Current recommendation (2026-10-18)

- **Positional (shape) differences: the positional F-test.**
  - It is the only calibrated test on human vs mouse.
  - Its AKI shape calls replicate in an independent IRI experiment.
  - Its ranking agrees with external shape at ρ 0.54.
- **Level differences (AKI): notebook 91's level screen, or clustering,** judged by external replication.
  - Never read species levels: the probe panels differ.
- **Interpretation:**
  - group calls with the overlap programs or the networks, shown with their settings sweep;
  - add PROGENy for signalling;
  - treat regulator names as hypotheses.
- **Validation:** for AKI, external replication (Kirita) is the arbiter, since the internal relabeling contains
  injury-like biology.

## Still to try (queue)

1. **Spatial repair niches.** Cell-cycle and transport-loss residuals cluster among neighbouring tubules beyond
   position and depth (96, secondary). A frozen test of niche structure, and of what neighbours share, is the follow-up.
