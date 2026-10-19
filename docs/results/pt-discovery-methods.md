# Ways of finding biology along the PT: the running scoreboard

The user's request (2026-10-18): "keep analyzing different ways of finding biological information in our data that we
will eventually write into the paper … just keep trying stuff to decide on the best method."

This page keeps one row per way of finding biology that has been tried. Each row records what it is good for, how it
did on the shared yardsticks, and its verdict. It is updated after each experiment, and the evidence lives in the
linked results notes. Notebook numbers are internal; shareable pages use "reconstructed pseudospace" and no notebook
names.

## The yardsticks (in order of weight)

0. **Planted truth** (notebooks 100, 102, 103): power and false calls on our own data with known pathway-level differences. This is the only yardstick with a ground truth; it shows what each method can see, not whether a real call is right.

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
8. **Contamination exposure** (notebook 104): the share of a method's driver genes that are edge-enriched (neighbouring
   structure RNA), against expression-matched genes.

## Calling methods (which sets or genes differ)

| Method | Question | Specificity | External | Other | Verdict |
|---|---|---|---|---|---|
| **Positional F-test** (NB T_spatial, moderated F, cameraPR; 71) | Does the shape of the difference along the PT differ? Planted: finds shape, S3-only and S2-bump changes at 0.5 log2 (≥ 0.96), blind to pure level (≤ 0.06); 5–12 BH-driven clean false sets on the AKI base (which had 29 marginal calls), 0.6–1.9 on human vs mouse; profile generalises across datasets (100, 103). | hvm pass (0.071); AKI fails calibration (0.120; 76: 0.144) | AKI shape calls replicate 6.8–8.0× the matched null; set ranking ρ 0.54 with IRI shape (94). Human AKI: partly (genes 1.4×, sets not; location unmeasurable in biopsies, 98). hvm atlas replication not set-specific (77). | Reconstruction Jaccard 0.84 / 0.72; probe exposure 0.17; hub robustness 0.92 / 0.99 (93); contamination exposure 1.44× hvm (least), 1.18× AKI (n.s.) (104) | **Primary for positional differences** |
| Segment F-test (S1/S2/S3 steps; 74, 75) | The same, with steps. Planted: matches the positional F-test, S2 bump included (100). | hvm pass (76 calls, 71 shared); AKI fails | – | The spline adds 12 sets, the steps 5 | Discrete twin; continuous adds little set-level |
| **Clustering** (package gene groups and clusters, matched ORA; 71) | Which sets sit in up, down or shape gene groups? Planted: finds level at 0.5 log2 (0.92), shape only at 1 log2; 0–5 clean false sets (100). | passes both | AKI: replicates as **level** (1.5×), not as shape (0.5×) (94); sets do not generalise to human AKI (1.06×, 98) | Reconstruction Jaccard 0.42 / 0.66; probe exposure 0.69; contamination exposure 1.69× hvm, 1.04× AKI (104) | **Complementary level layer**; not for shape |
| **Level / location screen** (L, M, A; 91) | Whole-PT level change, location beyond level, gradient survival | relabeling bounds per call (level 40, 3 against 3,159 genes) | AKI level sets the most consistent level list: 0.57, 1.7×, ties clustering. Location sets 5.3× at 12 h (94). Best of our lists in human AKI (sets 1.55×, 98). | Filters show the male-biased and fasting-down programs fall in IRI too (94) | **Best AKI level description** |
| Position-split ORA (scoranet-style row N; 93) | Gene-first ORA, split up/down × segment | passes | AKI location 5.0× (94) | Reconstruction Jaccard 0.45 / 0.43; hub robustness 0.43 / 0.72; contamination exposure 1.89× hvm (104) | Not better than the F-test (93) |
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
| Groupings scored against a planted truth (RED, OVL, NET; 102) | Planted programs spill into 200–690 called sets; NET returns 8–10 programs with ≈ the baseline's false programs and keeps 1–3 planted programs separate (0.88–0.96), but merges 15 strong programs (0.5); OVL separates always but leaves ~6 false programs per planted one; RED removes only ~40% of spill | **NET is the best grouping at realistic doses; OVL beside it for very strong contrasts; never count raw sets** |
| Gene shape classes against segment steps (90 A) | Continuous favoured for 28% of non-flat genes against 3% for steps | Supports continuous over discrete at gene level |
| Positional information decoders (90 C) | 26 single pathways localise a unit better than its label | Novel methods number; keep |
| Change density along the axis (90 E) | Change concentrates beside, not at, the label transitions | Descriptive; keep |
| Level-free species comparison (90 D, 91 S) | 74 reversed gradients; about a quarter of the gradient shared; no pathway-level species program | Keep; species differences are gene-level |
| Positional programs against a cluster view (89) | Broad opposite gradients, no boundary programs | Keep for the continuous claim |
| Phase decomposition with an external time course (Kirita as a ruler; 99) | Acute/sustained present, repair partial, late absent: our IR ≈ 12 h–2 d; acute injury S3-weighted (H1) beyond the controls' spread | **New and useful**: dates a single-time-point dataset and locates each phase |
| Family-wide transporter handoff screen (label-permutation null, microdissection and atlas confirmation; 101) | Handoffs are not a family-wide principle (Q1 fails); 11 conserved, confirmed handoffs (amino acids, organic anions, SGLT2→SGLT5), overlapping and inside S2 | **Useful**: separates a general claim from specific, confirmed micro-gradients |
| Tissue location (depth, neighbours; 96) | Nephron position organises injury; depth adds nothing; injury graded per tubule | Supports the nephron coordinate as sufficient |
| Bin-level contamination audit (core rebuild against size-matched random bins; 104) | Neighbour RNA is 2–4× higher at unit edges and follows the anatomy along the PT (AL toward S3, DCT and glomerulus in early S1, infiltrate in late IR S3), but 96.5–98% of gradients and 98.5–99.7% of AKI calls survive the core rebuild; human vs mouse calls are 1.4–1.9× enriched for edge genes (F-test least), AKI calls not | **Keep as a standard check**; flags the few dozen neighbour-RNA genes |

## Current recommendation (2026-10-18)

- **Positional (shape) differences: the positional F-test.**
  - It is the only calibrated test on human vs mouse.
  - Its AKI shape calls replicate in an independent IRI experiment.
  - Its ranking agrees with external shape at ρ 0.54.
  - With a planted truth it finds shape, S3-only and S2-bump changes at 0.5 log2 in both datasets (100, 103).
- **Level differences (AKI): notebook 91's level screen, or clustering,** judged by external replication.
  - The positional F-test is blind to pure level by construction, and clustering finds planted level shifts (100).
  - Never read species levels: the probe panels differ.
- **Interpretation:**
  - read calls as driver-gene network programs (NET), with the overlap programs (OVL) beside them for very strong
    contrasts (102). Never count raw sets: one planted program makes 100–650 sets significant;
  - add PROGENy for signalling;
  - treat regulator names as hypotheses.
- **Validation:** for AKI, external replication (Kirita) is the arbiter, since the internal relabeling contains
  injury-like biology.

## Still to try (queue)

1. **Spatial repair niches.** Cell-cycle and transport-loss residuals cluster among neighbouring tubules beyond
   position and depth (96, secondary). A frozen test of niche structure, and of what neighbours share, is the follow-up.
