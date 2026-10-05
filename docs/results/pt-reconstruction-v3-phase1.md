# Workstream A (ws7), phase 1: a reconstruction method for PT cross-sections

Owner: ws7 (notebooks 51–55, `results/pt_reconstruction_v3/`). Phase 1 only: plan plus cheap,
descriptive feasibility checks. No coordinate was fitted and no benchmark metric was computed.

Sections:
1. Failure analysis
2. Literature
3. Problem structure
4. Feasibility protocol and results
5. Candidate methods
6. Pre-registered phase-2 evaluation
7. Recommendation
8. Side findings for the coordinator
9. Files

**Provenance.** The feasibility protocol in section 4 was written to this file before any check ran.
Its two addenda were each written after the preceding check and before the next one ran. Sections
1–3 and 5–9 were written afterwards.

---

## 1. Failure analysis: why scFates and DPT miss the within-segment and held-out-gene bars

This section uses existing tables (nb35/36, nb42, `phase1_diagnostics`) plus the feasibility checks
of section 4 (F1–F6).

### 1.1 What fails (existing tables)

| Bar | SCF13 | DPT13 | Segment labels / other | Source |
|---|---|---|---|---|
| External concordance P1, fold-protected (mouse male snRNA / human half B) | 0.74 / 0.25 | 0.74 / 0.29 | labels 0.73 / 0.27 | nb35 `p1_external_concordance.csv` |
| Within-segment concordance P1w (mouse / human) | 0.17 / 0.10 | 0.20 / 0.14 | CAC 0.29 / 0.16 | nb36 `summary.json` |
| Gene-fold within-segment agreement (S1 / S2 / S3) | 0.44 / 0.34 / 0.70 | 0.58 / 0.66 / 0.53 | CAC 0.09 / −0.11 / −0.47 | nb36 |
| Held-out genes P4: gain over segment means (mouse / human) | +0.09% / +0.02% | −0.43% / −0.01% | CAC +0.91% / +0.75% | nb36 `p4_heldout_prediction.csv` |
| P4: gain over the segment + coverage quadratic oracle | −1.2% / −1.3% | −1.8% / −1.3% | CAC −0.4% / −0.6% | same |
| Equal-depth agreement P5 (all; human S1) | 0.72; 0.36–0.46 | 0.80; 0.71–0.72 | – | nb35 `p5_exposure_invariance.csv` |
| Species registration gap across refits | 0.29 (S1<S2<S3 breaks in one refit) | 0.14 | CAC 0.52 | nb35/36 |
| Count split: coordinate from half A / half B vs saved | 0.68 / 0.49 | 0.82 / 0.82 | – | nb35 |
| Numerical floor: 10⁻¹⁶ input change | moves structures by up to 0.44 | median within-segment ρ 1.00 | – | nb35 |

Every candidate coordinate predicts held-out genes worse than a segment + coverage oracle (log
library, log spots, log detected genes). Coverage covariates cut held-out-gene error by about 1.3%
over segment means. That is about 15× SCF13's gain and 1.4× CAC's. P4 as defined is dominated by
how well the exposure dependence of log-normalized values is modelled; it says little about
position.

### 1.2 Measurement limits (F1, F2, F5, F5b, F6; section 4)

**Counts.** Mouse PT structures carry 5.5–9.4k UMI over 340–400 bins (15–22 UMI per 2 µm bin). Human
structures carry 2.3–7.8k UMI over 450–970 bins (4.5–7.7 UMI per bin). Human S3 is shallowest (2.3–2.6k).

**Depth is mostly cut size.** Polygon short side (tubule width in the section) correlates with log
library at ρ = 0.60–0.83 in every specimen × segment. "Depth dependence" is therefore largely
dependence on physical cut size, which varies by segment: human S1 short side is 60 µm, against 41 µm
in S3.

**The ceiling for a fixed positional score.** S = log(late/early landmark counts), with landmarks
taken from selection references only. With genes and reads both disjoint (F5b), within-segment
agreement is:

| | S1 | S2 | S3 |
|---|---|---|---|
| Mouse | 0.69–0.70 | 0.74–0.83 | 0.66–0.68 |
| Human | 0.20–0.22 | 0.47–0.48 | 0.35–0.41 |

So **in mouse, a reproducible within-segment axis exists in the counts, shared across independent
gene halves. In human S1 the measurement is the bottleneck.**

**Label noise.** Human label-transfer agreement is 0.71 (κ 0.54), against 0.87–0.92 in mouse (nb42).
Human external references resolve only convoluted PT versus S3, so within-early order has no human
external standard.

### 1.3 Method limits

- **The generic coordinates miss the reliable mouse axis in S1 and S3.** Within-segment Spearman of
  the coordinate with a landmark half-score (F5b):
  - SCF13: S1 0.13–0.33, S2 0.56–0.64, S3 −0.12 to 0.22;
  - DPT13: S1 0.05–0.20, S2 0.61–0.73, S3 0.50–0.56.

  The landmark score is reproducible at 0.66–0.83 across disjoint genes and reads. Within S1 (both
  coordinates) and S3 (scFates), the generic coordinates therefore order structures along something
  else. This is a method limit, not a count limit.
- **The representation leaks exposure and cut size.** The curve is fitted on Harmony PCs of
  log-normalized counts, and log1p(c/L·10⁴) depends on L for low counts (Townes et al. 2019).
  - Within segment, position correlates with log library at −0.31/−0.38 (SCF13, human S1) and
    +0.25 (SCF13, mouse S2).
  - The pseudostructure line saw within-segment order track detected-gene counts at ±0.94 in a
    log-normalized latent model; count-exposure models removed most of it (ledger, notebooks 16–17).
- **Projection instability.** The 30-node elastic curve is projected in 5 dimensions and can fold.
  Changes of order 10⁻¹⁶ move structures by up to 0.44, and gene thirds disagree within segments.
  The order is decided by small residual variation near nodes.
- **No explicit registration.** Species alignment is implicit, through Harmony and the joint curve,
  so it moves with the gene subset (gap 0.29). The gauge (curve arc length in Harmony space) is not
  physical, which makes peak positions coordinate-dependent.
- **Small anchor panels were CAC's problem, not anchors as such.** Notebook 36's 160 conserved anchors
  give a within-segment count-ratio half agreement of only 0.29–0.45 in mouse and 0.05–0.35 in human
  (F5c; median 0.31). That is comparable to CAC's 0.38 through Harmony and scFates. Larger,
  species-specific selection panels (515 mouse, 677 human genes) reach 0.66–0.83 in mouse (F5b).
  CAC's P4 gain (+0.9%, ten times SCF13's) and its P1w (best of all arms) show that the landmark
  information is useful. Its stability failed because the panel was too small.

**Diagnosis.**
- *Mouse:* method-limited. The counts carry a within-segment landmark axis, but generic,
  similarity-only trajectories on log-normalized Harmony space do not use it in S1 and S3, and they
  absorb cut size.
- *Human S1:* measurement-limited. Positional information from expression alone is weak; any gain must
  come from pooling information, for example across physical neighbours.
- *P4:* insensitive by construction. It needs an exposure-aware variant (section 6), or no coordinate
  can be seen to help.

---

## 2. Literature and first principles

Sources marked † were found and read (abstract, or full-text excerpts where noted) in this phase
through Europe PMC; WebSearch was capped for the session. Sources marked ‡ are listed in the repo
survey (`docs/research/pseudostructure-methods-survey.md`) and were not re-read here.

| Source | What structure it exploits | Lesson for us |
|---|---|---|
| † Halpern et al. 2017, *Nature* ([PMID 28166538](https://pubmed.ncbi.nlm.nih.gov/28166538/); full-text excerpts) | Landmark genes with smFISH-measured zonation. Per-cell posterior over lobule layers from the UMI likelihood, which accounts for each cell's total UMI. Prior ∝ layer area in a hexagonal geometry | The closest precedent: an external positional standard, an exposure-aware count likelihood and a geometric occupancy prior. It replaces similarity-based ordering |
| † Halpern et al. 2018, *Nat Biotechnol* ([PMID 30222169](https://pubmed.ncbi.nlm.nih.gov/30222169/)) | Paired-cell sequencing: the attached partner's expression places the pair | Physical partners carry position. This motivates neighbour coupling |
| † Pham et al. 2023, stLearn PSTS, *Nat Commun* ([PMID 38007580](https://pubmed.ncbi.nlm.nih.gov/38007580/); full-text excerpts) | Pseudo-time-space distance mixes expression and physical distance with weight ω. Validated partly by spatial variograms | A precedent for spatial pseudotime. Spatial continuity cannot validate a spatially regularized method (circular) |
| † Ren et al. 2022, SpaceFlow, *Nat Commun* ([PMID 35835774](https://pubmed.ncbi.nlm.nih.gov/35835774/)) | Spatially regularized graph embedding; pseudo-spatiotemporal map | Same, with a learned embedding |
| † Shang & Zhou 2022, SpatialPCA, *Nat Commun* ([PMID 36418351](https://pubmed.ncbi.nlm.nih.gov/36418351/)) | Spatial-kernel prior on latent factors; used for trajectories | Spatial priors on latent variables are established; ours acts on a 1D position |
| † Zhao et al. 2021, BayesSpace, *Nat Biotechnol* ([PMID 34083791](https://pubmed.ncbi.nlm.nih.gov/34083791/)) | Bayesian use of spatial neighbourhoods for clustering and resolution | A neighbourhood prior on discrete labels; ours is on ordered position |
| † Chitra et al. 2025, GASTON, *Nat Methods* ([PMID 39849132](https://pubmed.ncbi.nlm.nih.gov/39849132/)) | An "isodepth" coordinate with piecewise-linear expression: continuous gradients plus discontinuities | Matches our "continuous within discrete" goal, but its coordinate is tissue space, not a repeated object |
| † Ma et al. 2022, Belayer, *Cell Syst* ([PMID 36265465](https://pubmed.ncbi.nlm.nih.gov/36265465/)) | Piecewise-linear expression in relative layer depth, discontinuous at boundaries | Same idea in layered tissue |
| † Kueckelhaus et al. 2024, SPATA2 spatial gradient screening, *Nat Commun* ([PMID 39179527](https://pubmed.ncbi.nlm.nih.gov/39179527/)) | Gradients along histology-defined axes | A supervised physical axis; ours is unknown |
| † Zeira et al. 2022, PASTE, *Nat Methods* ([PMID 35577957](https://pubmed.ncbi.nlm.nih.gov/35577957/)); † Tang et al. 2024, CAST, *Nat Methods* ([PMID 39294367](https://pubmed.ncbi.nlm.nih.gov/39294367/)); ‡ PASTE2, moscot | Cross-slice alignment by expression plus physical distance (fused Gromov–Wasserstein OT; GNN search-and-match) | They align tissue slices. Our specimens are different kidneys, with no shared geometry to align, so these do not transfer except as a partial-overlap caution |
| † Campbell & Yau 2019, Ouija, *Bioinformatics* ([PMID 29939207](https://pubmed.ncbi.nlm.nih.gov/29939207/)) | Pseudotime from a few marker genes with switch-like or transient priors | A precedent for landmark-gene pseudotime |
| † Townes et al. 2019, GLM-PCA, *Genome Biol* ([PMID 31870412](https://pubmed.ncbi.nlm.nih.gov/31870412/)) | Multinomial count model; log-normalization creates false variability | Supports a count likelihood with an exposure offset |
| † Neufeld et al. 2023, count splitting, *Biostatistics* ([PMID 36511385](https://pubmed.ncbi.nlm.nih.gov/36511385/)) | Poisson splitting separates estimating the latent variable from inference | Used to calibrate temperature and coupling on training data only, and as a stability metric |
| † Nyengaard 1999, *JASN* ([PMID 10232698](https://pubmed.ncbi.nlm.nih.gov/10232698/)) | Stereology from section planes | Profile counts on sections relate to 3D length. This is the basis of an occupancy gauge (C3). The formula is to be checked in full text before use |
| † Zhai et al. 2003, *JASN* ([PMID 12595496](https://pubmed.ncbi.nlm.nih.gov/12595496/)) | 3D reconstruction of 160 mouse PTs: each occupies a separate cortical domain; superficial PTs have long straight parts converging in medullary rays; deeper ones are longer and more tortuous; no ultrastructural segment boundaries in that strain | Adjacent convoluted cross-sections often belong to the same nephron, which favours neighbour coupling. Straight parts run in parallel bundles. Mouse "segments" are a continuum, which supports a continuous coordinate |
| † Ransick et al. 2019, *Dev Cell* ([PMID 31689386](https://pubmed.ncbi.nlm.nih.gov/31689386/)) | Nephron organization and composition depend on the time of nephron specification | **Confound:** between-nephron (u_r) variation is spatially structured, so spatial coupling can import nephron type as "position" |
| † Polesel et al. 2022, *Nat Commun* ([PMID 36175561](https://pubmed.ncbi.nlm.nih.gov/36175561/); full-text excerpts) | Functional sub-segments along the PT; S1 segments seen leaving the glomerulus | Axial heterogeneity within segments is functional. The glomerulus is a physical start anchor |
| † Chen et al. 2021, *JASN* ([PMID 33769951](https://pubmed.ncbi.nlm.nih.gov/33769951/)); † Lake et al. 2023, *Nature* ([PMID 37468583](https://pubmed.ncbi.nlm.nih.gov/37468583/)) | Segment-resolved mouse and human references | The external positional standard is segment-level only, so the gauge inside segments is not externally pinned |
| ‡ DPT, scFates/ElPiGraph, Slingshot, PhenoPath, MEFISTO, MATCHER, cellAlign, Genes2Genes, novoSpaRc, Lause et al. 2021, Saelens et al. 2019 | Similarity-based curves; covariates; alignment; transport | Comparators and cautions: MATCHER's uniform occupancy is not anatomy; forced alignment collapses |

**First principles.**
1. A cross-section is one point (or a short interval) on one nephron's path. Expression similarity
   orders by the dominant variance, which within a segment is cut size, exposure, composition and
   nephron type, not position.
2. Position is identified only through genes whose external, segment-level behaviour is known
   (landmarks), or through physical constraints (adjacency, anchors, stereology).
3. Each landmark contributes evidence ∝ its counts. A likelihood that conditions on exposure uses this
   evidence without importing cut size.
4. Neighbours add independent reads about a shared or nearby position, but also about shared nephron
   type and tissue fields. Coupling helps only where neighbour position agreement exceeds the confound,
   and that must be tested with placebo graphs and injected spatial artefacts.

---

## 3. Problem structure that generic methods ignore

| Structure | Evidence (section 4) | Covered by closed hypotheses? | Use |
|---|---|---|---|
| **(a1) Physical neighbours** | 90–96% of PT polygons have another PT polygon ≤ 5 µm away; same-label share 0.73–0.74 (mouse) and 0.62–0.63 (human) against 0.33–0.41 by composition. Neighbour carries gene-, read- and polygon-disjoint landmark signal: 0.17–0.47 in mouse, 0.11–0.25 in human (random partner ≈ 0). The library-adjusted result is the same | **No.** Notebooks 15–29 excluded x/y by design | C2. Needs coordinator/user approval: it reverses the line's "x/y never builds the coordinate" choice, and then no spatial quantity may validate |
| (a2) Glomerulus anchors | 168–217 (mouse) and 127–202 (human) glomeruli; PT within 10 µm: 543 mouse, 237 human; S1 odds ratio 4.5 / 2.0 (nb42). Attached PT are only 59–71% S1 | No | Weak and noisy as a position-0 anchor; descriptive only. Depth never validates |
| (a3) Outer stripe / medullary-ray geometry | Mouse S3 polygons are more elongated (AR ≥ 2 in 41–46% of S3, against 24–29% of S1/S2), consistent with in-plane pars recta. Human sections are cortex only | No | Orientation weight for the occupancy gauge (C3). An along-chain direction check needs user approval |
| (b) Within-polygon bins | Bins reproduce `n_spots` exactly (956/956). Major-axis landmark gradients: mouse 22% of elongated polygons at \|z\| > 1.96, against 12–14% on the minor axis and 6–9% for placebo genes. This is 1.6–1.8× the controls, below the pre-registered 2×. Human shows none | No | **Rejected** at the pre-registered bar. A heavy tail in mouse (mean z² 4.3 against 1.8 and 1.2) may mark polygons spanning a boundary or merged tubules; out of scope |
| (c) Repeated-object prior | One shared rate path per species, specimen-balanced | Partly: equal-specimen marginal (16), hierarchy (25–26), initialization (28) and ordered registration (29) are closed. Occupancy/stereology is not tested | C1 uses a shared path without offsets or hierarchy. C3 uses occupancy only as a gauge, never as a constraint |
| (d) Count measurement model | Cut size drives depth (ρ 0.60–0.83). A landmark ratio is exposure-free in expectation | Partly: notebooks 16–17 (count atlas helps thinning; NB projection lost fine programs); notebook 30's measurement-invariance protocol never ran | C1: NB with library offset on **landmarks**. New relative to 17, which used marker-excluded genes and a frozen projection |
| (e) Multi-specimen and species joint inference | – | Yes: offsets, registration and hierarchy (24–26, 29) are closed | Not rerun. Species get separate rate paths and a gauge shared through external segment structure |

---

## 4. Feasibility checks: protocol (frozen before running) and results

### 4.1 Protocol (verbatim, written before any check)

Descriptive only. No metric from nb35/36 (P1, P1w, P4, P5, gene-fold agreement, registration gap)
is computed, and nothing is tuned. External evaluation references (male mouse snRNA, human snRNA
half B) are not read. Positional scores use **selection references only**:
- mouse anchors: |S3 − early| ≥ 1 with the same sign in mouse microdissection and female mouse snRNA;
- human anchors: |S3 − early| ≥ 1 in human snRNA half A.
Both also require detection ≥ 10% in our PT structures of that species. Script:
`workstream-A scratch: feasibility.py`; outputs in `results/pt_reconstruction_v3/ws7_phase1_feasibility/`.

- **F1 · Polygon geometry.** Inputs: all 12,866 PT polygons, joined through `feature_index` with a
  ≤ 2 px centroid guard. Per polygon: area, minimum rotated rectangle long and short side (µm), aspect
  ratio (AR). Reported per specimen × segment.
- **F2 · Bins and counts.** `n_spots`, total counts, counts per bin. Guard: reassigning the 2 µm bins to
  the sampled polygons must reproduce the stored `n_spots` exactly.
- **F3 · Within-polygon gradient.**
  - *Sample.* Per specimen, up to 150 polygons per segment label with AR ≥ 2.5, drawn with seed 0.
  - *Bins.* Bins inside each polygon, eroded by 4 µm to limit edge spillover. Positions are projected
    on the major and minor axes of the minimum rotated rectangle.
  - *Statistic.* U = Σ_b L_b·s_b − (Σ L/Σ n)·Σ_b n_b·s_b, where L_b and E_b are the late- and early-anchor
    counts in bin b and n_b = L_b + E_b. It is compared with 200 within-polygon permutations of bin
    positions to give z.
  - *Controls.* (i) The same score along the minor axis. (ii) A placebo score from two random sets of
    non-zonated genes with matched expression (|contrast| < 0.2 in the selection references), along
    the major axis.
  - *Decision.* Within-object gradients count as *detectable* if, in at least one specimen per
    species, the anchor/major fraction with |z| > 1.96 exceeds both control fractions by ≥ 2× and by
    ≥ 5 percentage points. Otherwise the within-object channel is recorded as infeasible at this depth.
  - *Secondary, descriptive.* Is the excess larger in polygons near a reviewed segment transition
    (SCF13 within 0.08 of the specimen's S1/S2 or S2/S3 label boundary) than in polygons mid-segment?
- **F4 · Spatial neighbours.** Boundary-to-boundary distance from every PT polygon to the nearest other
  PT polygon. Reported:
  - the fraction with a PT neighbour within 5 µm and within 20 µm;
  - the fraction of such neighbours with the same label, against the label-composition expectation;
  - candidate same-tubule pairs: gap ≤ 3 µm, both AR ≥ 2, major axes within 20°, and the
    centroid-to-centroid vector within 20° of both axes.
- **F5 · Information carried by neighbours, against the measurement ceiling.**
  - *Score.* Per polygon S = log((L + 0.5)/(E + 0.5)) from polygon-level anchor counts, residualized
    on specimen × segment.
  - *Split.* Binomial thinning (p = 0.5, seed 0) gives independent read halves A and B per polygon.
  - *Within-segment correlations (Spearman):*
    - r_self = corr(S_A(i), S_B(i)), the split-half reliability (measurement ceiling of a fixed
      positional score);
    - r_nb = corr(S_A(i), S_B(j)), with j the nearest same-segment PT neighbour within 20 µm;
    - r_rand, the same with a random same-segment polygon (null);
    - r_lib, the neighbour correlation of residual log library size (technical spatial
      autocorrelation);
    - the Spearman of S with log library within segment (depth leakage of a ratio score).
- **F6 · Morphology, descriptive.** Within segment, the Spearman of polygon short side (tubule width)
  and area with the residual S, per specimen.
- **F7 · Anchors.** Glomeruli and glomerulus-attached PT per specimen (read from notebook 42 outputs;
  no recomputation).

None of these checks selects a method parameter. They decide only which structural channels are worth
a phase-2 protocol.

**Addendum (written after F1–F7 had run, before F5b ran).** F5's read-split reliability cannot
separate positional signal shared across genes from signal specific to single genes; notebook 36's
conserved-anchor coordinate had high P1w but negative gene-fold agreement. F5b therefore splits the
anchor genes into two halves by sha256 hash and also splits the reads, so the halves share neither
genes nor reads:
- r_gene = corr(S_{genes A, reads A}(i), S_{genes B, reads B}(i));
- r_nb_gene, the same between i and its nearest same-segment neighbour j.

Everything else is unchanged. F5b is descriptive and selects nothing.

**Addendum 2 (after F5b, before F5c ran).** F5c repeats F5b's gene- and read-disjoint agreement
using notebook 36's own 160 conserved anchors, signed by mouse microdissection. Through the generic
pipeline (Harmony, then scFates), the anchor halves agreed at a median of 0.38. F5c asks whether the
same genes carry a more reproducible within-segment signal when read as a count ratio. Descriptive only.

### 4.2 Results

All outputs are in `results/pt_reconstruction_v3/ws7_phase1_feasibility/`. Landmark panels: mouse
178 late and 337 early genes; human 147 late and 530 early (`gene_counts.json`).

**F1/F2 · Geometry and counts** (`f1_f2_summary.csv`; medians per specimen × segment).

| | Area µm² | Short side µm | AR | AR ≥ 2 | AR ≥ 3 | Bins | UMI | UMI/bin | Landmark UMI |
|---|---|---|---|---|---|---|---|---|---|
| Mouse S1/S2 | 1,370–1,605 | 35–37 | 1.6 | 24–29% | 3–6% | 341–403 | 5.5–9.4k | 14.6–21.9 | 979–1,314 |
| Mouse S3 | 1,473–1,475 | 33 | 1.8–1.9 | 41–46% | 13–15% | 368–369 | 6.3–7.5k | 17–19 | 615–739 |
| Human S1 | 3,671–3,896 | 58–60 | 1.5 | 12–13% | 1% | 912–968 | 6.4–7.8k | 7.2–7.7 | 428–552 |
| Human S2 | 2,970–3,080 | 52–54 | 1.4–1.5 | 12–14% | 1% | 731–766 | 4.2–5.2k | 5.7–6.6 | 247–325 |
| Human S3 | 1,810–2,123 | 41–47 | 1.2–1.4 | 5–11% | 1% | 452–532 | 2.3–2.6k | 4.5–5.8 | 131–163 |

The bin guard passed: 956/956 sampled polygons reproduce `n_spots` exactly. Human areas use the
Space Ranger scale (section 8).

**F3 · Within-polygon gradient** (`f3_summary.csv`): fraction of polygons with |z| > 1.96.

| Specimen | Polygons (AR ≥ 2.5) | Landmark, major axis | Landmark, minor axis | Placebo, major axis |
|---|---|---|---|---|
| Ctrl1A2 | 368 | 0.217 | 0.122 | 0.060 |
| Ctrl1A4 | 372 | 0.223 | 0.137 | 0.091 |
| HUK1_COR1 | 93 | 0.118 | 0.140 | 0.086 |
| HUK1_MED1 | 119 | 0.050 | 0.084 | 0.092 |

**Decision: not detectable under the pre-registered rule.** The mouse excess is +8.6 to +9.5 points
but only 1.6–1.8× the minor-axis control, below the 2× bar; human shows no excess. The minor-axis
control itself sits above 5%, which points to non-axial within-polygon structure (spillover or
oblique-cut asymmetry). The transition secondary is inconsistent: near versus away from a transition,
0.33 versus 0.18 in Ctrl1A4 and 0.21 versus 0.22 in Ctrl1A2. The within-object channel is dropped.

**F4 · Neighbours** (`f4_neighbours.csv`).

| | Ctrl1A2 | Ctrl1A4 | HUK1_COR1 | HUK1_MED1 |
|---|---|---|---|---|
| Median gap to the nearest PT polygon (µm) | 0.14 | 0.12 | 1.6 | 2.1 |
| PT neighbour ≤ 5 µm / ≤ 20 µm | 95% / 98% | 96% / 99% | 91% / 98% | 90% / 98% |
| Same label among neighbours ≤ 5 µm (composition expectation) | 0.74 (0.33) | 0.73 (0.34) | 0.62 (0.41) | 0.63 (0.41) |
| Candidate same-tubule (collinear abutting) pairs | 67 | 58 | 5 | 4 |

**F5/F5b/F5c · Information content** (`f5_*.csv`, `f5b_gene_split.csv`, `f5c_cac_anchor_ratio.csv`;
Spearman within specimen × segment).

| | Mouse S1 | Mouse S2 | Mouse S3 | Human S1 | Human S2 | Human S3 |
|---|---|---|---|---|---|---|
| Read split-half, full panel | 0.86 | 0.87–0.91 | 0.82–0.83 | 0.37 | 0.64–0.68 | 0.58–0.63 |
| Gene- and read-disjoint (F5b) | 0.69–0.70 | 0.74–0.83 | 0.66–0.68 | 0.20–0.22 | 0.47–0.48 | 0.35–0.41 |
| Same, notebook 36's 160 CAC anchors (F5c) | 0.29–0.34 | 0.38–0.45 | 0.36–0.40 | 0.05–0.06 | 0.19–0.20 | 0.24–0.35 |
| Neighbour, read-disjoint (r_nb) | 0.17–0.42 | 0.40–0.47 | 0.51–0.59 | 0.15 | 0.28–0.32 | 0.19–0.36 |
| Neighbour, gene/read/polygon-disjoint | 0.17–0.35 | 0.37–0.43 | 0.41–0.47 | 0.11 | 0.20 | 0.12–0.25 |
| Random partner | −0.04 to 0.04 | ~0 | ~0 | ~0 | ~0 | ~0 |
| Neighbour log library (technical) | 0.10–0.16 | 0.05–0.10 | 0.14 | 0.11–0.20 | 0.23–0.25 | 0.29–0.30 |
| S vs log library | 0.09–0.24 | 0.06–0.08 | 0.05–0.27 | −0.18 / −0.25 | −0.18 / −0.21 | −0.23 / −0.01 |
| SCF13 vs landmark half-score | 0.13–0.33 | 0.56–0.64 | −0.12 to 0.22 | 0.02–0.13 | 0.37–0.49 | 0.05–0.34 |
| DPT13 vs landmark half-score | 0.05–0.20 | 0.61–0.73 | 0.50–0.56 | −0.08 to 0.06 | 0.45–0.55 | 0.30–0.48 |

**F6 · Morphology** (`f6_morphology.csv`). Short side versus the residual landmark score:
|ρ| ≤ 0.20 everywhere. Short side versus log library: 0.60–0.83. Morphology carries little
within-segment positional information; it drives exposure.

**F7 · Anchors (nb42).** Morphological glomeruli: 217, 168, 127 and 202. PT within 2 µm: 46 mouse and
9 human; within 10 µm: 543 mouse and 237 human. S1 share among attached PT is 59–71%, so the anchor is
weak.

---

## 5. Candidate methods

Comparator at zero cost, used in every arm: **LRS, the landmark ratio score**. S above, built from
the selection-reference panel; within each gene fold, the landmarks of the scored fold are withheld.
A new model must beat this trivial score to earn its complexity.

### C1 · Landmark count-position model (LCP). Recommended first

- **Mechanism.** Positions on a grid z ∈ {0, …, 1} with K = 50 points.
  - *Likelihood.* For structure i of species d and landmark gene g:
    y_ig ~ NB(ℓ_i·exp(η_dg(z_i)), θ_g), where ℓ_i is the all-gene library (an exposure offset).
    Tempered by 1/T, with T ≥ 1 calibrated by count-split predictive likelihood on the training
    specimen only.
  - *Rate paths η_dg.* Spline, df 6. Initialized and shrunk (fixed λ) toward an external segment-level
    shape: mouse microdissection S1, S2 and S3 means placed at equal gauge thirds; human early and S3
    only.
  - *EM.* The E-step gives the posterior over the grid. The M-step refits η from posterior-weighted
    counts with **equal specimen weight at each grid point**. This is the repeated-object prior: one
    shared path for all nephrons of a species. No offsets, hierarchy, covariance or flow.
  - *Labels.* Never used in fitting. They remain an independent check.
  - *Outputs.* Posterior mean, 80% interval and entropy. Downstream analyses average over assignments
    (notebook 27 precedent).
- **How structure enters.**
  - External segment atlases fix the axis and the gauge (landmarks).
  - The exposure offset removes cut size.
  - Specimen-balanced shared rates impose the repeated object.
  - Species registration is through conserved landmark shape, not Harmony.
- **Why it should fix named failures.**
  - The mouse S1/S3 within-segment order is present in landmark counts (F5b 0.66–0.83) but unused by
    SCF13/DPT13: P1w and gene-fold failures.
  - The offset targets P5 and the human S1 depth dependence.
  - An explicit shared landmark gauge targets the registration gap and peak dependence.
  - CAC's +0.9% P4 shows that landmark information predicts held-out genes.
- **Expected failure modes.**
  - Inherits external-atlas definitions: "position" means progress through the external S1→S3
    program, not physical length.
  - Within-segment landmark variation may partly be nephron type, cell state or spillover. The panel
    is not yet filtered for PT specificity; phase 2 applies notebook 35's PT-specific mask.
  - Human within-early order has no human external standard (early versus S3 only), and human S1 is
    measurement-limited (F5b 0.20).
  - The independent-gene NB likelihood is overconfident without T.
  - P1 and P1w are partly aligned with the construction; evaluation references must stay disjoint
    from selection references.
- **Cost.** One new module (`pseudospace/landmark_position.py`) with synthetic tests. Each fit takes
  seconds to minutes; about 40 fits for folds and perturbations.
- **Abandon if**, in mouse (where counts are adequate):
  - gene-fold within-segment agreement is not above DPT13 (0.59), **or**
  - P1w is not above DPT13 (0.20), **or**
  - P4c (section 6) is not > 0 in ≥ 5 of 6 cells, **or**
  - it fails to beat LRS on P4c and stability (the model's complexity is then unearned).

### C2 · Tissue-coupled LCP (LCP-T). Conditional on C1, and on approval to use x/y in inference

- **Mechanism.** C1 plus a robust pairwise prior over PT polygons ≤ 5 µm apart:
  - ψ(z_i, z_j) = (1 − π) + π·exp(−(z_i − z_j)²/2τ²);
  - mean-field updates on the polygon graph;
  - π, τ and the coupling weight set by count-split held-out-read likelihood on training data, or
    fixed a priori (π = 0.5, τ = 0.05).
- **How structure enters.** Physical adjacency: each mouse PT occupies its own cortical domain (Zhai
  et al. 2003), and pars recta run in bundles. Neighbours carry disjoint-read positional signal (F5b,
  0.11–0.47).
- **Why.** Human S1 is measurement-limited (0.20). Pooling neighbours with ≥ 0.11 disjoint
  information is the only available extra information source. It should also stabilize gene folds.
- **Expected failure modes.**
  - It imports nephron-type variation (Ransick et al. 2019) and spatial tissue fields (neighbour
    library correlation 0.05–0.30).
  - It smooths across segment boundaries.
  - It mechanically raises gene-fold agreement.
  - It makes every spatial validation circular.
- **Mandatory controls.**
  - *Placebo graph:* the same degree sequence with edges rewired within the same label pair, which
    removes physical adjacency but keeps label-level smoothing.
  - *Spatial-artefact injection:* a smooth random field, wavelength 500 µm and amplitude log 1.5,
    multiplied into 20% of landmarks balanced late/early, then refitted.
- **Cost.** One graph module plus mean-field. Cheap; graph construction is done (F4 code).
- **Abandon if:**
  - the placebo graph gives ≥ 70% of the real graph's gain over C1 on P4c or gene-fold agreement, **or**
  - injected-field absorption gives within-segment |ρ(Δz, field)| > 0.10, **or**
  - human S1 gene-fold agreement improves by < 0.10 over C1.

### C3 · Stereological gauge (SG). Cheap presentation layer, optional

- **Mechanism.** Map any coordinate to an estimated fraction of PT length. Within each specimen, take
  the cumulative occupancy of profiles weighted by aspect ratio (≈ 1/|cos φ|; the orientation
  correction rests on Nyengaard 1999 and is to be verified in full text before use). A
  Halpern-style geometric prior, but used only as a monotone gauge, never as a constraint on order.
- **Fixes.** Peak positions that depend on the coordinate's arc-length gauge; segment-unit reporting.
- **Failure modes.** Human sections are cortex only, so S3 is truncated; anisotropy is
  specimen-specific; the AR is affected by segmentation merges.
- **Abandon if** gauge-transformed peaks of the 22 primary pathways are no more concordant between
  SCF13 and DPT13 than segment-quantile units (Spearman gain < 0.05).

### C4 · Within-section gradients. Rejected in phase 1

F3 failed its pre-registered bar. It is not pursued.

---

## 6. Pre-registered phase-2 evaluation

**Arms.**
- Comparators: SCF13 and DPT13 (frozen nb13 coordinates and their nb36 refits), segment labels, the
  segment + coverage quadratic oracle, CAC (reported from nb36), and LRS.
- New arms: C1 (LCP); C2 (LCP-T), C2-placebo and C2-injected (only if C1 is not abandoned).

**Machinery.** `pseudospace/coordinate_evaluation.py` and the nb35/36 code paths: sha256 gene folds;
selection versus evaluation references exactly as in nb35; `fit_segment_covariate_oracle`.

**Construction rules.**
- Landmarks come from selection references only, with notebook 35's PT-specific mask.
- Scored genes are withheld from construction by fold.
- K = 50, df 6 and λ are fixed in advance.
- T and the C2 coupling are estimated only from count-split predictive likelihood within the training
  data.
- No depth, glomerulus distance or x/y-derived quantity validates or selects anything.

**Metrics.** Each is per species, with the specimen, fold and direction cells shown.

1. **P4**, as in nb36: held-out folds and cross-specimen prediction, gain over segment means.
2. **P4c (new, primary for "continuous beats discrete").** The gain from adding a within-segment
   spline of the coordinate to the segment + coverage quadratic oracle, on held-out-fold genes,
   training on one specimen and predicting the other. Also reported with the Poisson deviance and a
   library offset, as an exposure-aware sensitivity.
3. **P1** and **P1w**, fold-protected, evaluation references only.
4. Gene-fold within-segment agreement (median over pairs, per segment).
5. Specimen stability:
   - LOSO within-segment agreement;
   - cross-specimen consistency: rates fitted on specimen A, positions in B, agreement with B's own
     fit;
   - count-split half-A against half-B coordinates, within segment.
6. P5 (equal depth; human S1 shown separately), and the within-segment ρ(position, log library).
7. P2 injected-gradient absorption; for C2, also spatial-artefact absorption.
8. Registration gap across fold, depth and count-split refits; human S1 < S2 < S3 kept.
9. Downstream specificity: notebook 37's joint T_spatial test with relabeling (NCDR and the
   relabeled-call share), run only on the arm that passes 1–8. Overlap with the 22 primary pathways
   reported, not a criterion.
10. Uncertainty (reported, not gated):
    - the fraction of structures whose half-B MAP falls in the half-A 80% interval;
    - interval width by species and segment;
    - downstream T_spatial under posterior sampling.
11. Exploratory and not gating: within-segment rank distance of the 67 + 58 mouse collinear
    same-tubule pairs, against touching non-collinear neighbours. Expression-only arms only, because it
    is circular for C2.

**Adoption rule.** Adopt C1 (or C2) as the paper's reconstruction only if all of the following hold.

- **(A) Non-inferiority.**
  - P1 ≥ SCF13 − 0.03 in each species;
  - P5 ≥ 0.80 overall and ≥ 0.70 in human S1;
  - registration gap ≤ 0.15;
  - P2 median ratio ≤ 1.10;
  - downstream joint test NCDR < 0.25 and relabeled calls < 5% of tested pathways.
- **(B) The named failures improve.**
  - P4c > 0 in ≥ 5 of 6 mouse cells and ≥ 4 of 6 human cells, and above the best of SCF13, DPT13 and
    LRS in mean;
  - P1w ≥ DPT13 + 0.05 in mouse and ≥ DPT13 in human;
  - gene-fold within-segment median ≥ 0.70 in mouse and ≥ 0.59 in human;
  - LOSO ≥ 0.85.
- **(C) C2 only.** Its gain over C1 exceeds the placebo graph's gain over C1 by ≥ 0.10 in gene-fold
  agreement and by > 0 in P4c, with injected-field absorption ≤ 0.10.

**If nothing qualifies.** Keep the SCF13/DPT13 robust-list framing. Report the measurement analysis as
the result: within-segment order is identifiable in mouse, but not by generic trajectories, and is
measurement-limited in human S1.

**Order of work** (cheap first; stop at the first abandon rule that fires):
- **51:** landmark panel, LRS comparator and P4c for all existing arms. P4c is needed anyway for the
  outline's "continuous versus discrete" claim.
- **52:** C1 fits and metrics 1–8.
- **53:** C2 with its placebo and injection controls, only if C1 survives and x/y use is approved.
- **54:** the downstream joint test and uncertainty propagation on the adopted arm.
- **55:** C3 gauge and figures.

---

## 7. Recommendation

Start with **C1 (LCP), alongside the zero-cost LRS comparator and the P4c metric** (notebook 51,
then 52).

It targets the clearest diagnosed failure. In mouse the counts already contain a within-segment
landmark axis that is reproducible across disjoint genes and reads (0.66–0.83), and the generic
coordinates ignore it in S1 and S3 (SCF13 −0.12 to 0.33). It is cheap and label-free. It needs no
change to the repository's rules on x/y. It replaces Harmony and the curve, the components implicated
in exposure leakage and projection instability, and is not a rerun of any closed hypothesis.

C2 should follow only if C1 survives and the coordinator or user approves using tissue adjacency in
inference. It is the only route to the human S1 measurement limit. Its neighbour signal is real (≤ 0.47
disjoint, against 0 for random partners), but it carries the highest risk of confounding.

**Main risks.**
- Landmark position equals "external program progress", not physical length.
- Nephron-type and state variation can masquerade as position, in both C1 and C2.
- Human within-early order lacks an external standard.
- The paper's continuous claim may end up mouse-only.

---

## 8. Side findings for the coordinator (no existing file changed)

1. **Human pixel calibration.** The human GeoJSON polygons are in the Space Ranger full-resolution
   frame at 0.4408 µm/px. Polygon area divided by n_spots × 4 µm² is 1.000 at that scale, but 1.287 at
   the 0.5 µm/px that nb42 read from the QuPath measurements (`f2_pixel_calibration.csv`).
   - Human µm values in nb42 and the supplementary text are overstated: distances by about 13%
     (medians of 149/168/259 µm should be about 131/148/228), areas by 29%.
   - The "attached" thresholds were effectively 2.3 and 11.3 µm.
   - Ranks and odds ratios are unaffected. Suggested fix: in nb42 take µm/px from
     `scalefactors_json.json`.
2. **Library size is mainly cut size.** Polygon short side versus log library: ρ 0.60–0.83. This is
   worth a sentence where the paper discusses depth dependence and equal-depth thinning.
3. **P4 is dominated by coverage.** The coverage oracle beats every coordinate by 0.4–1.8%. The
   outline's bar "predicts held-out genes no better than segment labels" should be read with that,
   and P4c added.
4. **A heavy tail of strong within-polygon axial gradients in mouse** (F3; mean z² 4.3) may flag
   merged or boundary-spanning polygons. One line; out of scope.

## 9. Files created

- `workstream-A scratch: plan_phase1.md` (this file)
- `workstream-A scratch: feasibility.py`, `workstream-A scratch: feasibility_f5b.py`, `workstream-A scratch: feasibility_f5c.py`,
  `workstream-A scratch: feasibility.log`
- `results/pt_reconstruction_v3/ws7_phase1_feasibility/`: `f1_f2_polygon_table.csv`,
  `f1_f2_summary.csv`, `f2_nspots_guard.csv`, `f2_pixel_calibration.csv`, `f3_*.csv` (4 files),
  `f4_neighbours.csv`, `f4_candidate_pairs.csv`, `f5_neighbour_information.csv`,
  `f5b_gene_split.csv`, `f5c_cac_anchor_ratio.csv`, `f6_morphology.csv`, `gene_counts.json`
  (gitignored)
