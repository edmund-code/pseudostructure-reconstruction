# An independent IRI benchmark for the AKI calls (notebook 94)

**The question.** Our AKI calls cannot be judged by our own design. With two IR and two control mice:
- the screens fail calibration (0.120 in notebook 71, 0.144 in notebook 76);
- the relabeling aligned with the A2/A4 factor carries most false calls;
- transport loss, its progression and the location of the cell-cycle program were "not claimable" (notebook 76).

Which method's AKI calls hold in an independent ischemia–reperfusion experiment with more mice? The candidates are:
- the positional F-test;
- clustering;
- position-split ORA;
- notebook 91's level and location screen;
- the published programs.

This notebook builds that benchmark from Kirita et al. 2020 and scores every existing AKI call list on it. Later
notebooks score new methods on the same benchmark.

Protocol frozen on 2026-10-18, before any comparison of our calls with the external data. At freezing time the
external files had been downloaded, and only their metadata (groups, replicates, author cell types, nucleus counts)
had been tabulated. No external expression value had been compared with any of our statistics.

## 0 · The external experiment (read in full text, PMC7355049)

Kirita Y, Wu H, Uchimura K, Wilson PC, Humphreys BD. *Cell profiling of mouse acute kidney injury reveals conserved
cellular responses to injury.* PNAS 2020;117:15874 (PMID 32571916; GEO GSE139107).
- **Mice and injury:** male C57BL/6J mice, 8–10 weeks old. Bilateral IRI: the renal pedicle was clamped for 18 min at
  a body temperature of 36.5–37.5 °C. Controls had sham surgery.
- **Design:** sham, then 4 h, 12 h, 2 d, 14 d and 6 weeks after IRI. Three mice per group (Fig. 1A, "n = 3 mice per
  group").
- **Libraries:** GEO lists 24 libraries. Replicate labels such as `1_1` and `1_2` are two libraries of mouse 1, which
  gives 18 mice.
- **Assay:**
  - snRNA-seq, 10x 3′ v2;
  - zUMIs counts over exonic and intronic reads;
  - nuclei with 150–8,000 genes, mitochondrial share < 0.01 and ribosomal share < 0.06;
  - SoupX in the authors' analysis.
- **Author cell types used here:**
  - PTS1, PTS2 and PTS3: healthy segment identities;
  - NewPT1: "Injured PT", mostly the three acute injury states;
  - NewPT2: mostly failed-repair PT (FR-PTC).

  Injured nuclei therefore leave the segment labels.

**Why this dataset.** It has the same sex and species as our mice, and an IRI model. Mice, not nuclei, are the
replicates. It has segment labels, and enough mice for a moderated test. Notebook 80 used only its published gene
lists (Dataset S2), never its counts.

## 1 · Labels

- Notebook: `analysis/notebooks/94_pt_external_iri.ipynb`; logic version `94.external_iri.1`; seed 94.
- Results: `results/pt_external_iri/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`).
- New module: `pseudospace/external_replication.py`, with `tests/test_external_replication.py`.
- External files in `data/external/kirita2020_iri/` (outside Git). Their sha256:

  | File | sha256 |
  |---|---|
  | GSE139107_MouseIRI.metadata.txt.gz | ac607fa70c7ee0c6b269009893cb0ca06b463b3c846db613d7dec6a2b5ee9a11 |
  | GSE139107_MouseIRI_control.dge.txt.gz | 626a4a354bd0263cbf68f91cb2de1304e629a7bfa9664aa34cbb899fb3427628 |
  | GSE139107_MouseIRI_4hours.dge.txt.gz | 21ef8eae053a069f2e31a0cdf870aeeabf4a90d83b33d37ce06545bb39486321 |
  | GSE139107_MouseIRI_12hours.dge.txt.gz | 5bb3a2bae56c90fb854e655d13c1ce032f0cdf5c14ab97a06521ef863769f732 |
  | GSE139107_MouseIRI_2days.dge.txt.gz | 652b258719afea520a5c44b13f7bb4e4a964b549242e5ec921498c20c23dbdd6 |
  | GSE139107_MouseIRI_14days.dge.txt.gz | a0a249ddb5b646f882ef4e9955821c4e71e40b23ac9df2b8999780889151d6bc |
  | GSE139107_MouseIRI_6weeks.dge.txt.gz | 7b80febcc52a936ffa42e8930debed2d793fcfb4d2896de415f2b45a6792c3af |

- Our inputs are read only and never refit:
  - notebook 71's AKI tables and clustering stage cache;
  - notebook 76's tables;
  - notebook 91's gene and set tables;
  - notebook 93's call table;
  - notebook 80's published program members (`docs/results/known-programs/programs_aki_members.csv`).

## 2 · Guards (not cached)

- G1: the seven files match the sha256 above.
- G2: the metadata reproduces the frozen tabulation:
  - 126,578 nuclei;
  - 24 libraries and 18 mice;
  - PT nuclei per group:

    | Group | PTS1 | PTS2 | PTS3 | NewPT1 | NewPT2 |
    |---|---|---|---|---|---|
    | sham | 4,750 | 4,106 | 1,937 | 153 | 54 |
    | 4 h | 1,373 | 2,293 | 460 | 2,522 | 264 |
    | 12 h | 1,273 | 1,841 | 361 | 2,783 | 355 |
    | 2 d | 4,062 | 3,580 | 516 | 806 | 529 |
    | 14 d | 1,883 | 2,481 | 1,057 | 713 | 2,062 |
    | 6 wk | 2,897 | 2,812 | 1,511 | 189 | 552 |

- G3: every metadata barcode is a column of exactly one count matrix, and every matrix holds whole-number counts.
- G4: our read-only inputs exist, and their call counts equal the recorded ones:
  - notebook 71: F 103 AKI calls;
  - notebook 93: 459 clustering calls and 134 position-split ORA calls on AKI;
  - notebook 76: 106 F calls;
  - notebook 91: 3,159 level-change and 362 location-change genes, 538 level sets and 236 location sets.

## 3 · External statistics (fixed now)

**Pseudobulks.**
- PT nuclei are the five author types above.
- For each mouse × type, counts are summed over that mouse's libraries.
- Library size = the summed counts over every gene in the matrix.
- Whole-PT pseudobulk = the sum over the five types.
- log2 CPM = log2((count + 0.5) / (library + 1) × 10⁶).

**Gene universe E.**
- Candidates are the genes in at least one of our tested AKI universes:
  - notebook 71's AKI universe;
  - notebook 76's U76;
  - notebook 91's shared, induced and lost genes.
- A candidate needs an exact symbol match in Kirita's matrices.
- It also needs whole-PT CPM ≥ 1 in at least 3 of the 18 mice.
- Mitochondrially encoded genes are excluded, as in our analyses.

**Level, per time point t (4 h, 12 h, 2 d, 14 d, 6 wk).**
- Fit: limma-trend with robust eBayes (`limma_trend_species`) on whole-PT log2 CPM, 3 IRI-t mice against 3 sham mice.
- The coefficient is centred by its median over E, as notebook 91 centres by the median reference gene.
- The output is a centred t, a two-sided p, and BH q over E within t.

**Location, per t.**
- For each mouse, d = log2 CPM(PTS3) − log2 CPM(PTS1).
- d goes through the same limma-trend fit, IRI-t against sham.
- The output is a centred coefficient, t and q. Positive means the injury effect is larger, or less negative, in S3
  than in S1.
- S3 − S2 and S2 − S1 are fitted the same way and reported.
- This is the location of injury in identity-retaining nuclei: nuclei that became NewPT are not in a segment group.

**Location sensitivity (assigned injured nuclei).**
- Markers come from the sham mice: for each of S1, S2 and S3, the 50 genes in E with the largest one-vs-rest log2 CPM
  difference, among genes with CPM ≥ 10 in that segment.
- Each NewPT1 or NewPT2 nucleus goes to the segment whose sham centroid has the highest Spearman correlation with it.
  The profiles are log1p CP10k over the markers, and the centroid is the mean over that segment's sham nuclei.
- **Accuracy check:** leave one sham mouse out. Markers and centroids come from the other two sham mice, and the
  held-out mouse's PTS nuclei are classified. The accuracy is reported.
- Segment pseudobulks with the assigned nuclei added are refitted as above.

**External set statistics.**
- Set test: limma's cameraPR (`camera_pr`, inter-gene correlation 0.01, BH within t) on each external gene t vector:
  - the level t;
  - the location t;
  - |location t|.
- Sets: notebook 71's AKI gene sets (sizes 10–300 in E), notebook 76's sets and notebook 91's sets. They are
  library sets with the same identifiers. Each set is tested once on its members in E.
- A signed set z is Φ⁻¹(1 − p/2) carrying camera's direction. |location| has only an "Up" call.
- Programs (claim checks, S-P80) are tested on their members in E when they keep ≥ 5 of them (amended before any
  fit: the library-set size rule cannot apply to programs such as the 12-gene lost program).

**Time-point alignment.**
- For each t, Spearman's ρ over E ∩ notebook 91's shared genes, between our L (notebook 91, condition) and the external
  level coefficient.
- The primary time point t* has the largest ρ.
- The same is reported for location: our P_late = Δ late − L against the external location coefficient.
- **Caveat:** t* is chosen on levels, which slightly favours level calls. Every result is therefore also shown at
  every t.

## 4 · Our internal signs (the same for every list)

Every gene or set gets one internal level sign and one internal location sign from notebook 91. A list's own method
need not be signed.
- **Genes:** level sign = sign(L). Location sign = sign(Δ late − L), meaning the late-PT effect beyond the PT-wide level.
- **Sets:** level sign = sign of 91's set L. Location sign = sign of 91's set P averaged over the bins centred at
  ≥ 1.8 (91's LATE, the bins its genes' Δ late uses). *Amended before any fit: the frozen text said 2.1 and called it
  91's LATE, but 91's LATE is 1.8; the rule's intent (91's LATE, the same bins for genes and sets) is kept.*
- A gene or set absent from 91 is scored only on the unsigned location axis.

## 5 · The call lists (frozen)

| List | Source | Items | Primary axis |
|---|---|---|---|
| G-L91 | 91 `level change` genes | genes | level |
| G-M91 | 91 `location change` genes | genes | location |
| G-F71 | 71 condition gene q ≤ 0.05 (the positional F-test's gene step) | genes | unsigned location |
| G-C71 | clustering condition run (71's cache, as read in 93): selected genes with `direction` Up or Down (crossover excluded) | genes | level |
| S-F71 | 71 positional F-test calls (cameraPR Up, FDR ≤ 0.05) | sets | unsigned location |
| S-F76 | 76 positional F-test calls on U76 | sets | unsigned location |
| S-C71 | clustering calls (93's table, method `clustering`) | sets | level |
| S-N93 | position-split ORA calls (93's table) | sets | location |
| S-L91 | 91 `level set` | sets | level |
| S-M91 | 91 `location set` | sets | location |
| S-P80 | 80's published injury programs except Kirita's own (K-*), signed by the published direction | programs | level |

Kirita's own programs (K-INJ12, K-INJ3, K-SEV, K-REP, K-FR) are excluded, because they were derived from these data. If
a list's sign comes from its own method (S-P80's published direction, G-C71's group), that sign replaces the internal
level sign on the level axis.

## 6 · Scoring (fixed now)

**Replicated**, per item and axis at a time point:
- **level:** the external level q ≤ 0.05 (genes) or camera FDR ≤ 0.05 (sets), with the external sign equal to the
  item's level sign;
- **location:** the same, with the external location test and the location sign;
- **unsigned location:** genes need the external location q ≤ 0.05 in either direction; sets need the |location t|
  camera Up with FDR ≤ 0.05.

**Per list, axis and time point:**
- n: items scorable on that axis;
- k: replicated items; precision = k / n;
- **matched null:** 1,000 random lists of size n, drawn without replacement from the same pool. The pool is the list's
  own tested universe ∩ E (genes) or its method's tested sets ∩ the externally tested sets (sets).
  - Genes are matched on the decile of control expression (91's log2 CP10k; 76's or 71's mean count when absent).
  - Sets are matched on size tercile.
  - The null's mean, 2.5–97.5% band and empirical p = (1 + #{null ≥ k}) / 1,001.
- **Enrichment** = precision / null mean.

**Threshold-free agreement** (all items, not only calls): Spearman's ρ between an internal statistic and the external
statistic of its axis.
- Genes:
  - 91's L against the external level coefficient;
  - P_late against the external location coefficient;
  - 71's F z against |external location t|.
- Sets:
  - 91's z L against the external signed level z;
  - 91's z M · location sign against the external signed location z;
  - 71's F set z (from its one-sided p) against the external |location| z.

## 7 · Decision rules (fixed now)

1. **A list replicates externally** if, at t*, its precision on its primary axis exceeds the null's 97.5th percentile
   with n ≥ 10.
2. **Comparison on one axis.** List A beats list B if both have n ≥ 20 and A's precision exceeds B's by ≥ 0.10.
   Neither "beats" otherwise.
3. **A comparison is stable** if it holds at t* and at one or more other acute time points (4 h, 12 h, 2 d).
4. **The external verdict on the positional F-test** for AKI shape calls:
   - "replicates" if rule 1 holds for S-F71 and G-F71 on the unsigned location axis;
   - "partly" if one of them holds;
   - "does not replicate" if neither holds.
5. **The best AKI calling method on an axis** beats every other list with n ≥ 20 on that axis (rule 2) and is stable
   (rule 3). Otherwise there is no single best method on that axis.

## 8 · Claim checks (named now; a miss is reported, not tuned)

External camera on the axis named. "Holds" means FDR ≤ 0.05 in the stated direction at t*. Every t is shown.

| Check | Program | Axis, direction | Our status |
|---|---|---|---|
| K1 | 76's 18 injury markers (positive control) | level up | detected; S3-concentrated |
| K2 | 76's cell-cycle program C | level up | claimable (A4) |
| K3 | 76's transport program T_PT | level down | not claimable (A3) |
| K4 | T_PT, progressive | location down (more loss in S3) | not claimable |
| K5 | 91's S3-located losses: PPAR signalling, peroxisome, peroxisomal lipid metabolism, cholesterol biosynthesis, plasma lipoprotein assembly, Hallmark fatty-acid metabolism | location down | 91: location sets peaking in S3 |
| K6 | 91's S3-located gains: Hallmark G2-M, E2F, TNF-α/NF-κB, p53, apoptosis, inflammatory response, hypoxia | location up | 91: location sets peaking in S3 |
| K7 | Xiong 2023 male-biased PT program (91's definition) | level down | 91: strongest external signal |
| K8 | fasting-up and fasting-down programs (91's definitions, GSE267280) | level up / level down | 91: IR mice carry a fasting signature |
| K9 | 91's induced (28) and lost (12) programs | level up / level down | 91: on/off along the PT |

K1–K3 and K7–K9 are level checks of whole-PT pseudobulks. K4–K6 rest on identity-retaining nuclei. Notebook 91's
library-set names are matched by exact set identifier.

## 9 · Outputs

**Tables (`tables/`):**
- `nuclei.csv`: nuclei per mouse × type, and libraries;
- `universe.csv`: E, with the source universes;
- `external_genes.csv.gz`: per t, level and location (S3−S1, S3−S2, S2−S1) coefficients, t and q, plus the
  assigned-nuclei location;
- `external_sets.csv.gz`: per t, camera for level, location and |location|;
- `assignment_accuracy.csv`;
- `alignment.csv`;
- `list_scores.csv`: list, axis, t, n, k, precision, the null band, enrichment and p;
- `comparisons.csv`;
- `agreement.csv`;
- `claim_checks.csv`;
- `how_the_protocol_was_read.csv`;
- `run_record.json`.

**Figures (plots only):**
- 94a: alignment ρ by time point, for level and location;
- 94b: precision against the null band for every list on its primary axis at t*, with small multiples over t;
- 94c: claim checks, external z by time point;
- 94d: internal against external set statistics for the F-test's sets and 91's sets.

## 10 · Caveats that travel with every result

- **Platforms differ.** Kirita profiles nuclei (exonic and intronic reads, 3′ v2); we profile Visium HD probe counts
  in tubule units. Genes whose nuclear and cytoplasmic RNA differ will disagree for that reason alone.
- **Our time point is unknown.** t* is a data-driven alignment, not a fact about our mice.
- **Kirita's injured nuclei leave the segment labels,** so its location is the location of injury in
  identity-retaining cells. The assigned-nuclei sensitivity bounds this.
- **Different labs, surgeries and clamp times.** Agreement is a lower bound on shared biology. A miss does not show
  that a call is wrong in our mice.
- External replication tests whether an injury effect generalises. It does not test whether a method is calibrated in
  our design.

## Notebook 94 results

`analysis/notebooks/94_pt_external_iri.ipynb` (logic `94.external_iri.1`) writes to `results/pt_external_iri/`. It
followed the protocol frozen in 9da52b7 and amended before any fit in 1599681. Guards G1–G4 passed. A cold run takes
about 10 min, almost all of it parsing the six dense matrices. Module `pseudospace/external_replication.py` has
`tests/test_external_replication.py` (4 tests).

**How the protocol was read** (`how_the_protocol_was_read.csv`):
- the three pre-fit amendments;
- set members are the union of a set across notebooks 71, 76 and 91, in E;
- programs are tested as their own camera family;
- S-P80's null pool is 91's library sets;
- expression strata come from 91's control expression, else clustering's;
- K5 uses the Hallmark fatty-acid set as frozen. 91's location set of that name was Reactome's.

Two fixes after the first development run:
- clustering's direction labels are lowercase;
- notebook 76's `programs.csv` truncates members at 60 genes, so C (297 genes) and T_PT (93) were rebuilt from 76's
  sets and transporter table and checked against 76's counts.

After execution, the pseudobulk stage's cache key was corrected: a list was passed as `code`, so the key held a
memory address and never hit. Results are unchanged.

**The external data.**
- 18 mice. E has 10,733 genes and 1,476 library sets.
- At 12 h, 4,014 genes change level at q ≤ 0.05; at 4 h 1,829, at 2 d 1,554, at 14 d 689, at 6 weeks 3.
- **Location is weak externally.** Only 45 genes (12 h) pass q ≤ 0.05 for S3 − S1, and none at any other time.
  Identity-retaining nuclei, 3 vs 3 mice.
- Assigning the injured nuclei to segments works: leave-one-mouse-out accuracy is 0.90–0.92, with S2 recall 0.81–0.87
  the weakest. Most injured nuclei at 4–12 h go to S2 and S3.

**Alignment.** Our L against the external level coefficient:

| | 4 h | 12 h | 2 d | 14 d | 6 wk |
|---|---|---|---|---|---|
| level ρ | 0.431 | **0.524** | 0.407 | 0.302 | 0.128 |
| location ρ | 0.182 | 0.254 | 0.274 | 0.060 | 0.053 |

t* = 12 h. Our IR mice look most like acute IRI, 4 h to 2 d after reperfusion, which matches notebook 80's reading of
Kirita's state lists. The cell-cycle program (K2) points later, to 2 d.

**Lists on their primary axis** (replicated / scored; enrichment over the matched null; R = rule 1 holds):

| List | 4 h | 12 h (t*) | 2 d | 14 d |
|---|---|---|---|---|
| G-L91 · level | 789/3071 (1.8, R) | 1585/3071 (1.7, R) | 605/3071 (1.8, R) | 310/3071 (2.0, R) |
| G-C71 · level | 1043/5190 (1.4, R) | 2158/5190 (1.4, R) | 823/5190 (1.5, R) | 381/5190 (1.5, R) |
| S-L91 · level | 192/537 (1.6, R) | 306/537 (1.7, R) | 142/537 (1.8, R) | 107/537 (1.6, R) |
| S-C71 · level | 145/459 (1.3, R) | 238/459 (1.5, R) | 97/459 (1.3, R) | 71/459 (1.2) |
| S-P80 · level | 15/19 (2.9, R) | 16/19 (2.0, R) | 16/19 (4.5, R) | 16/19 (5.5, R) |
| G-M91 · location | 0/352 | 18/352 (7.9, R) | 0/352 | 0/352 |
| S-M91 · location | 10/236 (3.3, R) | 6/236 (5.3, R) | 1/236 | 0/236 |
| S-N93 · location | 8/134 (4.8, R) | 3/134 (5.0, R) | 1/134 | 0/134 |
| G-F71 · unsigned location | 0/4239 | 37/4239 (1.8, R) | 0/4239 | 0/4239 |
| S-F71 · unsigned location | 4/103 (6.5, R) | 9/103 (6.8, R) | 0/103 | 0/103 |
| S-F76 · unsigned location | 4/105 (6.2, R) | 10/105 (8.0, R) | 0/105 | 0/105 |

**Verdicts.**
- **Rule 1.** Every list replicates on its primary axis at t*.
- **Rule 4.** The positional F-test's AKI shape calls **replicate**, as sets and as genes.
  - Its sets have the highest enrichment on the unsigned location axis: 6.8× and 8.0×, against 1.8× for its genes.
  - Clustering's sets do not replicate on that axis: 3/459, 0.5×.
- **Rule 5: no single best method on any axis.**
  - On the level axis the stable wins (rules 2–3) are 91's level sets over the positional F-test's sets (0.570 against
    0.392 and 0.375) and over the F-test's and clustering's gene lists. 91's level sets do not beat clustering's sets:
    0.570 against 0.519, inside the 0.10 margin.
  - On the location axes no list can beat another by 0.10, because external location is weak (every precision is
    below 0.2). The frozen margin cannot separate methods there.

**Threshold-free agreement at 12 h** (Spearman over all items):

| Pair | ρ |
|---|---|
| sets · 91 z L vs external level z | 0.704 |
| sets · F z vs external \|location\| z | 0.537 |
| genes · L vs external level | 0.524 |
| genes · Δ late − L vs external location | 0.254 |
| sets · 91 z M × sign vs external location z | 0.247 |
| genes · F z vs \|external location t\| | 0.187 |

The positional F-test's set ranking predicts which sets change shape in another lab's IRI, at ρ 0.54 (0.34 at 4 h,
0.39 at 2 d). That is a stronger external result than its call list shows.

**Claim checks** (external camera z; holds at FDR ≤ 0.05 in the stated direction):

| Check | 4 h | 12 h | 2 d | 14 d | 6 wk | Our status → external |
|---|---|---|---|---|---|---|
| K1 injury markers up | 8.9 | 9.9 | 9.3 | 11.5 | 5.4 | holds at all times |
| K2 cell-cycle program up | 0.9 | 2.7 | **8.7** | 0.9 | 1.3 | holds at 12 h and 2 d |
| K3 transport (T_PT) down | −7.4 | **−11.1** | −5.2 | −6.0 | −1.8 | holds 4 h–14 d. 76: "not claimable" |
| K4 T_PT loss larger in S3 | −3.3 | −3.5 | −3.4 | 0.1 | −0.9 | holds 4 h–2 d. 76: "not claimable" |
| K5 peroxisomal lipid metabolism, S3 | −4.9 | −4.2 | −4.0 | −0.4 | −0.1 | holds 4 h–2 d |
| K5 peroxisome, S3 | −4.0 | −3.3 | −3.5 | 0.9 | 0.8 | holds at 4 h only (FDR) |
| K5 PPAR, cholesterol biosynthesis, lipoprotein assembly, Hallmark FA | | | | | | do not hold |
| K6 TNF-α/NF-κB, S3 | 5.4 | 3.7 | 2.8 | −0.3 | 1.5 | holds at 4 h and 12 h |
| K6 hypoxia, p53, S3 | 3.5, 3.8 | 3.4, 2.4 | | | | hold at 4 h only |
| K6 G2-M, E2F, apoptosis, inflammatory response, S3 | | | | | | do not hold |
| K7 male-biased PT down | −4.9 | −5.9 | −6.2 | −6.2 | −1.3 | holds 4 h–14 d |
| K8 fasting-down down | −6.8 | −7.8 | −4.5 | −4.6 | −0.5 | holds 4 h–14 d |
| K8 fasting-up up | 2.4 | 1.0 | 1.7 | 1.1 | 0.2 | holds at 4 h only |
| K9 induced up / lost down | 5.7 / −6.8 | 6.8 / −7.1 | 11.2 / −6.3 | 7.8 / −5.7 | 5.1 / −2.0 | hold at nearly all times |

**What this answers.**
1. **The positional F-test's AKI calls generalise,** although they fail our internal calibration (0.120).
   - Its shape calls replicate in an independent IRI experiment 6–8× above chance.
   - Its set ranking agrees with the external shape ranking at ρ 0.54.
   - Clustering's AKI sets replicate as level changes (1.5×) but not as shape changes.
2. **For level, notebook 91's level screen is the most externally consistent list.**
   - It beats the F-test sets on that axis; the F-test is not built for level.
   - It ties with clustering's sets (0.57 against 0.52).
   - Set level z agree at ρ 0.70.
3. **Two claims notebook 76 could not make hold externally.**
   - Transport-program loss: z −11 at 12 h, holding from 4 h to 14 d.
   - Its larger loss in S3: z −3.5 at 4 h to 2 d.
   - Our relabelings failed because of the A2/A4 factor, not because the effect is absent.
4. **91's "non-injury contributors" are part of IRI itself.** In Kirita's IRI, as in ours, the male-biased PT program
   and the fasting-down program fall from 4 h to 14 d. These programs belong beside the injury reading as features of
   the response, not as signs of a cohort artefact.
5. **Not supported externally:**
   - 91's S3-located losses in PPAR, cholesterol and lipoprotein handling;
   - its S3-located gains in cell cycle, apoptosis and inflammatory response.

   Of 91's S3-located sets, only peroxisomal lipid metabolism, TNF-α/NF-κB and (at 4 h) hypoxia, p53 and peroxisome
   hold.

**Limits.**
- External location rests on identity-retaining nuclei and 3 vs 3 mice. Only the 4–12 h contrasts carry location
  signal.
- The location verdicts are therefore early-injury verdicts, and the frozen 0.10 margin cannot rank methods on
  location.
- t* was chosen on levels.
