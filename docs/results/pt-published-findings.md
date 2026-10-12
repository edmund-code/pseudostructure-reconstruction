# Published human vs mouse kidney findings, re-read along the PT: notebooks 81 (our data) and 82 (their data)

Frozen on 2026-10-12, before any of the analyses below was computed. The user asked for this on 2026-10-12: "Based on
the findings of specific pathways, genes, conclusions and data they publish, reanalyze everything from a hypothesis
driven approach to see if we can elaborate on their findings or see things they missed using our method." The user
also said findings may be descriptive. An example: "we show this pathway is likely localized in a sharp switch at this
specific part of PT S2, unknown before."

**What of our data was seen before freezing.**
- Notebook 79's block-level results for the CellSpectra species blocks and the zonation programs (in
  `pt-known-programs.md`).
- Gene membership of notebook 71's universe.

Not computed before freezing: the segment-unique subsets, the per-gene classes, the CellSpectra function sets, the
Raghubar programs, the gene-level hypotheses, the sex-bias test, the continuous descriptors, and every Cheval and atlas
analysis.

**Frozen files.** They are committed in `docs/results/published-findings/`, and the notebooks assert their sha256:

| File | sha256 |
|---|---|
| `programs_lit.csv` (one row per program, with the gene accounting) | `feeee7503b3b522df5876b01c17855ff0e1bcd39430867cec7f2f623043a6e4e` |
| `programs_lit_members.csv` | `8bcbd2556f52b859ad4ef62ffd2c98c14d4821e784a68ec979ecfbebbab219c7` |
| `programs_lit_dropped.csv` (empty) | `01bbbc984159c871844813172e1421420f125d412cd5357d145e321231ec46a5` |
| `genes_lit.csv` (single-gene hypotheses with their predictions) | `111cca20e92708f7fbceccca6b5fa0d49453544bfb5a57c26336e8e3805a22d8` |

**External data** live under `data/external/` (public; not tracked), each with a `SOURCE.txt`:

| File | sha256 |
|---|---|
| `cheval2012/cheval2012_table_s1.csv` (Cheval 2012 Table S1, parsed; 3 typographic cells cleaned, listed in SOURCE) | `52ea96118eb2e3744c22515aa9c6a8bbee150a188d520d9e49c870100d6f7fb9` |
| `cheval2012/cheval2012_mouse_symbols.csv` | `33f65262448c0d53c52f5453bf84a3637b4834be41150a643fd2ad4142e4d30d` |
| `quickgo_2026-10-12/GO_0004467.tsv` | `736757d377bd3e0bb93cbabda7348a298efdb21298d17d75fd02f66cf81776ba` |
| `quickgo_2026-10-12/GO_0006525.tsv` | `14f4c1d77492109c06b4bbd69a24c51a350c78ab6ae85826a4699ba1f0339085` |
| `quickgo_2026-10-12/GO_0015849.tsv` | `70a2319b825e1d9f6c358040a730fb456cd86be47fc8b44b63ece830ee49f2af` |
| `quickgo_2026-10-12/GO_0030169.tsv` | `6dc83c33f6caa3ef481e30308557e8a3a9d7f7f5d44d88818e186fd31e662c3b` |
| `quickgo_2026-10-12/GO_0045056.tsv` | `afa07d490e0dd267853e78ababd7821cdd65799d1bd7cc26abb091d5cba36734` |
| `quickgo_2026-10-12/GO_0071071.tsv` | `1a0e08fa38a5914f6cd8f492de74c5e1e3be6dec0b8486183147653a58773510` |
| `xiong2023_sex_bias/PT-S1_fWTvsmWT.csv` | `4764f77435e1b232e2e686770776c16b3c9f93a085eee1b52d9914b116f00eee` |
| `xiong2023_sex_bias/PT-S2_fWTvsmWT.csv` | `59e9103e6210a43d6569166d8d3e40b3ec1e089e3caadcbe4244bd8a23c220a6` |
| `xiong2023_sex_bias/PT-S3_fWTvsmWT.csv` | `966286c7ae61b2f78db9b6f3aa4fe3d9c942d67527b0ca309f96ed8ded9cd9cd` |
| `census_pt_segments/mouse_pt_segment_pseudobulk_counts.csv.gz` | `d8bb110cd86c3d47c2d627de71edc17f5ad1634c6b9515e04162666adb0871bf` |
| `census_pt_segments/human_pt_segment_pseudobulk_counts.csv.gz` | `fcca3c628b94d19057d9c821748260900788b1f23d0365838b9efe95d75d5f47` |

**The builder** is `build_lists.py --scratch79 <scratch> --repo <repo>`. It stays in session scratch, because it reads
the 79 session's lit-full extraction, and it holds no absolute path. The committed CSVs are the source of truth.

## Why, and what our method can add

Each published comparison of human and mouse kidney is one of three things: a level difference inside one segment
(CellSpectra), a whole-tissue difference (Raghubar), or a disease response (Zhou, Hinze). Cheval 2012 measured S1 and S3
in both species but never compared the S1→S3 profile. Our data can add three things to a published claim:

1. **Reproduce** (level). Does the published direction hold in all 4 human × mouse specimen pairs?
   - This is descriptive. Per-gene probe efficiency differs between the species' panels, and the level t fails the
     specimen-swap calibration at every scope (notebook 79: 0.113–0.210 against 0.075).
2. **Place** (between segments). Does the difference change between S1, S2 and S3, and where is it largest?
   - A per-gene probe factor cancels here.
   - The segment F-test passes calibration (0.055), so this is the inferential part.
3. **Describe along the axis** (continuous, descriptive). This is the user's addition. In each species: where does a
   program switch, how sharply, does it peak inside a segment, and does it switch at the same place in human and mouse?

## Bounds from the record (they travel with every result)

**Level.** Probe-exposed; descriptive (above).

**Within-segment position is the weakest part of the record.** It is supported only where specimens agree.
- Notebook 72: the package coordinate is noise-dominated in 10 of 12 specimen × segment cells.
- Notebook 74, within-segment slope agreement between specimens beats the step nulls in mouse S1, S2 and S3 and in human
  S1 and S3. It does not in **human S2** (ρ 0.09; 632 units).
- Notebook 74, gene onset order in mouse S2 beyond the segment means: inconclusive (0.11). External onset
  replication: not met (0.30, interval touching 0).
- Rules that follow:
  - Continuous descriptions are made for **programs** (medians over many genes), not single genes, except the
    gene-level hypotheses (C2, D), which are labelled single-gene.
  - A position is stated only when both specimens of a species agree (section G).
  - Every position inside human S2 carries the flag "human S2: weak cell".
  - A landmark-coordinate check (LRS-Z, notebook 53's zone-landmark ratio, the best within-segment coordinate in
    notebook 72) runs on the matched units.

**Cross-species positions.** The species registration gap is 0.37 segment units (72, C8 not met). A human−mouse
difference in switch position below 0.37 segment units is therefore reported as "same place within registration error".

**Cohort.**
- Two male control mice and two cortex sections of one male human donor.
- Human S3 is cortical; mouse S3 includes the outer stripe.
- Everything is descriptive at n = 2 + 2.

## A · CellSpectra (Klötzer et al. 2025, Nat Genet; PMID 40775269)

### A1 · Are CellSpectra's segment-specific species genes segment-specific? (confirmatory; notebook 81)

**Why.**
- CellSpectra lists the top 99 species-difference genes per segment. Their blocks overlap heavily: 40 genes are in all
  three blocks of each direction.
- A gene listed in one segment only implies a difference specific to that segment. A gene listed in all three implies
  a PT-wide difference.
- A ranked top-99 list cannot separate "different here" from "ranked lower elsewhere". CellSpectra has markedly fewer
  human PTS3 cells, so S3 ranks are less powered. Our between-segment contrast tests it directly.

**Programs** (from 79's frozen blocks; already R1–R5 filtered):

| id | members | rule |
|---|---|---|
| U-H1, U-H2, U-H3 | 21, 5, 32 | human-high block of Sk, minus genes in the human-high blocks of the other two segments |
| U-M1, U-M2, U-M3 | 18, 8, 44 | the same for the mouse-high blocks |
| U-HC, U-MC (contrast) | 40, 40 | listed in all three blocks of that direction |

**Statistics.** Exactly notebook 79's:
- P_k, Q_k = d · P_k, pairs, relabelings;
- between-segment cameraPR on the segment F z (75), BH over the 8 programs;
- the localization class of 79 section 5.3.

**Readout.**
- **Supported (segment-specific):** U-xk is *localized to Sk*, or *graded* with its condition peak in Sk.
- **Not supported (PT-wide):** U-xk is *uniform within resolution*, or its peak lies in another segment with margin
  ≥ 0.25.
- **Otherwise:** not resolved.
- The core programs are expected to be uniform. A unique program that looks like the core is the "missed" reading: the
  segment-specific listing reflects ranking, not biology.
- Replication (level) is reported with 79's machinery as a descriptive column.

### A2 · Gene-level classes of every CellSpectra block member (descriptive; notebook 81)

For each member of S-H1…S-M3 (the claimed sign d, the listed segment k), with Q_j = d · P_j in the condition and in
each of the 4 pairs:

| Class | Rule |
|---|---|
| present in all segments | Q_j ≥ 0.25 for j = S1, S2, S3, in the condition and all 4 pairs |
| restricted to the listed segment | Q_k ≥ 0.25 in the condition and all pairs, and the condition Q_j < 0.25 for some other j |
| elsewhere only | Q_k < 0.25 (condition), and some other Q_j ≥ 0.25 in the condition and all pairs |
| opposite somewhere | some Q_j ≤ −0.25 in the condition and all 4 pairs |
| not reproduced | none of the above |

The rules are applied in this order. Counts are reported per block and per unique/core subset. Pathway labels are
attached to the "restricted" genes (section L).

### A3 · Do CellSpectra's conserved PT functions keep their place along the PT? (confirmatory; notebooks 81 and 82)

**Why.** CellSpectra reports PT functions as conserved across species: high cross-species coordination, a co-expression
property. It does not report where along the PT each function sits in each species. Coordination and placement are
different properties; this tests placement.

**Programs.** Mouse GO annotations with descendants (QuickGO), restricted to the universe, label genes removed.
- F-LCFA: long-chain fatty acid-CoA ligase activity, 10 genes;
- F-ARG: arginine metabolic process, 8;
- F-PLB: regulation of phospholipid biosynthetic process, 19;
- F-OAT: organic acid transport, 209;
- F-TCY: transcytosis, 18 (high leave-one-out coordination in their Supp Fig 10c);
- F-LDL: LDL particle binding, 13 (low in their Supp Fig 10c).

**Readouts per function.**
- **Within-species placement:**
  - the program median of each specimen's segment contrast (79 section 3.4);
  - the peak segment per specimen;
  - the section G descriptors per species.
- **Species × segment:** 79's localization class on the human − mouse difference (d = sign of the condition P at the
  function's largest |P|). Between-segment BH over the 6.
- **Class:**
  - *placement conserved*: the within-species peak segment is the same in all 4 specimens, and the localization class
    is *no difference* or *uniform within resolution*;
  - *placement differs*: the localization class is *localized* or *graded*, or the two species' peak segments differ
    in both specimens;
  - otherwise *not resolved*.
- **Notebook 82** repeats the within-species placement in the atlases (section B3).

### A4 · Not testable here (recorded)

- Fewer PTS3 cells in human: our human sections are cortex only, so the S3 fraction reflects sampling.
- Disease-specific coordination changes (their Supp Table 5): our human tissue is healthy.

## B · Cheval et al. 2012 (PLoS One; PMID 23056504): S1 and S3 microdissected in both species (notebook 82)

**Data.** SAGE tags per 10,000 per structure.
- One tag equals the column's smallest non-zero value, so counts are c = value / unit, with library N = 10,000 / unit.
- Eligible genes have c_S1 + c_S3 ≥ 10 in each species.
- Contrast d = log2((c_S3 + 0.5)/N_S3) − log2((c_S1 + 0.5)/N_S1), with Poisson SE² = (1/(c_S3 + 0.5) + 1/(c_S1 + 0.5)) / ln²2.
- *Zonated* in a species means |d| ≥ 0.5 and |d|/SE ≥ 2.
- Gene pairing follows the authors' table (`cheval2012_mouse_symbols.csv`).

### B1 · Their data, our question: is the S1→S3 contrast conserved, and does conservation rise with strength? (no data of ours)

**Why.** Cheval concluded that expression diverges: structures cluster by species before segment. They never compared
within-species contrasts, which cancel per-library and per-protocol factors. Our record (notebook 49, Lake and Census
check) says human and mouse agree mainly on their most strongly zonated genes, and human is flatter.

- **B1a.** Among genes zonated in mouse, the share with the same sign in human, by mouse |d| bins [0.5, 1), [1, 2),
  [2, ∞). Symmetrically, among genes zonated in human, by human |d| bins.
  - Binomial p against 0.5 per bin, and a Cochran–Armitage trend p across bins.
  - **Hypothesis:** agreement > 0.5 and rising with strength, in both selection directions.
- **B1b.** The slope through the origin of human d on mouse d (mouse-zonated genes), and of mouse d on human d
  (human-zonated genes), with a gene bootstrap (2,000; seed 81).
  - **Hypothesis:** the human-on-mouse slope is below the mouse-on-human slope (human flatter). Descriptive beside the
    atlas values of 0.12 and 0.62.

### B2 · Cheval against our data (notebook 82)

Genes: Cheval-eligible ∩ 71's universe, minus the 15 label/orientation genes. Our S3 − S1 per species is the mean of
the two specimens' pseudobulk r_S3 − r_S1 (79 section 3.2).
- **B2a (positive controls).** Spearman of Cheval mouse d with our mouse S3 − S1, and of Cheval human d with our human
  S3 − S1. The human one checks our cortex-section zonation against microdissected human S1 and S3.
- **B2b (confirmatory, one test).** Does our species × position interaction replicate in microdissected tissue?
  - x = Cheval (d_human − d_mouse) and y = our (human S3 − S1) − (mouse S3 − S1).
  - Spearman ρ, one-sided permutation p (10,000 gene permutations, seed 82).
  - **Beyond flattening:** the partial Spearman given both mouse contrasts (Cheval's and ours), with the same
    permutation, as in notebook 77's slope-free question. A flatter human makes both interactions ≈ −(1 − k) × mouse
    contrast, so the raw ρ alone could come from shared mouse zonation.
  - **Claims:** "replicates" when the raw ρ > 0 with p ≤ 0.05; "replicates beyond flattening" when the partial also
    has p ≤ 0.05.
- **B2c (79's transfer calls in Cheval).** Z-M1, Z-M3, Z-H1 and Z-H3; there is no S2 in Cheval.
  - Members ∩ eligible; ≥ 5 needed, else "not testable".
  - The own-segment sign: S3 programs use d; S1 programs use −d.
  - Per species: the median member contrast m, and the share s of members with a positive own-segment contrast.

  | Class | Rule |
  |---|---|
  | zonated | m ≥ 0.25, s ≥ 0.7 and sign-test p ≤ 0.05 |
  | reversed | m ≤ −0.25, s ≤ 0.3 and sign-test p ≤ 0.05 |
  | absent | \|m\| < 0.25 |
  | attenuated | otherwise |

  - A call **replicates** when Cheval's source-species class is *zonated* and its other-species class equals 79's
    transfer class. "Attenuated" and "absent" count as agreeing only with themselves.

### B3 · The same in single-nucleus atlases with donor replication (notebook 82)

**Donors** (77's helpers `species_calls.census_donors` and `lake_donors`):
- Census mouse: male donors at "3-6 month-old mature stage" and "12-month-old stage", segments S1, S2 and S3, at least
  50 cells per segment. This is notebook 74's six-donor set.
- Lake human: normal-disease cortex donors with at least 50 nuclei in each of S1, S2 and S3 (S2 optional). These are
  notebook 38's primary donors.
- Census human: 6 male donors, PCT (S1/S2) against S3.

**Per donor and program:** the median member contrast of the program's segment against the other segments (79
section 3.4). Census human uses S3 − PCT, sign-flipped for S1/S2 programs.

| Class (per species) | Rule |
|---|---|
| zonated | the donor median ≥ 0.25, ≥ 80% of donors ≥ 0.25, and Wilcoxon signed-rank p ≤ 0.05 |
| reversed | the mirror of zonated |
| absent | \|donor median\| < 0.25 |
| attenuated | otherwise |

**Programs.**
- 79's six Z programs. Positive control: the source-species atlas. Transfer: the other species.
- The A3 functions: the peak segment per species (Census mouse and Lake human) against our per-species peak.

**Readout:** the agreement of each 79 transfer class with the atlas class, per atlas.

## C · Raghubar et al. 2022 (Front Med; PMID 35872784) (notebook 81)

### C1 · Is their whole-cortex species difference PT-intrinsic, and is it uniform along the PT? (confirmatory)

**Programs.**
- **R-ENERGY** (mouse-high, 100): genes of their top-10 mouse-up GO terms named for respiration or oxidative
  phosphorylation, listed mouse-higher in their DE table, capped by |logFC|.
- **R-METAB** (mouse-high, 100): their other top-10 mouse-up terms, with energy genes removed.
- **R-STRUCT** (human-high, 100): their top-10 human-up terms.

**Readouts.**
- Replication at the PT scope with 79's family S machinery: level, calibrated empirically, descriptive. BH over 3.
- Localization class: between-segment BH over 3.
- Section G descriptors of the human − mouse program curve.
- For R-STRUCT, PT-intrinsic share: the share of members detected in at least 10% of units in both human sections,
  against the universe's share (descriptive). Their spots mix PT with glomeruli, vessels and interstitium.

**Hypotheses.**
- A metabolic-rate, age or platform explanation (theirs) predicts a uniform difference along the PT.
- A localized difference means they missed segment-specific biology.

### C2 · Their segment-relocation genes (single genes; descriptive)

F13b, Cela1, Acsm1, Anpep and S100g. Raghubar's Supp Material Table 8 places them in different segments or cell types
between species, from atlases. Predictions are in `genes_lit.csv`. Rule translations:
- "mouse: S1/S2 above S3": r_S3 − mean(r_S1, r_S2) ≤ −0.25 in both mice.
- "mouse: S3 above S1/S2": ≥ +0.25 in both mice.
- "human: low in PT": the centred human − mouse PT difference is ≤ −1 log2 in all 4 pairs (level; probe counts shown).
- "human: PT-expressed": detection in both human sections ≥ the universe median detection.
- Section G descriptors per species are labelled single-gene.

## D · Fernandes et al. 2026 (Acta Physiol; PMID 42494142): phosphate transporters (single genes; notebook 81)

Slc34a1, Slc34a3, Slc34a2, Slc20a1 and Slc20a2. Per species:
- detection;
- segment profile per specimen;
- section G descriptors.

**Readouts.**
- **Mouse Slc34a1 S1/S2 above S3** (their scRNA statement): the C2 rule.
- **Slc34a3 − Slc34a1 log ratio higher in human than in mouse** (their activity data): (human − mouse) of
  (r_Slc34a3 − r_Slc34a1) at the PT scope is ≥ 0 in all 4 pairs. This is level, probe-exposed, and both genes' probe
  counts are shown.
- **New information:** the S1→S3 profiles of SLC34A1 and SLC34A3 in human PT, which they lacked.

## E · Xiong et al. 2023 (Dev Cell; PMID 37673062): does mouse-only zonation follow male androgen programs? (confirmatory; notebook 81)

**Why.**
- Xiong maps mouse sex dimorphism mostly to PT, androgen-driven, with limited conservation in human.
- All our specimens are male, so sex is not confounded with species. A male human PT, however, need not carry the
  male-mouse androgen program.
- If mouse-only PT zonation is enriched for male-biased mouse genes, part of the human–mouse zonation difference is the
  mouse androgen program. Neither paper connects the two.

**Definitions** (our per-gene within-species segment contrasts c, 79 section 3.4, per specimen):
- **Mouse-zonated in Sk:** c ≥ 0.5 in both mice.
- **Conserved:** mouse-zonated, and human c ≥ 0.25 in both sections.
- **Mouse-only:** mouse-zonated, and human mean c < 0.25. The rest are excluded.
- **Male-biased in Sk:** in Xiong's PT-Sk table (F vs M), p_adj < 0.05 and log2FC ≤ −1. Female-biased is the mirror.

**Test.**
- Cochran–Mantel–Haenszel over quintiles of the gene's mean mouse rate, for OR(male-biased | mouse-only against
  conserved). One-sided p.
- **Primary:** pooled over segments, gene level (a gene counts as mouse-only if it is mouse-only in any segment and
  conserved in none). Per-segment ORs are descriptive.
- **Negative control:** female-biased OR, descriptive.
- Pathway labels (section L) for mouse-only × male-biased, mouse-only × unbiased, and conserved genes.

## F · Basit et al. 2019 (Drug Metab Dispos; PMID 31123036; abstract only) (single genes; descriptive; notebook 81)

- MDR1, OCT3 and OCTN1 are the only 3 of 19 transporters not more abundant in rodents.
- Testable here: Abcb1a and Slc22a4. Slc22a3 is not in the universe, and Abcb1b is not one-to-one.
- **Readout:** the centred human − mouse PT difference minus the median over the universe's Slc* and Abc* genes is
  ≥ 0 in all 4 pairs. This is level and probe-exposed.

## G · Continuous description of every program (the user's descriptive readout; notebook 81)

**Programs described:**
- A1, A3 and C1;
- 79's S and Z programs;
- the gene-level hypotheses (C2, D, F), labelled single-gene.

**Curves (package coordinate; no refit).**
- Input: 71's cached NB specimen curves (`curve|<specimen>`, natural-log rate on the 101-point grid).
- A member is usable in a specimen when it is converged, not separated, and detected in ≥ 5% of that specimen's
  analysed units.
- Each usable member's curve is centred by its mean over the grid and converted to log2.
- The **program curve** is the median over usable members, with ≥ 5 needed.
- Grid positions are in segment units from 72's pooled transitions (79's S_GRID).

**Descriptors per specimen curve.** Amplitude A = max − min. The shape class is the first that applies:

| Shape class | Rule |
|---|---|
| flat | A < 0.25 |
| interior peak or trough | `interior_extrema` with min_amplitude 0.25, min_prominence 0.125 and inner 0.2 |
| monotone rising or falling | \|Spearman(curve, grid)\| ≥ 0.9 |
| mixed | otherwise |

- For monotone curves: the **onset** (half-change crossing, in segment units, so 1.62 means 62% through S2) and the
  **width** (10–90% span), from `monotone_onsets`.
  - **Sharp switch:** width ≤ 0.5 segment units.
  - **Gradient:** width ≥ 1.0.
  - **Graded switch:** in between.
- For a peak or trough: its position.

**Species-level statement** (made only when both specimens agree).
- **Agreement:** the same shape class, and onsets (or peaks) within 0.25 segment units of each other. The stated
  position is their mean.
- If the specimens disagree: "specimens disagree", and no position.
- A human position in [1, 2): flag "human S2: weak cell".

**Cross-species.**
- When both species agree on a monotone class with the same direction, Δonset = human − mouse.
- |Δ| ≥ 0.37: "switches earlier or later in human".
- Otherwise: "same place within registration error".

**Species difference curve.** For the A1, A3 and C1 human − mouse curves, 79's program curves and `profile_summary`,
read with the same descriptors on the condition curve and the 4 pair curves.

**Landmark check** (matched set M of 72: package units mapped to old structures, LRS-Z from notebook 53).
- Unit-level program score: the mean over usable members of log2((count + 0.5)/library × the median library),
  centred per gene per specimen.
- `specimen_curves` (least squares, cubic B-spline) per specimen, on the package coordinate and on LRS-Z, each with 3
  interior knots at that coordinate's quartiles on set M.
- Each coordinate is converted to its own segment units (`logistic_transitions` of the package labels on that
  coordinate).
- **Agrees:** the same shape class and position within 0.25 segment units, per specimen.
- A species-level position is labelled "landmark-confirmed" when both specimens agree on both coordinates.

**What a descriptive finding looks like.** "In mouse, F-OAT rises sharply at S2 0.6 (width 0.3; both mice;
landmark-confirmed); in human it is a gradient (width 1.4); the species difference is localized to S3." Positions
without specimen agreement are not stated.

## L · Pathway labels (no tests)

A group of genes (A2 "restricted", the E1 groups) gets a label: the term of Hallmark 2020, KEGG 2019 Mouse or Reactome
2022 (71's library) containing the most of its genes (≥ 3; ties by library order, then name). The genes are listed
beside the label. No p-value is computed.

## Multiplicity and roles

| Family | Question | BH over |
|---|---|---|
| A1 | between-segment | 8 programs (confirmatory) |
| A3 | between-segment | 6 (confirmatory) |
| C1 | replication | 3 (confirmatory; level is descriptive at our calibration) |
| C1 | between-segment | 3 (confirmatory) |
| E1 | pooled CMH | 1 test |
| B2b | interaction | 1 test (plus the beyond-flattening partial) |

Class-rule readouts (B1, B2c, B3), gene-level hypotheses (C2, D, F) and section G are descriptive.

## What would make us say the method adds nothing here

- **A1:** every unique program is uniform or not resolved. Then the per-segment listings are PT-wide differences, and
  position adds nothing to CellSpectra beyond 79.
- **B2b:** the raw ρ fails. Then our species × position interaction does not replicate in microdissected tissue at SAGE
  depth, and every positional species claim must say so.
- **G:** no program reaches a specimen-agreed shape in either species. Then the continuous description adds nothing
  beyond the segment summaries.

## Outputs

- **Notebook 81** (`analysis/notebooks/81_pt_published_findings.ipynb`, logic `81.published_findings.1`) writes
  `results/pt_published_findings/tables/`:
  - `programs.csv`, `program_tests.csv`, `classes.csv`;
  - `gene_classes_cellspectra.csv`;
  - `gene_hypotheses.csv`;
  - `sex_bias.csv`;
  - `descriptors.csv` (per program × species × specimen) and `descriptors_species.csv`;
  - `landmark_check.csv`, `labels.csv`, `summary_by_paper.csv`;
  - `run_record.json`.
- **Notebook 82** (`analysis/notebooks/82_pt_published_findings_external.ipynb`, logic
  `82.published_findings_external.1`) writes `results/pt_published_findings_external/tables/`:
  - `cheval_genes.csv`, `cheval_strength.csv`, `cheval_vs_ours.csv`, `cheval_programs.csv`;
  - `atlas_donor_programs.csv`, `atlas_classes.csv`;
  - `run_record.json`.
- New reusable code goes in `pseudospace/published_findings.py`, with `tests/test_published_findings.py` (synthetic):
  - Cheval counts and contrasts;
  - the strength-bin agreement and trend;
  - the B2c class;
  - the gene class (A2);
  - the CMH test;
  - the continuous descriptors and species agreement (section G);
  - the pathway label.
- Notebook 79's functions are reused through `program_localization`.

## Notebook 81 results

`analysis/notebooks/81_pt_published_findings.ipynb` (logic `81.published_findings.1`) writes to
`results/pt_published_findings/`. It runs in about 1 min and fits no model.

**Guards.** All passed:
- the export contract and fingerprint `efafc9ac2e99`;
- 71's universe (11,055 genes), sets (1,513) and units (11,383);
- 71's cached model reproduces its T_spatial and z;
- the 8 frozen files match their sha256;
- 72's transitions reproduce;
- the pseudobulks and calibration equal notebook 79's.

**Calibration** (equal to 79's): every level t fails (0.113–0.210); the segment F z passes (0.055); the positional F z
passes (0.071).

### A1 · CellSpectra's segment-specific species genes

| Program | n | P S1 / S2 / S3 (log2) | Between q | Localization (peak; margin) | Readout |
|---|---|---|---|---|---|
| U-H1 human-high, S1 only | 21 | 2.50 / 2.36 / 1.45 | 7e-5 | graded (S1; 0.14) | supported |
| U-H2 human-high, S2 only | 5 | 2.06 / 2.28 / 3.02 | 0.017 | localized to S3 (0.74) | not supported: largest in S3 |
| U-H3 human-high, S3 only | 32 | 1.33 / 1.50 / 2.48 | 1e-8 | localized to S3 (0.97) | supported |
| U-HC human-high, all three | 40 | 2.97 / 3.25 / 3.73 | 0.018 | localized to S3 (0.48) | core, not uniform |
| U-M1 mouse-high, S1 only | 18 | −1.97 / −1.47 / −0.53 | 4e-6 | localized to S1 (0.51) | supported |
| U-M2 mouse-high, S2 only | 8 | −1.78 / −2.14 / −0.64 | 0.14 | not resolved (S2) | not resolved |
| U-M3 mouse-high, S3 only | 44 | −1.70 / −2.35 / −2.16 | 1e-8 | graded (S2; 0.18) | not resolved |
| U-MC mouse-high, all three | 40 | −1.48 / −1.54 / −1.39 | 0.013 | graded (S2; 0.06) | core, near uniform |

- **Supported:** 3 of the 6 segment-unique programs (U-H1, U-H3, U-M1) differ most in their own segment. The
  "adds nothing" rule is not met.
- **Not supported:** U-H2 (5 genes) differs most in S3, not S2. Every difference is present (≥ 0.25 log2) in all three
  segments. Segment specificity is a matter of degree (0.5–1 log2 on 1.5–3 log2), not presence.
- **The human-high core is not uniform.** It grows toward S3 (2.97 → 3.73 log2). Caveat: human S3 is cortical; mouse
  S3 includes the outer stripe.
- **Relabelings:** the relabeled ratios are ≤ 0.08 everywhere.
- **Level (descriptive):** human-high programs replicate; mouse-high programs hold direction only (79's relabeling-1
  pattern).

### A2 · The genes behind the CellSpectra blocks

Of 564 block members:

| Class | Genes | Share |
|---|---|---|
| present in all segments | 420 | 74% |
| restricted to the listed segment | 37 | 7% |
| opposite somewhere | 56 | 10% |
| not reproduced | 42 | |
| elsewhere only | 9 | |

- Unique members are restricted to their listed segment in 11% of cases; core members in 1%.
- **What our data add:** CellSpectra's per-segment species lists mostly describe one PT-wide difference. The segment
  listing is reproduced as segment-restricted for few genes:
  - human: Mme and Slc16a9 (S1/S2); Rora, Nebl, Stox2, Tinag, Itpr2 and Ppara (S3);
  - mouse: Dab2, Plcb1, Enpp2, Gldc, Maf and Dpyd (S1/S2); Osbpl8, Chrm3, Pipox and Kcnc3 (S3).
- **Mouse-high genes that are human-high in at least one segment** in all 4 pairs include Pde4d, Fgfr2, Abcc2, Ptprd,
  Nhs and Itch. This is a level reading, so probe-exposed.

### A3 · CellSpectra's conserved PT functions: placement

| Function | n | Peak segment (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1) | Species × segment | Between q | Readout |
|---|---|---|---|---|---|
| long-chain fatty acid-CoA ligase | 10 | S1, S2, S1, S1 | not resolved | 0.44 | not resolved |
| arginine metabolism | 8 | S1, S1, S1, S3 | localized to S3 | 7e-6 | placement differs |
| regulation of phospholipid biosynthesis | 19 | S2, S2, S3, S3 | not resolved | 0.058 | placement differs |
| organic acid transport | 209 | S2, S2, S3, S2 | graded (S2) | 6e-6 | placement differs |
| transcytosis | 18 | S3, S2, S3, S2 | not resolved | 0.28 | not resolved |
| LDL particle binding | 13 | S3, S3, S2, S2 | not resolved | 0.44 | placement differs |

- **No function meets "placement conserved".**
- Two show a calibrated species × segment change:
  - **Arginine metabolism:** the human − mouse difference is largest in S3. In mouse the program dips in S2 (−0.40
    and −0.27); in human it is flat.
  - **Organic acid transport:** in mouse it is highest in S2 in both mice; in human it is nearly flat.
- Within-species amplitudes are small (mostly |contrast| < 0.4 log2), so "placement" here means weak zonation.
- **What our data add:** a function can be conserved as a PT function (CellSpectra's coordination) while its axial
  placement differs between species. Both cases with a calibrated difference are mouse-zonated and flat in human.

### C1 · Raghubar's whole-cortex species difference, inside the PT

| Program | Q at PT (all 4 pairs) | P S1 / S2 / S3 | Level | Localization (between q) | Detected ≥ 10% in both human sections |
|---|---|---|---|---|---|
| R-ENERGY (mouse-high) | 1.30 | −1.34 / −1.58 / −1.38 | direction only | uniform within resolution (0.67) | 0.88 (universe 0.66) |
| R-METAB (mouse-high) | 4.18 | −3.90 / −3.88 / −3.48 | direction only | graded, S1 (0.0067) | 0.36 |
| R-STRUCT (human-high) | 3.05 | 3.06 / 3.04 / 3.19 | direction only | uniform within resolution (0.36) | 0.65 |

- All three published directions hold inside PT units, in all 4 pairs. They fail the relabeling specificity, as
  CellSpectra's mouse-high blocks do.
- **Energy production:** the mouse-high difference is the same size along the PT. This fits Raghubar's own reading (a
  global metabolic-rate or age difference) and equally a platform offset; it is not segment biology.
- **Structural genes:** they are PT-expressed at the universe rate (0.65 against 0.66), so their human-high signal is
  not only non-PT spots. It is uniform along the PT.
- **Metabolic genes:** largest in S1–S2, and only 36% are detected in human PT. Many are mouse-only metabolic genes.

### C2, D, F · Single-gene hypotheses

| Gene | Published | Ours |
|---|---|---|
| F13b | mouse PT S1/S2; human DCT | mouse S1/S2 > S3 (both mice); human PT detection 0–1% (−7.0 to −7.9 log2); consistent |
| Cela1 | mouse PT S1/S2; human DCT/LOH | mouse S1/S2 ≫ S3 (−4.2 log2); human PT detection 0%; consistent |
| S100g | mouse pCT; absent in human | mouse S1/S2 ≫ S3 (−3.6 to −4.1); human 0%; consistent |
| Acsm1 | mouse PT S3; absent in human | human absent (consistent); mouse **not** S3-high (+0.16, +0.20) |
| Anpep | mouse PT S3; human PT | human PT 98% (consistent); mouse **not** S3-high |
| Slc34a1 | mouse mostly S1/S2 | mouse S3 − S1/S2 −2.97 and −2.91 (consistent). **New: human almost flat (−0.24, −0.70)** |
| Slc34a3 | larger share in adult human | human PT 45% and 38% detection (mouse 24–25%); Slc34a3 − Slc34a1 higher in human by 2.9–3.2 log2 in all 4 pairs (level, equal probe counts); both S1/S2-biased in both species, weaker in human |
| Abcb1a (MDR1) | not lower in human | +4.3 to +4.9 log2 above the transporter median, all 4 pairs (level); human detection 57%/48% against mouse 4–7% |
| Slc22a4 (OCTN1) | not lower in human | not consistent (−0.1 to −0.6) |

- **The new finding of this block:** the phosphate transporter NaPi-IIa (SLC34A1) is S1/S2-zonated in mouse (about
  8-fold) and nearly uniform along the human PT. SLC34A3 is relatively higher in human. Together these give a
  positional picture of Fernandes's functional result: the human PT spreads phosphate uptake along the tubule and over
  two transporters.
- **Raghubar's segment relocations:** the S1/S2 ones (F13b, Cela1, S100g) are reproduced. The S3 ones (Acsm1, Anpep)
  are not S3-zonated in our mice.

### E1 · Mouse-only zonation and male androgen programs (Xiong)

- **Primary, pooled:** male-biased share 0.066 among mouse-only zonated genes (n 1,570) against 0.089 among conserved
  genes (n 529); OR_MH 0.77, p 0.92. **Not supported.**
- **Per segment:** OR 1.28 (S1), 0.97 (S2), 0.65 (S3). The female-biased control is null.
- **What this means:** the human–mouse zonation difference in our male cohort is not explained by the male-mouse
  androgen program. Mouse-only zonated, male-biased genes carry a Xenobiotic Metabolism label (Abcc2, Asl, Bcat1, Crot,
  Fmo3, Shmt2, Vnn1 and others).
- **Conserved zonated genes** are labelled SLC-mediated transport. This is consistent with 49's "transporters
  conserved".

### G · Continuous description

**What the first run showed.** Most programs' within-species curves had the same human trough at mid-S1 and peak at
S3 0.5–0.75. A post hoc check found that size-matched random gene sets do the same. The curves are now centred by each
specimen's median usable gene at each position (deviation 1 below). After that, size-matched random programs are mostly
flat:
- at sizes 20–50, human 79–99% and mouse 69–97% flat;
- at size 10, 32–35% flat;
- about 1–2% reach an agreed non-flat species statement.

This post-fix null comes from `diag_null.py` in session scratch.

**Agreed, program-specific (null match 0 of 200) and confirmed on LRS-Z: two statements, both in mouse.**
- **CellSpectra mouse-high S1 genes (S-M1):** a gradual fall along the mouse PT, half done at early S3 (S3 0.18;
  gradient). Flat in human.
- **CellSpectra mouse-high, S3-only genes (U-M3):** a gradual rise in mouse, half done in late S1 (S1 0.80; gradient).
  Flat in human.

**Agreed and program-specific, not confirmed on LRS-Z (descriptive only):**
- mouse: Z-M3 (Chen PTS3), an interior peak at S3 0.66; U-M1, a fall half done at the S2/S3 boundary;
- human: Z-M1 and Z-H1, peaks at S1 0.48; Z-H3 and F-PLB, troughs at S1 0.5–0.6; S-H3, a rise half done at S1 0.91.

The human mid-S1 positions sit where the package labels place 47–60% of the old human S2 units (issue 6). They may
mark that hidden boundary rather than biology.

**No sharp switch** (width ≤ 0.5 segment units) is confirmed in either species. The user's example kind of finding ("a
pathway switches sharply at a specific part of S2") is therefore not supported here. In particular, human S2 positions
are where the record is weakest (632 units; 74's slope agreement 0.09).

**Species difference curves.** Size-matched random programs give interior extrema in about half of draws, clustered
in S1 and at S3 0.5–0.75. Five difference curves exceed their null's 95th amplitude with ≥ 3 of 4 pairs:
- U-H1, peak S1;
- U-H3, peak S3;
- U-M1, trough S3;
- F-LCFA, peak S1;
- F-ARG, peak S3.

Their segments repeat the between-segment results; their within-segment positions are not interpretable.

**The curve figure** (`fig81b_program_curves`) shows the zonation transfer of notebook 79 directly. The human S1
program (Z-H1) falls the same way in both species. The mouse S1 and S3 programs (Z-M1, Z-M3) and the human S3 program
(Z-H3) move only in their own species.

### Readout against the "adds nothing" rules

- **A1:** 3 of 6 unique programs supported. Not met.
- **G:** 2 confirmed positions, both mouse gradients. Not met for mouse, met for human.
- **B2b:** notebook 82; it replicates (ρ 0.31, p 0.0001) and beyond flattening (partial ρ 0.19, p 0.0045). Not met.

### Recorded deviations (all after the first run, before writing this section)

1. **G within-species curves** are centred by the specimen's median usable gene at each grid point (`reference_curve`).
   The protocol centred each gene only by its own mean. Size-matched random programs reproduced the shapes most
   programs showed, so those shapes were shared, not program-specific.
2. **The landmark-check unit score** subtracts each unit's median usable gene. The 0.5 floor otherwise carries the unit
   library into every program score.
3. **Post hoc null for agreed statements:** 200 size-matched random programs per agreed statement, with
   "program-specific" meaning a null match ≤ 0.05.
4. **Post hoc null for difference curves:** 300 size-matched random programs, with "reportable" meaning amplitude
   above the 95th percentile and ≥ 3 of 4 pairs sharing the shape.
5. **E1 excludes the 15 label/orientation genes** (protocol R3's reason). The result did not change (OR 0.75 → 0.77).

**Tables:**
- `classes.csv`, `program_tests.csv`, `a3_placement.csv`, `calibration.csv`;
- `gene_classes_cellspectra{,_counts}.csv`, `gene_hypotheses.csv`;
- `sex_bias{,_genes}.csv`, `labels.csv`;
- `descriptors{,_species,_difference,_confirmed}.csv`, `landmark_check.csv`, `curves.npz`;
- `summary_by_paper.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.

**Figures:** `fig81a_cellspectra_unique_core`, `fig81b_program_curves`, `fig81c_gene_hypotheses`, `fig81d_sex_bias`
and `supp81_cellspectra_unique_curves`.

## Notebook 82 results

`analysis/notebooks/82_pt_published_findings_external.ipynb` (logic version `82.published_findings_external.1`) writes to
`results/pt_published_findings_external/`. It followed the protocol in `pt-published-findings.md` (commit 0a44dec);
section 14 of the notebook lists how the protocol was read.

**Run.**
- 19 s wall. No model is fitted and no stage cache is used.
- Every rule comes from `pseudospace/published_findings.py`. The only notebook-local rule is the Cheval modal-unit
  sensitivity (below).
- The test suite passes (555) in the activated project environment. Two `test_program_localization.py` tests need R limma and fail rather than skip when R cannot find it (an unactivated shell); this is recorded as an issue.

**Guards.** All passed:
- the export contract and the pass-1 fingerprint `efafc9ac2e99`; 71's universe (11,055) and units (11,383); the 15
  label/orientation genes equal the export's panels;
- the recomputed specimen × scope pseudobulks equal 79's `pseudobulk_rates.csv.gz` (max difference 4e-15);
- the program medians of our segment contrasts reproduce 79's source and transfer contrasts for the six Z programs;
- the program files (79's `programs_hvm`, this protocol's `programs_lit`) and the Cheval and Census count files match
  their frozen sha256; every member is in the universe and none is a label/orientation gene;
- Census mouse: 6 male donors at 3–6 or 12 months, 330–1,513 cells per segment;
- Lake: notebook 75's S1/S2/S3 normal-cortex cache, mapped to our symbols by 77's rule, equals the protocol's npz
  exactly on all 26 normal-donor S1/S3 groups and 11,053 shared genes (sums, totals and nuclei). The 7 donors with
  ≥ 50 nuclei in S1, S2 and S3 are notebook 38's primary donors (3 male), all "normal".

### Cheval 2012 data

| | Value |
|---|---|
| Table rows | 4,644 |
| Rows without a mouse symbol (left out) | 258 |
| Paired genes | 4,317 (69 listed in several rows) |
| Eligible (≥ 10 tags over S1 and S3 in each species) | 270 |
| Eligible and zonated in mouse / in human | 101 / 71 |
| Eligible ∩ 71's universe, minus label genes (B2) | 185 |

One tag, by the protocol's rule (the smallest non-zero value): human S1 0.18, human S3 0.21, mouse S1 0.10, mouse
S3 0.10 per 10,000. The most frequent non-zero value is 0.23, 0.22, 0.12 and 0.13. Only 12, 2, 4 and 1 values lie
below it. So the smallest value is not one tag in any column. The counts are inflated by 5–33%, which shifts the
eligibility line and the SE. The modal-unit sensitivity reruns every Cheval readout; it changes no call.

### B1 · Is the S1 → S3 contrast conserved in Cheval, and does conservation rise with strength? (no data of ours)

**B1a.** The share of zonated genes with the same sign in the other species.

| Selected in | [0.5, 1) | [1, 2) | [2, ∞) | Pooled | Trend p (rising) |
|---|---|---|---|---|---|
| mouse | 22/32 = 0.69 (p 0.025) | 26/42 = 0.62 (p 0.082) | 16/27 = 0.59 (p 0.22) | 0.63 (p 0.005) | 0.78 |
| human | 9/15 = 0.60 (p 0.30) | 25/41 = 0.61 (p 0.11) | 10/15 = 0.67 (p 0.15) | 0.62 (p 0.028) | 0.35 |

p is the one-sided binomial p against 0.5.

- **Agreement above 0.5:** yes, in both directions (pooled 0.63 and 0.62).
- **Rising with strength:** no, in either direction. In mouse-selected genes the share falls (0.69 → 0.59).
- The hypothesis (above 0.5 and rising) is therefore half supported.

**B1b.** Slopes through the origin (2,000-draw gene bootstrap, seed 81):

| Slope | Genes | Estimate (95% interval) | Atlas value |
|---|---|---|---|
| human d on mouse d (mouse-zonated) | 101 | 0.09 (−0.12 to 0.31) | 0.12 |
| mouse d on human d (human-zonated) | 71 | 0.11 (−0.24 to 0.41) | 0.62 |

- The human-on-mouse slope is below the mouse-on-human slope, but only by 0.02, and both intervals include 0 and
  overlap.
- Both slopes are near 0.1. The mouse-on-human slope is far below the atlas value of 0.62. At SAGE depth, the
  magnitude of one species' contrast predicts little of the other's in either direction. The data cannot separate
  "human flatter" from noise. Under the modal unit both slopes are −0.01.

### B2 · Cheval against our data

185 genes (Cheval-eligible ∩ 71's universe, minus the label/orientation genes).

| Readout | ρ | p | Claim |
|---|---|---|---|
| B2a mouse: Cheval d vs our mouse S3 − S1 | 0.46 | 6e-11 (Spearman) | positive control holds |
| B2a human: Cheval d vs our human S3 − S1 | 0.14 | 0.057 (Spearman) | weak |
| **B2b interaction (raw)** | **0.31** | **0.0001** (one-sided permutation, 10,000, seed 82) | **replicates** |
| **B2b beyond flattening** (partial, given Cheval's and our mouse contrasts) | **0.19** | **0.0045** | **replicates beyond flattening** |

- **Mouse positive control.** Our mouse S3 − S1 agrees with microdissected mouse S1 and S3 (ρ 0.46).
- **Human positive control.** Our cortex-section human S3 − S1 agrees only weakly with microdissected human S1 and S3
  (ρ 0.14, p 0.057). Human zonation is the weaker measurement in both data sets.
- **The confirmatory test.** Our species × position interaction replicates in microdissected tissue (ρ 0.31), and the
  replication survives the partial on both mouse contrasts (ρ 0.19). So it is not only shared mouse zonation times a
  flatter human.
- Once both mouse contrasts are partialled out, what remains of each interaction is mostly its human contrast. So
  "beyond flattening" here means that the human-specific part of the S1 → S3 contrast agrees between Cheval's human
  tissue and our human sections.
- **Context.** The human-on-mouse slope through the origin on these 185 genes is 0.11 in Cheval and 0.14 in our data.
  Both data sets show a flatter human.
- The abandonment rule for B2b (the raw ρ fails) is not met.

### B2c · Notebook 79's transfer calls in Cheval

| Program | Members eligible | Mouse: class (m, s, p) | Human: class (m, s, p) | 79 transfer | Call |
|---|---|---|---|---|---|
| Z-M1 (Chen PTS1) | 7 of 100 | zonated (1.69, 1.00, 0.008) | attenuated (1.41, 0.71, 0.23) | attenuated | **replicates** |
| Z-M3 (Chen PTS3) | 4 of 68 | not testable | not testable | absent | not testable |
| Z-H1 (Hinze S1) | 3 of 39 | not testable | not testable | zonated | not testable |
| Z-H3 (Hinze S3) | 0 of 61 | not testable | not testable | absent | not testable |

m is the median own-segment member contrast (log2), s the share of positive members, p the sign-test p.

- Only Z-M1 has the 5 eligible members needed. Its 7 eligible members are Gpx3, Slc4a4, Miox, Ass1, Slc39a5, Dcxr
  and Lgmn.
- Z-M1 is zonated in Cheval's mouse (the positive control) and attenuated in Cheval's human, as 79 found. Its call
  replicates.
- Caution: the human median is large (1.41). "Attenuated" here comes from the sign test (5 of 7 members positive,
  p 0.23), not from a small median. In Cheval's human the S1 program is inconsistent across members rather than
  smaller.
- The other three calls cannot be tested at SAGE depth.

### B3 · Notebook 79's zonation calls in the atlases (donors as replicates)

Donors: Census mouse 6 (male, 3–6 or 12 months); Lake human 7 (normal cortex, ≥ 50 nuclei in S1, S2 and S3); Census
human 6 (≥ 50 cells in both labels; 2 male, 4 female; S3 − PCT).

| Program | Census mouse | Lake human | Census human | 79 (source; transfer) |
|---|---|---|---|---|
| Z-M1 (S1, mouse) | **zonated** 2.05 (6/6) | absent 0.08 | absent 0.24 (3/6) | zonated; attenuated |
| Z-M2 (S2, mouse) | **zonated** 0.48 (6/6) | absent 0.02 | absent −0.25 | zonated; absent |
| Z-M3 (S3, mouse) | **zonated** 0.85 (6/6) | absent 0.09 | absent −0.00 | zonated; absent |
| Z-H1 (S1, human) | zonated 1.17 (6/6) | **zonated** 1.39 (7/7) | **zonated** 2.13 (6/6) | zonated; zonated |
| Z-H2 (S2, human; 6 genes) | zonated 0.59 (5/6) | **zonated** 0.57 (6/7) | **zonated** 1.31 (6/6) | zonated; absent |
| Z-H3 (S3, human) | absent 0.04 | **zonated** 0.97 (7/7) | **zonated** 1.35 (6/6) | zonated; absent |

Bold is the source-species atlas (the positive control). Values are the donor median of the own-segment contrast
(log2), with the donors ≥ 0.25. Every zonated class has Wilcoxon p 0.008–0.016.

- **Positive controls.** All 12 source-species cells are zonated, in every atlas.
- **Transfer: agreement with 79's class is 2 of 3 in each atlas** (6 of 9 in all).
  - Agreeing: Z-H1 conserved in mouse (zonated); Z-H3 absent in mouse; Z-M2 and Z-M3 absent in both human atlases.
  - **Z-M1** (mouse S1 program in human): absent in Lake (0.08) and in Census human (0.24, just below the 0.25 line).
    79 called it attenuated (our two human sections 0.25 and 0.24). Both atlases put it at or below the border our
    sections sit on.
  - **Z-H2** (human S2 program in mouse, 6 genes): zonated in Census mouse (0.59, 5 of 6 donors). 79 called it absent
    because our two mice disagree (−0.02 and 0.44).
- The atlases confirm the main transfer pattern of 79. The mouse S1, S2 and S3 programs are not zonated in human. The
  human S1 program is zonated in mouse. The human S3 program is not zonated in mouse.

### A3 (atlas half) · Do CellSpectra's conserved PT functions peak in the same segment?

Program median of the segment contrast (log2). "Range" is max − min over S1–S3.

| Function | Species | Ours S1 / S2 / S3 (range) | Our peak (per specimen) | Atlas S1 / S2 / S3 (range) | Atlas peak (donors) | Same peak |
|---|---|---|---|---|---|---|
| F-LCFA | mouse | 0.08 / 0.08 / −0.11 (0.19) | S2 (S1, S2) | 0.07 / 0.21 / −0.21 (0.42) | S2 (4 of 6) | yes |
| F-LCFA | human | 0.06 / −0.05 / −0.17 (0.22) | S1 (S1, S1) | 0.13 / −0.01 / −0.07 (0.20) | S1 (6 of 7) | yes |
| F-ARG | mouse | −0.03 / −0.34 / −0.06 (0.31) | S1 (S1, S1) | −0.11 / −0.27 / 0.04 (0.31) | S3 (5 of 6) | no |
| F-ARG | human | 0.07 / 0.03 / −0.01 (0.09) | S1 (S1, S3) | −0.24 / 0.05 / 0.14 (0.38) | S3 (5 of 7) | no |
| F-PLB | mouse | −0.15 / 0.09 / −0.10 (0.24) | S2 (S2, S2) | 0.03 / 0.05 / −0.21 (0.26) | S2 (5 of 6) | yes |
| F-PLB | human | −0.17 / −0.05 / 0.18 (0.35) | S3 (S3, S3) | 0.10 / −0.14 / −0.05 (0.24) | S1 (6 of 7) | no |
| F-OAT | mouse | −0.11 / 0.07 / 0.01 (0.17) | S2 (S2, S2) | −0.05 / 0.04 / 0.01 (0.09) | S2 (5 of 6) | yes |
| F-OAT | human | 0.01 / 0.04 / −0.01 (0.05) | S2 (S3, S2) | −0.01 / −0.00 / 0.01 (0.02) | S3 (3 of 7) | no |
| F-TCY | mouse | −0.02 / 0.03 / 0.02 (0.05) | S2 (S3, S2) | −0.05 / −0.06 / 0.05 (0.10) | S3 (5 of 6) | no |
| F-TCY | human | 0.02 / −0.05 / −0.04 (0.07) | S1 (S3, S2) | −0.06 / −0.04 / 0.09 (0.15) | S3 (6 of 7) | no |
| F-LDL | mouse | 0.02 / −0.34 / 0.26 (0.60) | S3 (S3, S3) | 0.06 / −0.11 / 0.03 (0.17) | S1 (3 of 6) | no |
| F-LDL | human | −0.08 / 0.19 / −0.32 (0.51) | S2 (S2, S2) | 0.02 / −0.23 / 0.04 (0.27) | S3 (3 of 7) | no |

- The atlas peak equals our peak in **4 of 12** function × species cells: F-LCFA in both species, and F-PLB and F-OAT
  in mouse.
- Most of these functions barely change along the PT. In 15 of the 24 profiles the range is below 0.25 log2. F-OAT
  (209 genes) and F-TCY are flat in every data set, so their peaks are not interpretable.
- Where both data sets show a range ≥ 0.25, the peaks disagree (F-ARG mouse, F-LDL human). F-LDL is the clearest
  case: our data show it in S3 in mice and in S2 in human, with both specimens agreeing in each species (ranges 0.60
  and 0.51). The atlases show at most 0.27, with no consistent peak (3 of 6 and 3 of 7 donors).
- Between species: the atlas peaks are the same for F-ARG and F-TCY (S3 in both). Our peaks are the same for F-ARG
  (S1) and F-OAT (S2).
- So the atlases do not confirm most of our within-species placements of these functions. The functions that
  CellSpectra calls conserved are mostly flat along the PT in both kinds of data.

### Sensitivities

| Sensitivity | Change |
|---|---|
| Cheval modal unit (1 tag = the column's most frequent non-zero value; 236 eligible genes, 157 in B2) | No call changes. B1a pooled 0.61 and 0.64 (p 0.031 and 0.027); trend p 0.81 and 0.80. B1b slopes −0.01 and −0.01. B2a 0.44 and 0.11 (p 0.18). B2b 0.27 (p 0.0005) and partial 0.17 (p 0.019). Z-M1 still replicates |
| Census human, the 6 male donors with no cell floor (S3 cells 2–287) | Z-M1 attenuated (0.32, 4 of 6) instead of absent, so it agrees with 79. Z-H3's positive control drops to attenuated (1.17, 4 of 6; p 0.31). Other classes unchanged |
| Lake, the npz as given (S1 and S3 only; 12 normal donors with ≥ 50 nuclei in each) | Same classes as the primary Lake arm for Z-M1 (absent, 0.10), Z-M3, Z-H1 and Z-H3. S2 programs not testable |

### Readout

- **B1 (Cheval alone).** Human and mouse S1 → S3 contrasts agree in sign more often than chance (0.62–0.63), but the
  agreement does not rise with gradient strength. Neither slope separates from 0 at SAGE depth.
- **B2 (confirmatory).** Our species × position interaction replicates in microdissected tissue (ρ 0.31, p 0.0001),
  and beyond flattening (partial ρ 0.19, p 0.0045). The mouse positive control is strong (0.46). The human one is weak
  (0.14), which bounds any human-specific positional claim.
- **B2c.** One of four 79 transfer calls is testable in Cheval (Z-M1), and it replicates.
- **B3.** Every published zonation program is zonated in its source-species atlas. 79's transfer class matches the
  atlas class in 6 of 9 transfer cells. Both mismatches are at the border of our own data: Z-M1 in human (our sections
  0.25 and 0.24) and Z-H2 in mouse (our mice disagree).
- **A3 (atlas half).** Atlas and our peaks agree in 4 of 12 cells. Most CellSpectra conserved functions are flat along
  the PT, so "placement" is not resolvable for them in either data set.

### How the protocol was read

- **Census human donors.** The protocol names "6 male donors". The file has 6 male donors, but only 2 have ≥ 50
  cells in both labels; the others have 2–28 S3 cells. The record's "Census human 6" (pt-gene-model-nb.md,
  pt-package-pseudospace.md) is the 6 donors with ≥ 50 cells in both labels (2 male, 4 female), which 77's
  `census_donors` rule gives. That set is primary. The 6 male donors with no floor are a sensitivity. The 2 male
  donors with ≥ 50 cells are below `donor_class`'s minimum of 3.
- **Lake.** The protocol's npz holds only S1 and S3 (groups `donor|PT-segment`, all diseases, no sex). The protocol's
  donors need ≥ 50 nuclei in S2 as well (notebook 38's primary donors). So the primary arm uses notebook 75's
  S1/S2/S3 normal-cortex cache of the same Lake pseudobulk, which equals the npz exactly on every shared group and
  gene. "S2 optional" is read as: S2 programs are tested where the atlas has S2. The npz-only reading is a
  sensitivity.
- **Cheval pairing.** Rows without a mouse symbol in the authors' table are left out, B1 included. Within a gene and
  species, rows with the same 8-structure values are one tag listed twice and count once; distinct rows are summed.
  The unit (one tag) is each column's smallest non-zero value over all 4,644 rows.
- **B1a.** "Agreement > 0.5": every non-empty bin above 0.5 and the pooled one-sided binomial p ≤ 0.05. "Rising": the
  one-sided Cochran–Armitage p ≤ 0.05.
- **B1b.** It is read on point estimates; the intervals are reported.
- **B2.** Our S3 − S1 is centred per specimen by its median over the universe (79 section 3.2). This does not change
  any Spearman statistic. The partial uses rank residuals on both mouse contrasts with residual permutation
  (`permutation_spearman`).
- **B3.**
  - Atlas contrasts are centred per donor by their median over the universe genes measured in that atlas (79 section
    3.4).
  - Members are 79's members measured in the atlas, with no detection filter.
  - Agreement is equality of classes; attenuated and absent agree only with themselves.
- **A3.** Our peak is the argmax of the mean of the two specimens' program medians. The atlas peak is the argmax of
  the per-segment donor median. Census human (two labels) is not used for A3. The range is shown as a descriptive size
  only; the protocol fixes no peak threshold.

**Tables:** `cheval_genes.csv`, `cheval_strength.csv`, `cheval_vs_ours.csv`, `cheval_programs.csv`,
`atlas_donor_programs.csv`, `atlas_classes.csv` (B3 and A3 rows), `how_the_protocol_was_read.csv`; `run_record.json`.

**Figures:** `fig82a_cheval_strength` (B1a by strength bin, both selection directions), `fig82b_cheval_vs_ours` (B2a
controls, B2b raw and beyond flattening), `fig82c_program_classes` (Z programs × data set, classes; A3 functions'
segment contrasts and peaks).

## Combined readout by paper (notebooks 81 and 82)

| Paper | Published | What our method adds, or finds they missed | Support |
|---|---|---|---|
| CellSpectra 2025 | per-segment human vs rodent DE lists | most listed genes differ in all three segments (74%); segment-only listings are segment-largest for human S1, human S3 and mouse S1, not human S2 | between-segment test (calibrated) |
| CellSpectra 2025 | conserved PT functions | conserved as functions, but mostly flat along the PT; arginine metabolism and organic acid transport are mouse-zonated, flat in human; atlas peaks agree in 4 of 12 cells | calibrated between-segment, but small amplitudes |
| Cheval 2012 | S1/S3 expression diverges between species | in their own data, S3 − S1 contrasts agree in sign (0.62–0.63), without a strength trend; our species × position interaction replicates in their microdissected tissue, beyond flattening | B2b confirmatory |
| Raghubar 2022 | mouse cortex higher in energy and metabolism, human in structure | holds inside PT units and is uniform along the PT; F13b, Cela1 and S100g relocations reproduced; Acsm1 and Anpep not S3 in our mice | level (descriptive) and between-segment |
| Fernandes 2026 | human phosphate uptake relies on SLC34A3 more than mouse | SLC34A1 is S1/S2-zonated in mouse (about 8-fold) but nearly flat in human; SLC34A3 relatively higher in human | single genes, descriptive |
| Basit 2019 | MDR1, OCT3, OCTN1 not lower in human | MDR1 (Abcb1a) much higher in human; OCTN1 not | single genes, level |
| Xiong 2023 | mouse PT sex programs, limited in human | mouse-only zonation is not enriched for male-biased genes (OR 0.77) | E1 (negative) |
| Chen 2021, Hinze 2022 (via 79) | segment programs per species | the transfer pattern (mouse S1–S3 programs absent in human; human S1 conserved in mouse; human S3 absent in mouse) replicates in Lake and Census in 6 of 9 cells; Z-M1 in Cheval | atlas donors |

**Continuous description.** Two gradual mouse gradients are specimen-agreed, program-specific and landmark-confirmed. No
sharp switch and no human position is confirmed.

**Bounds.**
- The human positive control against microdissected tissue is weak (ρ 0.14).
- The human-specific part of our zonation is therefore the least certain measurement. Even so, it is the part that
  carries B2b beyond flattening.
