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
