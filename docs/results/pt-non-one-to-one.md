# Genes outside one-to-one orthology along the elastic curve (notebook 86)

Frozen on 2026-10-13, before any curve, zonation call, family test or AKI statistic of notebook 86 was computed. The
user asked: "Analyze the non 1:1 orthologs with new pseudospace. Test as many things as you can think of to find
interesting results. Look at the literature to see questions we could ask."

**What I had seen before freezing.**
- Notebook 62's results on the retired scFates coordinate (summarised below), and notebooks 64 and 68's family results.
- 73's paralog-swap verdict on the package units (M8).
- Notebook 85's results, including the thin-limb units at the end of mouse S3 (issue 27).
- A data check without any fit (2026-10-13):
  - stage 01's per-sample unit files hold each species' full panel (mouse 19,059 genes, human 18,132). Joined to the
    export by `polygon_id`, every export unit is found and every shared gene's counts are identical;
  - the export's `primary_identity` per unit, and its spread along the mouse curve: almost every unit beyond mouse S3
    0.25 is Thin_Limb;
  - after the common-support trim, the Proximal_Tubule-identity units are Ctrl1A2 1,737, Ctrl1A4 1,844, HUK1_COR1
    2,625 and HUK1_MED1 2,636;
  - on those units, the PT-expressed genes (detection ≥ 2% in both specimens) by notebook 62's class and probe status.
    Mouse, probe-clean: 1:1 10,110; 1:many 172; many:1 180; many:many 148; no consensus 119; no ortholog 11. Human,
    probe-clean: 1:1 10,643; 1:many 178; many:1 138; many:many 158; no consensus 195; no ortholog 116.
- The literature agent's report (section L), written without access to our data or results.
- After the report: whether its named genes are on our panels, with their probe counts and 62's class (no fit). Mouse
  Gsta1, Nat8f6, Nat8f7 and Apol9b and human UGT2B15 and AKR1C1 are absent from the panels.

## What notebook 62 found, and why it is redone

Notebook 62 asked these questions on SCF13, the retired scFates coordinate:
- non-1:1 genes were 11% of strongly zonated mouse genes against 6% of expressed ones (odds ratio 2.1), not enriched in
  human (0.86);
- they were mostly rodent-expanded drug-handling families, partly a male-biased program (27% against 16%, post hoc);
- paralogs co-positioned (peak separation 0.23 against 0.36 for matched random pairs), and segment means found 16 of
  the 17 separations;
- family units differed between species no more often than 1:1 genes (6 of 20 against 7.3).

72's C14 and 73's decision retired 62's coordinate results (`pt-package-species-calls.md`, section 4). The user's rule
since then: the package's elastic curve is the only coordinate. Notebook 86 asks 62's questions again on that curve,
adds the questions 62 could not ask, and takes its hypotheses from the literature (section L).

## Labels

- Notebook: `analysis/notebooks/86_pt_non_one_to_one_orthologs.ipynb`
- Logic: `86.non_one_to_one.1`
- Module: `pseudospace/non_one_to_one.py`, tested in `tests/test_non_one_to_one.py`
- Results: `results/pt_non_one_to_one_package/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`)

## Inputs (read only)

| Source | What 86 reads |
|---|---|
| 70 | the human vs mouse export (`hvm/pt_units.h5ad`), the AKI export (`aki/pt_units.h5ad`), stage 01's unit files (`runs/hvm/01_units/<sample>_units.h5ad`) and the AKI QC units (`runs/aki/02_qc/units_qc.h5ad`) |
| 71 | `hvm_inputs` (as 71, 81, 84 and 85); its universe, knots and grid; its cached species and relabeling models, denominators and gene statistics (human vs mouse) |
| 72 | `transitions.csv` (mouse units) |
| 76 | its inputs (rebuilt by its own call), cached gene models, denominators, gene statistics, segment cuts |
| 62 | the coordinate-free columns of `gene_inventory_{mouse,human}.csv` only: gene, gene_id, component, n_human, n_mouse, ortholog_class, class_support2/3/5, paper_partner, probe_status, offtarget_genes, probe_clean, probe_clean_strict |
| External | Xiong 2023 (female vs male mouse PT S1/S2/S3); Chen 2021 (male PTS1–3, GSE150338) and Chen 2023 (male and female PTS1–3, GSE212213); Ensembl 98 within-species paralogs (new download, `data/external/ensembl_compara_paralogs`) and human–mouse homologies; Ensembl 98 mouse and r110 human cDNA headers (gene positions) |

62's old-coordinate columns (detection, q, amplitude, peaks, zonated, Z, segment means, flags) are never read.

**Guards (not cached).**
- G1: `hvm_inputs` reproduces 71's universe, units and gene sets; 71's cached knots and grid are its stored ones.
- G2: 72's transitions are reproduced (as 85's G3).
- G3: every export unit is found in its stage-01 file by `polygon_id`, and the shared genes' counts are identical
  (mouse symbols for mouse; `human_symbol` for human).
- G4: 62's inventories list exactly the stage-01 panels, in order.
- G5: 76's inputs are rebuilt by its own construction; its cached condition model's T_spatial is reproduced from its
  table (rtol 1e-12).
- G6: the sha256 of every external file is recorded; the Ensembl paralog files must match the values in
  `data/external/ensembl_compara_paralogs/SOURCE.txt`'s download (recorded in the run record).

## Coordinate, units and genes

- **Coordinate:** the package's elastic curve, as `hvm_inputs`' support-scaled position (and `aki_inputs`' for AKI).
- **Mouse units:** the mean of the two mice's own logistic label transitions (85's S0), piecewise linear.
- **Human:** curve position only. The human transitions cross (issue 6), so human segment units are undefined. A human
  gene's segment is read from label pseudobulks (W5).
- **Units.**
  - U_PT: the common-support units whose `primary_identity` is Proximal_Tubule. This removes the thin-limb units at the
    end of mouse S3 (issue 27) and the DCT, TAL and CNT units that carry neighbour signal (issue 24). Primary for every
    within-species question.
  - U_all: every common-support unit (71's units). Primary for the human vs mouse family test (H7), so that 71's 1:1
    statistics are its comparator; a sensitivity elsewhere.
- **Library:** `hvm_inputs`' library (the shared orthologs, mt-free), the offset of 71 and 85. Sensitivity (H1): the
  native-panel total without mt genes.
- **Genes:** each species' full panel minus mt genes (mouse `mt-`, human `MT-`). PT-expressed: detection ≥ 2% in both
  specimens on U_PT (62's rule).
- **Classes:** 62's HCOP orthogroups at support ≥ 3 (`class_support3`): 1:1, 1:many, many:1, many:many,
  no_consensus, no_ortholog. "Non-1:1" means 1:many, many:1 or many:many; no_consensus and no_ortholog are reported
  apart. The support ≥ 2 and ≥ 5 classes are sensitivities.
- **Probes:** a gene-level claim needs a probe-clean gene (62's clean or near). Family units use family-clean members
  (`non_one_to_one.family_clean`): probe-clean, or ambiguous only through members of the same orthogroup.

## W · Within-species curves (both species, U_PT)

- **W1 · Fits.** Per species, one ML dispersion per gene (design: specimen intercepts and the shared basis). Per
  specimen, an NB GLM with 71's 6-df cubic B-spline (71's knots), evaluated on 71's 101-point grid, NaN outside the
  specimen's 1st–99th percentile range. φ = log2(CP10k curve + 0.05) (62's pseudocount).
- **W2 · Typical shape.** Per specimen, the median over usable genes (detection ≥ 5% in that specimen, converged, not
  separated) of their mean-centred φ. 81 and 85 showed that raw within-species curves share a shape that random gene
  sets reproduce.
- **W3 · Test.** Per specimen, a likelihood-ratio test of the spline against the typical shape: both models carry the
  offset log(library) + ln 2 · shape(position), and the null has an intercept only (6 df, chi-square). BH per specimen
  over the PT-expressed genes.
- **W4 · Calls** (62's rule on the centred curves c = φ − mean φ − shape, log2):
  - zonated: q < 0.01 in both specimens, amplitude (max − min of c) ≥ 1 in both, Pearson r of the two specimens'
    c ≥ 0.5, converged and not separated in both;
  - strong: zonated with amplitude ≥ 2 in both; Z = the smaller amplitude;
  - peak and trough: argmax and argmin of the specimens' mean c (mouse units for mouse, curve position for human).
  - Mouse peak segment: S1 for a peak < 1, S2 for < 2, otherwise S3.
- **W5 · Label view.** Per specimen and label, log2((Σ counts + 0.5) / Σ library) on U_PT, centred by the median
  PT-expressed gene per label (`program_localization.segment_pseudobulk`, `segment_profile`'s centring). A gene's
  label peak is the label with the largest mean of its two specimens. This is the segment-resolution comparator, and
  the human peak segment.

## L · What the literature says (the agent's report)

A separate agent searched the literature without access to our data or results
(`non-one-to-one/literature_report.md`; every PMID checked against PubMed). What it changes here:
- **A positive control is needed first** (its Q1). Ancient 1:1 paralogs split the PT: SGLT2 (Slc5a2) in S1/S2 against
  SGLT1 (Slc5a1) later (human IHC, PMID 25304002), b(0,+)AT (Slc7a9) in S1 against rBAT (Slc3a1) in S3 (PMID 26739563).
  If the curve does not separate them, a "paralogs co-position" result cannot be read. This is H0.
- **Cross-species agreement of 1:1 genes is only partial** (liver: 111 of 313 zonated 1:1 orthologs, PMID 30297808).
  Every cross-species paralog statement is made against the 1:1 baseline (H8's B).
- **Duplicate-gene theory gives competing predictions:** dosage sharing (summed duplicates ≈ singleton; PMIDs
  20708291, 27199432) against divergence after duplication (15466287, 28743766); co-regulation of tandem duplicates
  (27199432, 27634932) against neighbourhood effects that are not specific to paralogs (18273387); expression
  divergence growing with duplicate age (12446139, 12840042, 26173681). No mammalian paper analyses paralogs along one
  tissue axis.
- **Sex:** mouse PT sex differences sit in S2 and S3, S1 is nearly flat (PMIDs 36758122, 37673062, 31689386). Whether
  sex-biased kidney genes are enriched in rodent-expanded families has never been tested; the one gene-age study found
  mouse kidney sex-biased genes are mostly old genes with 1:1 orthologs (PMID 37917687). H3d tests it.
- **Probe panels:** in our panels mouse Gsta1, Nat8f6, Nat8f7 and Apol9b, and human UGT2B15 and AKR1C1, have no
  included probe (structural zeros, confirmed in the stage-01 panels). Family sums are the defensible unit; human genes
  with one probe may look flat for lack of sensitivity (H1b's probe sensitivity).
- **Named claims** (H12): `non-one-to-one/literature_claims.csv` (60 rows: segment peaks, paralog orders, presence,
  flatness; mostly Chen 2023's deposited male microdissection and primary ISH/IHC papers). AKI directions:
  `non-one-to-one/literature_aki.csv` (19 rows).
- **Thin limb:** late-S3 "peaks" can be thin-limb composition (its Q16). U_PT removes those units; H12 checks DTL marker
  genes.

## Hypotheses and tests

Each test is fixed here, with its direction. "Two-sided" tests have competing predictions in the literature. The
primary p-values listed in the multiplicity section are also reported after BH.

### H0 · Positive control: the curve separates ancient 1:1 paralog pairs

- **Pairs:** Slc5a2 before Slc5a1, and Slc7a9 before Slc3a1 (mouse; claims C01–C02); SLC5A2 in S1/S2 and SLC5A1 in S3
  (human label view; C45–C46).
- **Passes** (mouse) when both pairs are separated by H4b's rule (both specimens, ≥ 0.5 units, the stated order).
- If it fails, H4–H6 and H8's positional statements are reported as "uninterpretable on this curve".

### H1 · Mouse PT zonation is enriched for genes outside 1:1 orthology

- **Prediction** (62; rodent expansions of xenobiotic families, section L): non-1:1 genes are over-represented among
  strongly zonated PT-expressed probe-clean mouse genes. Not in human.
- **Test.** Mantel–Haenszel odds ratio, non-1:1 against 1:1, of being strong, over deciles of log2 mean CP10k (the
  species' U_PT pseudobulk, specimens weighted equally). One-sided (OR > 1) for mouse, two-sided for human.
- **Holds** when the mouse OR ≥ 1.5 with p < 0.05.
- **Reported beside it:** each class separately; zonated in place of strong; support ≥ 2 and ≥ 5; the strict probe rule;
  U_all; the native-panel library; the share of strong non-1:1 genes outside the shared space (not an export gene).
- **Probe count** (human): zonation shares of human non-1:1 genes by probe count (1, 2, ≥ 3).

### H1b · The human counterparts of that zonation are flat

- **Prediction** (62: human members flat or absent): the human side of a strongly zonated probe-clean mouse non-1:1
  gene is zonated less often than the human partner of a strongly zonated probe-clean mouse 1:1 gene.
- **Counterpart.** Non-1:1: the human members of the gene's orthogroup; 1:1: its HCOP partner. Status, in this order:
  zonated (any PT-expressed probe-clean member zonated), PT-expressed but flat, not PT-expressed, no member on the
  panel.
- **Test.** Fisher, zonated against the rest, one-sided (lower for non-1:1). Several human members can only raise the
  non-1:1 side's chance, so the test is conservative.
- **Sensitivity** (probe sensitivity, the literature's Q17): only human members with ≥ 2 probes count as measured.

### H2 · Where the non-1:1 zonation sits

- **Test.** Among strongly zonated probe-clean mouse genes: χ² of peak segment (S1, S2, S3) by class (non-1:1 against
  1:1), two-sided; and the Mann–Whitney test of peak position, two-sided.
- Descriptive beside it: the peak-position density of each class, and the peak segment of the label view (W5).

### H3 · Sex: the non-1:1 zonation is partly a male program, and sex moves positions

- **H3a.** Among strongly zonated probe-clean mouse genes, the share that is male-biased in Xiong 2023 (BH-adjusted
  p < 0.05 and female − male log2FC < −1 in any of PT-S1, S2, S3) is larger for non-1:1 than for 1:1 genes. Fisher,
  one-sided. 62 found 27% against 16% post hoc; here it is pre-registered.
- **H3b · Position moves with sex.** Among strongly zonated probe-clean mouse genes with TPM ≥ 1 in some PT segment of
  both sexes in Chen 2023, the share whose argmax segment differs between male and female microdissection is larger for
  non-1:1 than for 1:1 genes. Fisher, one-sided.
- **H3c** (descriptive): the peak positions of male-biased and female-biased zonated genes against unbiased zonated
  genes (Mann–Whitney, two-sided), and section L's segment claims about sex.
- **H3d · Never tested in kidney** (literature Q15). Over all PT-expressed probe-clean mouse genes measured by Xiong:
  sex-biased (male or female) non-1:1 against 1:1, Mantel–Haenszel over deciles of log2 CP10k, **two-sided** (the one
  gene-age study predicts no enrichment). Male and female bias are also reported separately.

### H4 · Paralogs co-position rather than partition the PT

- **Units.** Pairs of zonated probe-clean paralogs in the same orthogroup and species.
- **Statistics.** Peak separation (mouse units; human curve position), and 1 − r of their mean centred curves (rows
  standardised over the grid points both specimens cover).
- **Null.** 1,000 random pairs per observed pair, matched member by member on quintiles of log2 mean CP10k and of Z
  among zonated probe-clean genes, drawn outside the family (`paralog_zonation.matched_pair_draws`; seed 86).
- **Tests.** The mean over observed pairs against the 1,000 draw means, two-sided (`global_two_sided`); per family, the
  mean percentile against 0.5 (Wilcoxon, two-sided).
- **Reading.** Co-positioning when both statistics are below the null at p < 0.05; partitioning when both are above;
  otherwise neither.
- **Stop rule** (62's): the test runs in a species with ≥ 8 families holding ≥ 2 zonated probe-clean paralogs.
- **H4b · Do label means already separate them?** A pair is separated when, in both specimens, the peaks are ≥ 0.5
  units apart in the same order (human: 0.15 curve position). For each separated pair: do the label peaks (W5) differ
  in both specimens?
- **H4b′ · Continuous against discrete** (literature Q6, D16): the number of pairs separated on the curve (H4b) against
  the number whose label peaks differ in both specimens; the pairs found only on the curve are what the curve adds.
  Also the share of zonated paralog pairs with the same peak segment against the chance share from the marginal peak
  segments of zonated genes (D10; descriptive).
- **H4c · How the handover looks** (mouse, descriptive). For separated pairs, the log2 ratio of the two genes on 40
  equal-count U_PT bins per mouse (85's P1–P2), fitted with 85's edges (E1–E5) and agreed between the mice (R1–R2). The
  main edge's class (sharp, graded, gradual) and location. 85's Q3 label-step width, recomputed on U_PT bins, is the
  reference at the boundaries.

### H5 · More diverged paralogs are more positionally diverged

- **Prediction** (duplicate-gene literature, section L): expression divergence grows with sequence divergence.
- **Pairs.** Every Ensembl 98 within-species paralog pair (within_species_paralog or other_paralog) of PT-expressed
  probe-clean genes in which at least one member is zonated.
- **Statistics.** Protein identity, the mean of the two directional identities; similarity r, the Pearson correlation
  of the two genes' mean centred curves.
- **Test.** Spearman ρ(identity, r), one-sided (ρ > 0). 95% interval by resampling orthogroups (Compara connected
  components of the pairs; 2,000 resamples, seed 86). Mouse primary, human reported.
- **Beside it:** within_species_paralog pairs only; pairs inside one HCOP orthogroup only; dS < 2 in place of identity;
  the floor: the median r of random pairs matched member by member on quintiles of log2 CP10k and Z among PT-expressed
  probe-clean genes, outside the family (20 draws per pair).

### H6 · Genomic neighbours share positions

- **Pairs.** H5's pairs.
- **Tests.** Tandem pairs (same chromosome, gene midpoints < 1 Mb apart) have larger r than the others (Mann–Whitney,
  one-sided). Among tandem pairs, r falls with distance (Spearman, one-sided).
- **Beside it:** a rank regression of r on identity and tandem status, to show whether either survives the other.
- **H6b · Paralogs against non-paralog neighbours** (D06 against D07). Neighbour pairs: PT-expressed probe-clean genes
  on the same chromosome < 1 Mb apart that are not Ensembl paralogs and not in one orthogroup, with at least one member
  zonated. Mann–Whitney of r, tandem paralog pairs against neighbour pairs, **two-sided**; the difference in medians is
  the reported effect size. Mouse primary.
- Positions: Ensembl 98 cDNA headers (mouse), r110 (human).

### H7 · Family units across species, with the 1:1 analysis's own test

This is the user's earlier question: merge each paralog family's counts and run the same gene analysis as the 1:1
orthologs.
- **Family unit.** An orthogroup (support ≥ 3) that is not 1:1, with ≥ 1 family-clean panel member in each species.
  Its counts are the sum of those members' counts in each species' units (`family_sums`). Eligible when the sum's
  equal-specimen detection is ≥ 2% on U_all (71's universe rule).
- **Model.** 71 unchanged: `standard_gene_model` on U_all with 71's library, position and grid; the species split and
  its two balanced relabelings (dispersion from the species fit); `specimen_shape_lr`; `moderated_f` with covariate log
  mean count. The variance prior is estimated over 71's universe and the family units together; guard: the same
  moderation on 71's genes alone reproduces 71's F (rtol 1e-6). q over the family units alone is reported.
- **Called** (for the comparison): BH q ≤ 0.05 over the union of 71's universe and the family units, one multiplicity
  frame for both sides.
- **Test.** The share of family units called against 71's genes that are 1:1 in both inventories:
  Mantel–Haenszel odds ratio over deciles of log mean count, **two-sided**. Duplicates diverge more in expression (a
  higher share), or dosage keeps the family sum conserved (no difference or lower).
- **Calibration.** Family calls are claimable only if the pooled share of relabeled family p ≤ 0.05 is ≤ 0.075 (71's
  gene-level rule).
- **H7c · Divergence grows with the expansion** (D12, literature Q4; descriptive): across family units, Spearman of
  the moderated-F z against log2(mouse panel members / human panel members), and against the total member count.
- **H7b · Which 1:1 species calls depend on the paralog kept.** For each 71 gene whose mouse or human member lies in a
  multi-member orthogroup: its 71 call (q ≤ 0.05) against its family unit's call, and the correlation of the two
  centred species-difference curves (Δ). Reported as a table; no test.

### H8 · The single-copy gene against its co-orthologs

- **Families.** One species holds exactly one PT-expressed probe-clean member (the single-copy side S), and the other
  holds ≥ 2 PT-expressed family-clean members (the expanded side E). The S gene must be zonated (W4).
- **Statistics** on the shared grid, from W's centred curves (each species' own fit; the family sum of E is fitted with
  W's model as one more column): r(S, each E member), r(S, E sum).
- **H8a · Subfunctionalisation.** The E sum matches S better than the best single E member: r(S, sum) > max r(S, Eᵢ).
  Sign test against one half, two-sided, over families where ≥ 2 E members are zonated.
- **H8b · Ortholog conjecture.** r(S, best Eᵢ) is lower than r for 1:1 pairs. Each family's r(S, best Eᵢ) is placed as
  a percentile among the 1:1 pairs of the same S species whose S-side gene is zonated and in the same Z quintile
  (quintiles over that species' zonated genes); the percentiles against one half, Wilcoxon one-sided (lower). 1:1
  pairs: export pairs whose genes are both PT-expressed, probe-clean and 1:1 in their inventory. Curves on the grid
  points all four specimens cover.
- **B · The 1:1 baseline** (D15, descriptive): the distribution of cross-species r for 1:1 pairs with either gene
  zonated, and the share with the same peak half (curve position below or above 0.5).
- **H8d · Dosage sharing in level** (D03, descriptive): for families with one member on one side, each gene's
  within-species percentile of log2 CP10k among PT-expressed genes. Sign test over families: the single-copy gene's
  percentile above the expanded side's top member, and the expanded side's sum closer to it than the top member is.
- **H8c** (descriptive): is the positionally best Eᵢ the sequence-closest one (Ensembl human–mouse identity) and the
  reciprocal-best kept partner? Against chance 1/k.

### H9 · Genes without an ortholog (descriptive)

PT-expressed probe-clean no_ortholog and no_consensus genes in each species: zonation, Z, peak, abundance. Kap is
expected (section L).

### H10 · Probe-based paralog measurement agrees with microdissection as well as 1:1 genes do

- **Test.** Strongly zonated probe-clean mouse genes with TPM ≥ 1 in some PT segment of Chen 2021 (male): agreement of
  our peak segment (W4) with Chen's argmax segment (mean of replicates), non-1:1 against 1:1. Fisher, two-sided.
- **Beside it:** the label view's peak (W5) in place of the curve peak; probe-ambiguous zonated genes, which carry no
  claim.

### H11 · AKI: genes outside the 1:1 pool after injury

76 tested only the 1:1 pool; 468 detected genes outside it were left out.
- **Inputs.** 76's construction unchanged (`aki_inputs` on U76, its library, position and partitions). The non-1:1
  genes are the AKI export's genes outside the pool, minus mt genes, with equal-specimen detection ≥ 2% over the four
  AKI mice on 76's units. Their counts come from the export.
- **Model.** 76's: `standard_gene_model`, `specimen_shape_lr`, `moderated_f` with the prior estimated over U76 and the
  new genes together. q is BH over the new genes. Gene claims need the new genes' gene-level calibration share ≤ 0.075.
- **H11a.** The share of positional calls of non-1:1 genes (1:many, many:1, many:many) against U76 genes that are 1:1
  in the mouse inventory: Mantel–Haenszel over deciles of log mean count, two-sided. Called: BH q ≤ 0.05 over the union
  of U76 and the new genes; guard: the moderation on U76 alone reproduces 76's F (rtol 1e-6).
- **H11b · Loss.** Among genes PT-enriched in the control mice (76's rule: log2 PT/non-PT rate ≥ 1 in both controls, on
  the QC units with 76's library), the share with replicated loss is larger for non-1:1 than for 1:1 genes (Fisher,
  one-sided). Replicated loss: P̄ ≤ −0.25 in the condition and P̄ < 0 in all four IR × control pairs. P̄ is 76's level
  (the grid mean of the centred Δ, log2), centred by 76's m(x) from U76. Profile specificity (|P̄ relabeled| < 0.25 ×
  |P̄ condition| in both relabelings) is reported beside it.
- **H11c · Where.** Among genes with replicated loss, the share whose maximal loss lies in S3 (76's segment cuts):
  non-1:1 against 1:1, Fisher two-sided.
- **H11d** (descriptive): section L's named AKI genes (`literature_aki.csv`), in or outside the pool: the condition
  P̄, pair agreement and q, with a verdict "consistent" when P̄ has the expected sign with |P̄| ≥ 0.25 and all four
  pairs agree.
- **Bounds:** the AKI screen fails calibration (issue 9) and the package drops much of the injured S3 (issue 7).

### H12 · Named claims from the literature

`non-one-to-one/literature_claims.csv`, each claim read as it is stated:
- **peak** (mouse): the gene is zonated and each mouse's own peak segment is among the expected ones; (human) the label
  peak is among the expected ones;
- **order** (mouse): the genes' peaks increase in the listed order, and each pair is separated by H4b's rule;
- **present**: PT-expressed (detection ≥ 2% in both specimens) or not, as expected;
- **flat** (mouse): PT-expressed and not zonated;
- **shoulder** (Kap): CP10k ≥ 1 in the S1 label pseudobulk of both mice, with the peak in S3.

A verdict is "consistent", "not consistent", "not zonated" (a peak claim on a flat gene) or "not measurable" (not on
the panel, not PT-expressed where a profile is needed, or probe-ambiguous).
- **Thin-limb check** (literature Q16): Aox3, Sult1a1, Gstm2, Fst and Gdf15 (thin-limb genes in microdissection) on
  U_PT and on U_all: detection and peak. On U_all they should peak at the end of mouse S3; on U_PT they should not be
  zonated toward it.

### H13 · Dominance inside expanded families (descriptive)

For families with ≥ 2 PT-expressed members in a species: the share of the family's U_PT counts carried by the top
member, mouse against human.

## Multiplicity

The primary p-values are H1 (mouse), H1b, H2 (χ²), H3a, H3b, H3d, H4 (two statistics per species), H5 (mouse), H6 (two), H6b,
H7, H8a, H8b, H10, H11a, H11b and H11c. Each is reported raw and after BH over this list.

## Bounds

- Two male control mice and one male human donor (two sections). Every species statement is descriptive; the tests
  are across genes, not across animals.
- Human positions are the weakest measurement (issue 21) and the human transitions cross (issue 6).
- C14: mouse positions inside S1 are noisy.
- Probe panels differ by species; level comparisons across species are not made.
- All specimens are male (memory: all specimens male): sex statements rest on external data.
- AKI: issues 7 and 9.

## What would make us say the non-1:1 genes add nothing on this curve

- H1 fails (OR < 1.5 or p ≥ 0.05): the zonation that the 1:1 space misses is not enriched.
- H7 shows no difference from 1:1 genes and H7b finds no 71 call that changes with the family unit.
- H4, H5 and H6 all give "neither".
- H11a and H11b give no difference.

## Figures (plots only)

- 86a: strong-zonation share by class and species; peak-position densities by class (mouse units).
- 86b: family heatmaps on the curve: mouse members, human members and family sums (centred curves).
- 86c: paralog co-positioning against the null; identity against r; tandem against dispersed pairs.
- 86d: family units across species: call shares against matched 1:1 genes and the relabelings; example Δ curves.
- 86e: sex: male-biased shares by class; example genes whose microdissection peak moves with sex.
- 86f: AKI: P̄ of PT-enriched non-1:1 and 1:1 genes; example curves.
- 86g: section L's named genes on the curve, both mice.
- 86h: H0's positive-control pairs and the thin-limb check genes on the curve (U_PT and U_all).

## Notebook 86 results

`analysis/notebooks/86_pt_non_one_to_one_orthologs.ipynb` (logic `86.non_one_to_one.1`) writes to
`results/pt_non_one_to_one_package/`. It followed the protocol frozen in 6f54b8b. Every guard passed (G1–G6, plus the
H7 and H11 moderation guards: 71's and 76's F are reproduced from their caches). Units: U_PT holds Ctrl1A2 1,737,
Ctrl1A4 1,844, HUK1_COR1 2,625 and HUK1_MED1 2,636 PT-identity units. PT-expressed genes: mouse 10,773 (zonated 1,291,
strong 440), human 11,574 (zonated 555, strong 153).

### Primary tests

| Test | Result | p | BH |
|---|---|---|---|
| H0 · 1:1 paralog control separates (mouse) | Slc5a2 → Slc5a1: Δ 2.05, 2.11 units; Slc7a9 → Slc3a1: 2.28, 2.24 | passes | |
| **H1 · non-1:1 enriched in strong mouse zonation** | 8.2% against 3.8%; MH OR **2.24** | < 10⁻⁵ | < 10⁻⁴ |
| H1b · human counterparts flat | zonated 8/41 against 82/380; OR 0.88 | 0.47 | 0.60 |
| H2 · peak segment by class | non-1:1 S1/S2/S3 13/8/20 against 179/35/166 | 0.050 | 0.13 |
| H3a · male bias among strong genes | 10/41 against 63/380; OR 1.62 | 0.15 | 0.25 |
| H3b · sex moves the microdissection peak | 14/38 against 104/359; OR 1.43 | 0.20 | 0.31 |
| **H3d · sex bias across PT-expressed genes** | 16.3% against 10.7%; MH OR **1.50** (male 1.90, female 1.16) | 0.001 | 0.004 |
| H4 · paralog peak separation (mouse) | 0.73 against 1.08 units | 0.008 | 0.024 |
| H4 · paralog shape dissimilarity (mouse) | 0.75 against 0.99 | 0.056 | 0.13 |
| H5 · identity against similarity (mouse) | ρ 0.016 | 0.10 | 0.19 |
| **H6 · tandem paralogs more similar** | median r 0.31 against 0.07 | 10⁻⁴ | 4 × 10⁻⁴ |
| H6 · similarity falls with distance (tandem) | ρ −0.20 | 0.001 | 0.004 |
| **H6b · tandem paralogs against non-paralog neighbours** | median r 0.31 against 0.04 | < 10⁻⁴ | 4 × 10⁻⁴ |
| H7 · family units called against 1:1 genes | 35.4% (139/393) against 32.6%; OR 0.96 | 0.74 | 0.83 |
| H8a · family sum beats best member | 0 of 1 | 1 | 1 |
| H10 · agreement with microdissection | 26/40 against 254/363 | 0.59 | 0.71 |
| H11a · AKI positional calls | 65/144 against 4,327/10,286; OR 1.42 | 0.081 | 0.16 |
| H11b · AKI replicated loss | 33/54 against 596/896; OR 0.79 | 0.83 | 0.88 |
| H11c · AKI loss in S3 | 24/33 against 476/596 | 0.37 | 0.52 |

### What holds

1. **H0 · The curve separates paralogs when they really split.** SGLT2 (Slc5a2) peaks in S1 and SGLT1 (Slc5a1)
   rises toward S3, and b(0,+)AT (Slc7a9) falls while rBAT (Slc3a1) rises, in both mice by more than two segment
   units. Human SLC5A2 and SLC5A1 peak in S1 and S3 labels. The co-positioning results below are therefore readable.
2. **H1 · The mouse PT zonation that the 1:1 space misses is real and replicates 62 on the new curve.**
   - Strong-zonation shares: 1:1 3.8%, 1:many 9.9%, many:many 14.2%, but many:1 only 1.7%. The excess sits in mouse
     genes of **mouse-expanded** families, not in mouse single-copy genes whose human side expanded.
   - MH OR 1.84–4.05 under every sensitivity (support ≥ 2 or ≥ 5, strict probes, U_all, native library).
   - 41 strong probe-clean non-1:1 mouse genes, 25 of them outside the shared space. The largest Z: Cyp2d12, Nat8f5,
     Slc22a19, Scd1, Cyp2d26, Ugt2b38, Serpina1a, Clec2h, Akr1c14, Ceacam2, Akr1c18, Ugt2b37, Ces1e, Ugt2b34, Cyp4a14,
     Cyp4a12a, Slco1a6, Cyp4a10, Cyp2d9, Cyp2j13 and Slco1a1 (`w_calls_mouse.csv`).
   - Human: OR 1.57, p 0.17 (not enriched at support ≥ 3; 1.85 and 2.69 at support ≥ 2 and ≥ 5).
3. **H3d · Sex-biased mouse PT genes are enriched in non-1:1 families** (the literature's open question; never tested
   in kidney). Over all PT-expressed probe-clean genes, against an expression-matched 1:1 background: any sex bias OR
   1.50 (p 0.001), **male bias 1.90 (p 10⁻⁴)**, female bias 1.16 (p 0.40). The one gene-age study predicted no
   enrichment (PMID 37917687). Sex bias comes from Xiong 2023's female-vs-male PT data, so our all-male cohort does not
   confound it.
   - H3c (descriptive): in our male mice, male-biased zonated genes peak late (median 2.17 units, early S3) and
     female-biased ones early (0.57, S1), against 1.48 for unbiased genes (p 10⁻⁴ and 2 × 10⁻⁷). This fits the literature's S2/S3
     location of the mouse sex program.
4. **H6/H6b · Tandem paralogs share positions, and it is specific to paralogs.**
   - Tandem pairs (< 1 Mb): median r 0.31 against 0.07 for dispersed paralogs and 0.04 for non-paralog neighbours
     within 1 Mb (19,651 pairs). Human: 0.48 against 0.07 and 0.04.
   - Within tandem pairs, similarity falls with distance (ρ −0.20 mouse, −0.21 human).
   - In mouse, tandem status predicts shared position and sequence identity does not (rank regression t 3.7 against
     0.5). Of the literature's competing hypotheses, co-regulation of tandem duplicates (D06) wins over a plain
     neighbourhood effect (D07).
5. **H4 · Paralogs sit closer together than matched random genes** (0.73 against 1.08 units, p 0.008), but their
   shapes are not significantly more alike (p 0.056), so the frozen reading is "neither". 59% of zonated paralog pairs
   share a peak segment against 40% by chance (D10, descriptive). Human: not run (2 families with ≥ 2 zonated
   paralogs).

### What does not hold

- **H1b:** the human counterparts of strongly zonated mouse non-1:1 genes are zonated about as often as the human
  partners of strong 1:1 genes (20% against 22%). Restricting to human members with ≥ 2 probes changes nothing.
- **H5:** across 5,959 Ensembl paralog pairs, protein identity does not predict positional similarity in mouse
  (ρ 0.016). It does among young pairs (dS < 2: ρ 0.17, interval 0.003–0.32) and in human (ρ 0.06, p 7 × 10⁻⁴;
  dS < 2 0.52). Descriptive.
- **H7 · Family sums through the 1:1 test are called as often as 1:1 genes** (35% against 33%). The family test also
  fails calibration (relabeled share 0.089 > 0.075), so individual family calls are not claimable. H7c: no relation to
  the expansion ratio (ρ 0.06).
- **H7b · The 1:1 species calls barely depend on the paralog kept.** Of 122 genes 71 calls inside multi-member
  orthogroups, 116 are called at the family level too; 6 are lost (Akr1c18, Gnai1, Zrsr2, Nr2f1, Ier3ip1, Slc52a2) and
  32 families are called only as sums.
- **H8:** only 4 families have a zonated single-copy gene and ≥ 2 PT-expressed co-orthologs; H8b was not run (n < 5).
  Baseline B: 1:1 pairs with either gene zonated have median cross-species r 0.39, and 60% share a peak half.
- **H10:** the curve's peak segment agrees with male microdissection for 65% of strong non-1:1 genes and 70% of 1:1
  genes (p 0.59): paralogs are measured as reliably as 1:1 genes.
- **H11 · AKI:** the 468 genes outside 76's pool fail gene-level calibration (0.114), so no gene claim. Non-1:1 genes are
  not lost more often than 1:1 genes after injury (61% against 67% of PT-enriched genes), nor more often in S3.
- **H13:** the top member carries a similar share of family PT counts in both species (median 0.72 mouse, 0.69 human).

### Positions (H4b, H4c) and named claims (H12)

- **Continuous against discrete (H4b′).** 16 of 39 zonated mouse paralog pairs separate on the curve; the S1/S2/S3 label
  peaks separate 15 of the same 16. The curve adds one pair: Ugt2b38 (1.29 units) before Ugt2b37 (1.81) inside S2,
  which matches male microdissection (Ugt2b38 S1+S2, Ugt2b37 S2).
- **Handovers (H4c, descriptive).** Of the 16 separated pairs, most hand over gradually or graded: Ugt2b34 to
  Ugt2b37/38 at the S1/S2 boundary (widths 0.65–0.93), Cyp2j13 to Cyp2j11 and Ceacam1 to Ceacam2 at S2/S3, Cyp4a12a to
  Cyp4a10/14/31 inside S2 (0.60–0.71). Two are sharp: Acaa1b/Acaa1a inside S1 (0.44) and Ces1e/Ces1f in S3 (0.21).
  The U_PT label steps are 0.18–0.31 wide.
- **H12 · 48 of 60 literature claims are consistent, 4 are not, 8 are not measurable.**
  - Consistent: Kap peaks in S3 with an S1 shoulder (the intact-male androgen pattern), Slco1a1 and Cyp2j13 in S2,
    Akr1c18 rising into S3 in males, the Cyp2d, Ugt2b and Ces1 paralog orders, Nat8f1/Nat8f5/Nat8 (no segment data in
    any publication before), Slc22a19 late, Gstm1 S3, Mt1 and Gsta4 S1, Slc7a12 absent in male PT, and on the human
    side SLC5A2/SLC5A1, SLC22A24 in S3, CES2 in S1/S2, SLC22A6/A8 in S1/S2, SLCO1A2 and SLC22A9/A10 not PT-expressed.
  - **Not consistent: the Cyp4a family.** Microdissection puts Cyp4a10 and Cyp4a14 in S2 and Cyp4a12a in S3; our curve
    and our labels put Cyp4a12a in early S2 (1.19 units) and Cyp4a10/14/31 at the S2/S3 boundary (2.08–2.20). 62 found the
    same order on the old coordinate, so the conflict is not the coordinate's. The male microdissection levels of
    Cyp4a12a are low (9/26/55 TPM) and its probe status is "near" (one probe within 4–5 mismatches of a paralog).
  - Not measurable: probe-ambiguous genes (Ces2c, Slc22a30, Gsta2, CYP4A11, CYP4F2), genes below detection (Slc22a29,
    Slc22a28) and Gsta1 (no probe in our panel).
- **Thin-limb check.** On U_all, Gstm2 and Gdf15 are zonated toward the end of mouse S3 (peak 2.87); on U_PT neither is
  zonated. U_PT removes the thin-limb end, as intended.
- **H11d · named AKI genes** (descriptive): consistent with the literature are Kap loss (P̄ −2.0, all four pairs),
  Slc22a30, Acsm2 and Slc7a13 loss, and Mt1, Mt2, Cxcl1, Havcr1 and Ifitm3 gain. Slc7a12 is not gained (−0.73),
  Slc5a12 and Ifi27 are not replicated, and Saa1/3, Ccl2/7, Serpina3n and Sprr2f are below detection.

### Recorded deviations (after the first run)

1. **H11's usable genes.** The first run's P̄ included genes whose specimen fits hit zero counts, giving P̄ values up to
   ±10³ (for example Cxcl2, Vmn1r19, Nat8f3). 76 uses only converged, non-separated genes. The rerun applies 76's rule,
   extended to every specimen fit: U76 10,597 of 10,701, new genes 457 of 468. H11b moved from 37/59 against 606/916
   (OR 0.86, p 0.76) to 33/54 against 596/896 (OR 0.79, p 0.83); every verdict is unchanged.
2. **H8b** needed ≥ 5 families for the Wilcoxon test; 4 qualified, so it was not run. The protocol named no minimum.
3. **Figures:** Wilson intervals replace normal-approximation intervals in 86a; 86c gained a fourth panel (H6b's three
   groups); 86f clips P̄ at ±8 for display.

### Bounds

Two male control mice, one male human donor, two male IR mice. The tests run across genes, not animals. Human
positions are the weakest measurement, and the human transitions cross. Probe panels differ by species, so no level
comparison across species is made. The sex statements rest on external female-vs-male data.

### Tables and figures

- Tables (`tables/`): `w_calls_{mouse,human}.csv`, `h0_positive_control.csv`, `h1_enrichment.csv`, `h1_by_class.csv`,
  `h1_human_by_probe_count.csv`, `h1b_*`, `h3_tables.csv`, `h3b_moved_genes.csv`, `h3c_*`, `h3d_*`, `h4_*`,
  `h4c_handovers.csv`, `h5_*`, `h6_*`, `h6b_*`, `h7_*`, `h7b_*`, `h8_*`, `h8d_*`, `h9_no_ortholog.csv`,
  `h10_microdissection.csv`, `h11_*`, `h11d_named_aki_genes.csv`, `h12_literature_claims.csv`,
  `h12_thin_limb_check.csv`, `h13_dominance.csv`, `primary_tests.csv`.
- Figures (`figures/`): 86a inventory, 86b families, 86c paralogs, 86d family units, 86e sex, 86f AKI, 86g literature
  genes, 86h controls.
