# SAMap and other many-to-many cross-species methods for genes outside 1:1 orthology (2026-10-14)

The user asked: "what about using samap to analyze non 1:1 orthologs". This note answers from the full text of SAMap
and of the four papers that test it or its successors, against what notebooks 64, 86 and 88 already did. It changes no
result and freezes no protocol.

## Read status

| Paper | What was read |
|---|---|
| Tarashansky et al. 2021, *eLife*, PMID 33944782 (SAMap) | Full text: main text, methods, figure and figure-supplement legends, decision letter and author response. Supplementary files 1–7 (tables) not read. SAMap README v3.0.0 read. |
| Song et al. 2023, *Nat Commun*, PMID 37838716 (BENGAL benchmark) | Full text: main text, methods, figure and table legends. Supplementary Information not read. |
| Liu et al. 2023, *Genome Res*, PMID 36526433 (CAME) | Full text and the Supplemental Material PDF (legends; Figs S1, S3, S4 as images). GO-term tables and code not opened. |
| Rosen et al. 2024, *Nat Methods*, PMID 38366243 (SATURN) | Full text and the Supplementary Information PDF (Notes 1–6, Figs 1–5, Tables 1–3). Reporting Summary (scanned) not read. |
| Zhong & Han 2025, *Nucleic Acids Res*, PMID 39778870 (benchmark across 20 species) | Full text; all 55 supplementary figure legends; Tables S1–S6 and the human–mouse and ablation rows of Tables S9–S10. Other supplementary panels not inspected. |

Per-paper notes: `full-text-notes/`.

## What SAMap does

- **Homology:** reciprocal BLAST of the two proteomes or transcriptomes (E < 10⁻⁶). Every homolog enters, weighted by
  bit score, so one gene can link to many.
- **Cells:** each species gets its own kNN manifold. A gene's expression in the other species is the weighted mean
  of its homologs. Cells are linked by mutual nearest neighbourhoods in the joint space.
- **Refinement:** gene-pair weights are replaced by their expression correlation over the joint manifold (negative
  correlations set to 0). This repeats three times.
- **Outputs:**
  - cell-type alignment scores;
  - gene pairs that drive each mapping (filtered by a within-species Wilcoxon over cells);
  - **paralog substitutions**: a gene whose correlation with a paralog beats its correlation with its ortholog by
    more than 0.3.
- **Calibration:** the 0.3 threshold is the call. There is no null distribution, no p-value and no replicate
  structure; cells are the observations.

## Do genes outside 1:1 orthology help at human–mouse distance?

None of the five papers finds a clear gain.

- **SAMap itself.** It was never tested on mouse vs human. Even at frog vs zebrafish, SAMap restricted to 1:1
  orthologs gave "similar, albeit weaker" mappings. The authors credit SAMap's performance mainly to its
  "atlas stitching approach", not to the many-to-many homology.
- **BENGAL.** Differences came from the integration algorithms, "not homology methods". Adding in-paralogs raised
  scores only for frog vs zebrafish, where paralogs are a large share of the genes. Human vs mouse pancreas used
  11,248 one-to-one and 664 other homologs (Table 1). For close species the authors recommend scVI/scANVI or Harmony,
  and keep SAMap for whole-body atlases or species with poor homology annotation.
- **CAME** reran its model with 1:1 homologs only:
  - zebrafish pairs lost 6.3% accuracy on average (p 7.8 × 10⁻²³);
  - every other pair showed no difference (p 0.67). These include the human–mouse kidney pairs, whose points lie on
    the diagonal (Supplemental Figs S3D, S4A).
  - Non-1:1 genes made up 60–75% of the informative genes in zebrafish pairs and 15–40% elsewhere.
- **Zhong & Han 2025.**
  - Human–mouse scores were bunched across all nine methods. SAMap is weak for close species, which the authors
    attribute to "the insensitivity of BLAST-derived features to taxonomic close species".
  - In their human–mouse ablation (one task), restricting SAMap to 1:1 orthologs cost 4.9%.
  - SATURN with all homologs scored 0.539, against 0.533 for the best 1:1-only method (Table S10).
- **SATURN.** Its human–mouse comparison was against SAMap only, and SAMap's median was equal or slightly higher
  (Supplementary Fig. 3b). No 1:1-only baseline was run on any mammalian pair.
- **Our panels** (86's inventory, PT-expressed probe-clean genes): 94% of mouse genes and 93% of human genes are 1:1.

## SAMap's analyses against what notebooks 64, 86 and 88 already did

| SAMap analysis | Our version | Result |
|---|---|---|
| Joint space not restricted to 1:1 genes | 64: the main PT paralog put in place of the kept partner (92 pairs), or the 1,051 pairs inside multi-member orthogroups dropped | Indices moved by ≤ 0.007; class counts by < 10% |
| Many-to-many homology graph | 86: HCOP orthogroups (support ≥ 3; ≥ 2 and ≥ 5 as sensitivities), Ensembl 98 paralogs, family sums of probe-clean members | Mouse PT zonation enriched for genes in mouse-expanded families (OR 2.24); family units called as often as 1:1 genes |
| Paralog substitution (paralog r − ortholog r > 0.3, no null) | 88 P1: gain ≥ 0.5 and r_alt ≥ 0.5 on the PT curve, against random genes matched on expression and zonation | 6 of 32 pairs against a null mean of 7.5 (p 0.85): no more than chance |
| Expression-refined pairing inside families | 88 P4: positional pairing against sequence pairing | Same pair in 1 of 12 against 4.4 by chance (descriptive) |
| Species-specific genes kept in the joint space | 86 H9 and 88 P6 | P6 cannot resolve partners: random human–mouse pairs reach \|r\| 0.83 on the curve |

## Why SAMap fits our design poorly

1. **It has no null for substitution.** 88's P6 found that profiles along the curve are so smooth that random gene
   pairs often correlate strongly. A fixed 0.3 correlation gap would therefore call many substitutions by chance. 88's
   matched null is what turned 6 "substituted" pairs into a negative result. SAMap correlates genes across cells, not
   along a curve, but in PT-only data most unit-to-unit variation is position, so the same problem should apply. This
   is our inference; it was not tested.
2. **It blurs continuous structure.** BENGAL found SAMap had low trajectory conservation, with developmental stages
   partly mixed. Zhong & Han found SAMap and SATURN preserved trajectories worse than scGen and Seurat v4 CCA. Our
   question is about position along the PT.
3. **Its joint graph cannot be our coordinate.** The elastic curve is the only primary coordinate.
4. **It has no replicate structure.** Units are observations; we have two mice and one donor.
5. **Probe panels.** BLAST links paralogs that share or lose probes. A run would first need 86's probe-clean or
   family-clean rule; SAMap has no mechanism for this.
6. **The refinement is circular.** Correlations are computed on a manifold built from the same weights. A reviewer
   called this a possible "tautology". The authors' 1:1-only rebuild recovered 70% of the substitutions.

## What remains

- **A citation.** 88's substitution test is SAMap's paralog-substitution idea, moved onto a tissue axis and given a
  matched null. If 88 enters the paper, cite Tarashansky et al. 2021 there. Song et al. 2023 is already cited.
- **One question PT-only analyses cannot answer:** when the human member of a mouse PT-zonated family is not
  PT-expressed (86's H1b status), is it expressed in another segment instead? SAMap's whole-body setting points to
  this question, but segment pseudobulks of the probe-clean family members would answer it without SAMap. Not
  planned; recorded here.
- **Verdict.** Running SAMap would not add a testable result for genes outside 1:1 orthology in this cohort. Its
  gene-level question is already answered with a calibrated test (88: no substitution beyond chance). The literature
  finds no gain from many-to-many homology at human–mouse distance.
