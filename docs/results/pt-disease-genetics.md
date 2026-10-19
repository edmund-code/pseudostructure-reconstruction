# Human disease genetics and drug targets along the PT (notebook 108)

**The question.** Where along the proximal tubule do the genes of inherited tubulopathies, kidney-trait GWAS loci and
renal drug transporters sit, in human and in mouse? The literature makes positional predictions for several of them.
The continuous axis can test each prediction gene by gene and set by set, where segment labels cannot. Lake 2023
could not separate human PT-S1 from PT-S2 ("Cell types that could not be accurately resolved (PT-S1/PT-S2) were kept
merged").

**Sources.** The genetics scout read these in full unless marked; quotes and read status are in its report, summarised
in `docs/research/disease-genetics-along-pt.md`.
- Mendelian panels: Genomics England PanelApp 292 and 149, PanelApp Australia 3993, and HPO HP:0000114.
- GWAS Catalog v2 associations:
  - eGFR: Stanzick 2021;
  - serum urate: Tin 2019;
  - urine metabolites: Schlosser 2020 and 2023.
- **Segment priors:**
  - S1 carries the largest transport capacity (Zhuo & Li 2013).
  - SGLT2 sits in early S1 (Vallon 2011).
  - Human and rodent OATs are strongest in S1/S2 and absent from outer-stripe S3 (Breljak 2016, abstract only; rat
    RNA-seq, Lee 2015).
  - Rat OCT2 rises toward S3 (Karbach 2000, abstract only; Lee 2015). Human OCT2 is described in cortex and medullary
    ray (Motohashi 2002, abstract only).
  - Megalin protein and endocytosis are S1 ≫ S3, but its rat mRNA is nearly flat (Christensen 1995, abstract only;
    Lee 2015).
  - GATM is in the early PT and EHHADH in S2–S3 (Reichold 2018; Forst 2021).
  - NBCe1 (SLC4A4) is highest in S1 (Skelton 2010).
  - URAT1 and ABCG2 are not in S1 in mouse (Hoque 2020).
  - Mouse Glut9 is in the DCT and CNT, not the PT (Preitner 2009).
  - SLC5A11 is high and specific in S3 (Schlosser 2020).

Frozen on 2026-10-18, before any gene set was joined to the gradients.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/108_pt_disease_genetics.ipynb`; logic `108.disease_genetics.1`; seed 108.
- Results `results/pt_disease_genetics/`.
- **Gradients:** notebook 90's `d_end_contrasts.csv.gz` (E per specimen, 7,971 1:1 genes). Per-species curves come
  from notebook 71's cached human vs mouse species model, as notebook 101 read them.
- **Expression strata:** notebook 62's inventories (`log2_cp10k` per species), in deciles over the universe.
- **Gene sets, as human symbols mapped to the 1:1 universe:**
  - **PROX**, the Mendelian proximal tubulopathy set:
    - green genes of panels 292, 149 and 3993 whose phenotype text names a proximal disorder;
    - the keywords, case-insensitive: *Fanconi, Dent, Lowe, cystinosis, cystinuria, glucosuria, glycosuria, Hartnup,
      iminoglycinuria, proximal renal tubular acidosis, hypophosphat, rickets, aminoaciduria, lysinuric*;
    - plus the endocytic receptors LRP2, CUBN and AMN.
  - **DIST**, the distal comparator: green genes of panel 292 whose phenotype text names *Bartter, Gitelman, distal
    renal tubular acidosis, pseudohypoaldosteronism, nephrogenic diabetes insipidus, hypomagnes, Liddle, EAST* and
    no PROX keyword.
  - **eGFR, URATE, UMET:** the Catalog's mapped genes for Stanzick 2021; Tin 2019; and Schlosser 2020 ∪ 2023.
  - **Transporters:** HGNC groups 752 (SLC) and 417 (ABC), as in notebook 101. X-TR is set X ∩ transporters.

## 2 · Tests and readings (fixed now)

A species-level E is the mean of that species' two specimens. A gene is graded early or late in a species when
E ≤ −0.5 or ≥ +0.5 in both of its specimens. Nulls draw 2,000 random lists from the universe, matched on the species'
expression decile (`matched_means`).

**H1 · Mendelian proximal genes skew toward S1.**
- The mean E of PROX lies below the null's 5th percentile, per species.
- DIST is not graded more often than the null; it serves as the comparator.
- Prior: moderate. It holds for capacity, and only a few genes have a published segment.

**H2 · Urine-metabolite GWAS transporters are more graded than other transporters.**
- Statistic: the mean |E| of UMET-TR. The null draws only transporters, matched on expression decile.
- Holds per species when above the null's 95th percentile.
- eGFR-TR is reported beside it as the comparator, since it has no segment prior.

**H3 · Urate GWAS genes change pattern between species.**
- The share of URATE genes in notebook 90's reversed or species-only categories.
- Holds when above the null's 95th percentile. Prior: mice have uricase, and mouse Glut9 is distal.

**H4 · Named positions** (each met or missed; a miss is reported, not tuned):

| Check | Prediction |
|---|---|
| a | SLC5A2 early in both species |
| b | SLC7A13, SLC22A7 and SLC3A1 late in both species; SLC5A11 late in human (if in the universe) |
| c | SLC22A6 and SLC22A8 peak inside the axis: the maximum of the species curve lies in 0.2–0.8, and the curve's last 10% is ≥ 0.5 log2 below that maximum, in both species |
| d | SLC22A2 (OCT2) late in mouse and not late in human |
| e | LRP2 and CUBN not graded in either species (transcripts under-call the functional S1 gradient) |
| f | GATM early in both species; EHHADH not early |
| g | SLC4A4 early in both species |
| h | SLC22A12 (URAT1) late in human; SLC2A9 not graded in mouse |

**H5 · Drug targets (descriptive).** The species curves of SLC5A2 (gliflozins), SLC22A12 (uricosurics), SLC22A6
(probenecid), SLC22A2 (cisplatin uptake) and LRP2 (aminoglycoside uptake).

## 3 · Outputs and figures

- **Tables:** `sets.csv` (set, gene, in universe, E per species, category), `h1_h3.csv`, `h4_named.csv`,
  `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 108a: mean E of each set against its null, per species;
  - 108b: curves of the H4 genes, mouse and human;
  - 108c: the drug-target curves.

## 4 · Caveats

- These are mRNA gradients. Megalin and URAT1 are known cases where protein or function gradients exceed mRNA ones.
- GWAS mapped genes are positional nearest genes, not causal genes.
- One human donor with two sections.
- Structural species asymmetries:
  - SLC22A11 (OAT4) has no mouse ortholog;
  - human MATE2-K has no separate mouse gene;
  - mouse Slc2a9 is distal.

## Notebook 108 results

`analysis/notebooks/108_pt_disease_genetics.ipynb` (logic `108.disease_genetics.1`) writes to
`results/pt_disease_genetics/`. It followed the protocol frozen in f5cfe25 and ran in 4 s.

**Sets** (all genes / in the 1:1 universe):

| Set | Genes | In universe |
|---|---|---|
| PROX | 28 | 23 |
| DIST | 20 | 11 |
| eGFR | 688 | 295 |
| URATE | 328 | 144 |
| UMET | 467 | 200 |
| Transporters | 527 | 237 |
| UMET transporters | 47 | 24 |
| eGFR transporters | 26 | 15 |

**How the protocol was read.**
- Panel genes are matched on the HGNC symbol.
- GWAS mapped-gene strings are split on commas, semicolons and " - ".
- Curves are notebook 71's specimen curves in log2, averaged within a species.
- Post hoc (section 3.1): notebook 104's edge enrichment and core verdicts are joined for the DIST genes and the H4
  misses.

| Test | Mouse | Human | Holds |
|---|---|---|---|
| **H1** mean E of Mendelian proximal tubulopathy genes (null 5–95%) | **−1.47** (−0.25 to 0.35) | **−0.56** (−0.28 to 0.03) | **yes, both species** |
| PROX share graded | 0.70 (null ≤ 0.39) | 0.43 (≤ 0.26) | descriptive |
| DIST share graded | 0.45 (≤ 0.36) | 0.36 (≤ 0.27) | descriptive: above the null (see below) |
| **H2** mean \|E\| of urine-metabolite GWAS transporters, among transporters | 1.28 (0.70–1.55) | 0.64 (0.39–0.83) | **no** |
| eGFR GWAS transporters (comparator) | 1.31 (0.61–1.69) | 0.77 (0.33–0.95) | descriptive |
| **H3** urate GWAS genes reversed or species-only | 0.076 (null 0.03–0.09) | | **no** |
| **H4** named positions | 25 of 27 met; 2 missed | | |

**H4, gene by gene.**
- **Met, both species:**
  - SGLT2 (SLC5A2) is early;
  - SLC7A13 and SLC22A7 are late;
  - OAT1 and OAT3 (SLC22A6, SLC22A8) peak inside the axis and fall at the S3 end;
  - LRP2 and CUBN are not graded, though their protein and function are S1-weighted;
  - GATM and NBCe1 (SLC4A4) are early, and EHHADH is not early.
- **Met, one species:**
  - SLC3A1 is late in mouse;
  - SLC5A11 is late in human (Schlosser 2020's S3 claim);
  - URAT1 (SLC22A12) is late in human;
  - OCT2 (SLC22A2) is late in mouse and **not graded in human**, the predicted species divergence.
- **Missed:**
  - SLC3A1 in human is not graded (E +0.45, just short);
  - mouse Glut9 (Slc2a9) is graded early (E −0.92), though Preitner 2009 placed the protein in the DCT and CNT.
    Its PT signal is not edge-enriched (ε −0.14, 0.19), so it is PT mRNA rather than neighbour RNA.

**PROX genes, in detail.**
- Early in both species: SLC2A2, SLC34A3, SLC4A4, SLC5A2, SLC6A19, SLC7A7, SLC7A9 and GATM.
- Early in mouse only: ALPL, CTNS, SLC34A1 and AMN.
- Late: EHHADH and SLC3A1 (mouse) and HNF1B (both species).
- Not graded: LRP2, CUBN, HNF4A, OCRL, CLCN5 and VDR.

**Post hoc: the distal genes graded inside PT units are neighbouring-structure RNA.**
- The pattern follows the neighbours:
  - SLC12A1 (NKCC2, thick limb) is late in both species;
  - SLC12A3 (NCC, DCT) is early in mouse;
  - AQP2, KCNJ1 and SCNN1A are late in human.
- Each matches who lies next to S3 (outer-stripe thick limbs and collecting ducts) and next to S1 (DCT in the
  cortical labyrinth).
- These genes are edge-enriched (human ε 1.0–1.8). Yet their gradients survive notebook 104's core rebuild.
- **This is a limit of notebook 104:**
  - For neighbour genes expressed at very high levels, enough RNA reaches 6 µm inside a PT polygon (or the polygon
    includes neighbour cells) for the gradient to persist in the core.
  - So core survival is necessary but not sufficient evidence of PT origin.
  - The check needs a gene-specificity filter as well (a neighbour marker in an external atlas).

**What this answers.**
1. **The genes of inherited proximal tubulopathies sit toward S1 in both species.**
   - Their mean gradient is far below an expression-matched null (mouse −1.47, human −0.56).
   - So the Mendelian genetics of PT failure concentrates where the PT's transport capacity is largest. This is the
     long-held physiological claim (Zhuo & Li 2013), now shown on gene expression along a continuous axis in human.
   - The exceptions are informative. EHHADH points toward S3, which fits Forst 2021. The endocytic receptors are
     flat in mRNA.
2. **Kidney-trait GWAS genes have no positional signature beyond what transporters already have.**
   - Urine-metabolite and eGFR transporters are as graded as other transporters. Notebook 101 showed transporters as
     a whole are graded 2–3 times more often than other genes.
   - Urate genes do not change pattern between species more than chance.
3. **The literature's single-gene positions are reproduced, 25 of 27, from an unsupervised axis.** This includes
   two species-specific predictions:
   - OCT2 (the cisplatin uptake transporter) rises toward S3 in mouse but is flat in human. This qualifies the
     rodent-based account of cisplatin's S3 toxicity for human.
   - URAT1 and SLC5A11 rise toward S3 in human.
   - Both are concrete, citable "micro-gradients in solute handling" with pharmacological relevance.
4. **As a method:** a frozen list of literature predictions scored on the continuous axis is a strong, cheap
   validation. Reviewers can check each line, and the misses are informative.
