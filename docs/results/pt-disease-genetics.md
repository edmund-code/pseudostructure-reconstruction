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
