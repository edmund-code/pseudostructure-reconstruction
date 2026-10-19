# Human disease genetics and pharmacology mapped to positions along the PT (scout, 2026-10-18)

A literature and resource review for notebook 108. Every claim comes from the full text unless tagged *abstract
only*. The raw API downloads used by notebook 108 are in `data/external/panelapp_2026-10-18/` and
`data/external/gwas_catalog_v2_2026-10-18/`.

## 1 · Mendelian tubulopathy gene lists

**API notes.**
- The PanelApp API is reachable without login: `.../api/v1/panels/<id>/`, green genes = `confidence_level` 3.
- `?search=` is ignored, so filter on the client.
- PanelApp Australia uses the same API shape (`panelapp-aus.org`).

| Panel | ID | Version | Genes | Green |
|---|---|---|---|---|
| Renal tubulopathies (GE) | 292 | 6.10 | 67 | 52 |
| Nephrocalcinosis or nephrolithiasis (GE) | 149 | 6.3 | 52 | 35 |
| Renal Tubulopathies and related disorders (AUS) | 3993 | 2.10 | 134 | 117 |

**Missing panels.** Neither PanelApp has a "Renal Fanconi syndrome" or "Proximal tubulopathy" panel. HPO HP:0000114
("Proximal tubulopathy") lists 59 genes, about 40 of them mitochondrial complex-I genes.

**Proximal genes in panel 292** (green; the scout's reading of the phenotype text):
- CTNS, CLCN5, OCRL, GATM, HNF4A, SLC2A2, SLC4A4, SLC5A2, SLC22A12, SLC2A9, FAH, CYP24A1;
- VPS33B, VIPAS39, RMND1, SARS2;
- CA2 is mixed proximal and distal.

**Elsewhere:**
- Amber or red on panel 292: EHHADH, SLC34A1, SLC34A3, SLC9A3R1, XPR1, ABCG2.
- Cystinuria (SLC3A1, SLC7A9), Hartnup (SLC6A19) and iminoglycinuria (SLC36A2, SLC6A20) are on panels 149 and 3993,
  not 292.
- LRP2 is on no tubulopathy panel.

**The distal remainder:**
- thick limb: SLC12A1, KCNJ1, CLCNKB, BSND, CLDN16, CLDN19, UMOD, CASR and others;
- DCT: SLC12A3, WNK4, KLHL3, CUL3, TRPM6 and others;
- collecting duct: SCNN1A/B/G, NR3C2, AVPR2, AQP2, ATP6V0A4, ATP6V1B1, SLC4A1 and others.

## 2 · Kidney-trait GWAS

**API notes.**
- The legacy GWAS Catalog API returns 410.
- v2 works: `https://www.ebi.ac.uk/gwas/rest/api/v2/associations?pubmed_id=<PMID>&size=200&page=N`.

| Paper | PMID | Associations | Mapped genes |
|---|---|---|---|
| Stanzick 2021, eGFR | 34272381 | 927 | 688 |
| Tin 2019, urate | 31578528 | 604 | 328 |
| Schlosser 2020, urine metabolites | 31959995 | 286 | 132 |
| Schlosser 2023, plasma and urine metabolites | 37277652 | 1,299 | 450 |

**Key quotes:**
- Schlosser 2020: "gene expression data from 43,745 single murine kidney cells and from 4,524 nuclei of adult human
  kidney cells were used to identify cell-type-specific enrichment that was confined to different segments of the
  proximal tubule". Also: "SLC5A11 showed high, specific expression in the S3 segment of the proximal tubule".
  - This is the only GWAS paper here that makes a PT-subsegment claim.
- Tin 2019: the proximal tubule is "the major site of urate reabsorption and excretion".
- Stanzick 2021: enrichment "for three proximal tubule clusters" (PT level only).
- Sheng 2021: "We found the strongest enrichment for kidney function (eGFR) and serum urate level in proximal
  tubules." It has one PT cluster and no S1/S2/S3.

## 3 · Where along the PT these genes act

**Anchors.**
- Rat 14-segment RNA-seq (Lee 2015), RPKM in S1/S2/S3:
  - SGLT2 90.9/8.2/0.5;
  - Slc7a13 0.3/137/780;
  - Slc3a1 50/114/710;
  - Slc22a12 1.6/21.6/120;
  - Slc22a8 80.7/152.3/0.7;
  - Oct2 31.1/96.5/293.7;
  - Lrp2 19/34/15;
  - Cubn 8.9/7.7/8.9.
  - Zeros can be annotation dropouts.
- Limbutara 2020 maps S1 to the cortical labyrinth, S2 to the cortical medullary ray and S3 to the outer stripe.
- Lake 2023 (human):
  - curated PT-S3 markers: SLC22A7, MOGAT1, SLC5A11, SLC22A24, SLC7A13, SLC5A8, ABCC3, SATB2;
  - PT-S1 and PT-S2 "could not be accurately resolved" and were merged.

**Urate.**
- Human snRNA per segment (Sakaguchi 2024):
  - SLC22A12 0.36/0.92/1.42;
  - SLC2A9 4.68/4.53/7.35;
  - SLC17A1 falls slightly along the PT;
  - ABCG2 is "rarely expressed".
- Halperin Kuhns 2021 puts ABCG2, SLC2A9 and SLC22A12 highest in S3. It conflicts with Sakaguchi on OAT1 and OAT3.
- Hoque 2020 (mouse): "neither ABCG2 nor URAT1 strongly localize in the early proximal tubule (S1)".
- Preitner 2009: mouse Glut9 is in the DCT and CNT, "no detection in the PCT". Mice have hepatic uricase, so urate is
  30–50 µM against 240–360 µM in humans.
- SLC22A11 (OAT4) has no mouse ortholog.

**Nephrotoxicity.**
- **Cisplatin injures S3 in rat** (Dobyan 1980, abstract only).
  - Uptake is through OCT2; mice need an Oct1/Oct2 double knockout (Filipski 2009).
  - Rat Oct2 rises toward S3. Human OCT2 is described in cortex and medullary ray (Motohashi 2002, abstract only).
  - So "OCT2 explains S3" is rodent-supported and unproven in human.
- **Aminoglycosides:** megalin, S1 and S2. Labelling density is "markedly higher in segments 1 and 2" in rat, rabbit
  and human (Christensen 1995, abstract only), but rat Lrp2 mRNA is nearly flat.
- **Tenofovir:** OAT1 uptake and MRP4 efflux (Kohler 2011). Human OATs are stained "more intensely in S1/S2
  segments … whereas the S3 segment in the outer stripe remained unstained" (Breljak 2016, abstract only).
- **Mercury:** injury in S2 and S3. Losing Mrp2 moves the lesion to S1 (Zalups & Bridges 2014).
- **Contrast media:** the primary target is the medullary thick limb (Heyman 1991, abstract only).

**Fanconi and the S1 load.**
- "The S1 segment possesses the highest capacity for sodium, solute, amino acid, and fluid transport among all renal
  tubular segments" (Zhuo & Li 2013). No source gives a numeric S1 against S2 split except glucose: SGLT2 "accounts
  for 97% of the net glucose reabsorption" (Curthoys & Moe 2014).
- Gene-level evidence:
  - GATM is in "the early segments" (Reichold 2018).
  - EHHADH is in S2 and S3 (Forst 2021).
  - NBCe1-A is highest in S1 (Skelton 2010).
  - rBAT (SLC3A1) is highest in S3 and b0,+AT (SLC7A9) in S1 (Nagamori 2016, abstract only).

## 4 · Notes for protocol design

- **Declare the structural species asymmetries up front:**
  - SLC22A11 has no mouse ortholog;
  - mouse Slc2a9 is distal;
  - ABCG2 is protein-abundant but transcript-scarce in human;
  - human MATE2-K has no separate mouse gene.
- **The "S2" boundary differs between rodent microdissection (medullary ray) and human histology.** Some apparent
  species disagreement is nomenclature.
- **ChEMBL drug mechanisms per transporter:**
  - SLC5A2: gliflozins;
  - SLC22A12: uricosurics;
  - SLC22A6: probenecid.

**Abstract-only papers that carry load:** Motohashi 2002, Christensen 1995, Dobyan 1980, Breljak 2016, Karbach 2000,
Nagamori 2016 and Gopal 2007.
