# Figs. 5 and 6: which results to highlight (2026-11-09)

The user asked for the best results to show in the last two figures, chosen separately for mouse vs human (Fig. 5)
and control vs AKI (Fig. 6). Preference goes to clinically relevant genes, genes the literature discusses, and drug
targets. Subtle transporter handoffs are an example of the kind of result wanted.

**Rule.** Every candidate comes from a calibrated screen call: `fig5-human-vs-mouse-along-pt.md`,
`fig6-aki-along-pt.md`, or the one new read below. Clinical relevance only ranks calls the screens already made; it
never creates one.

**Annotation** (`analysis/paper/fig56_annotate.py`, run 2026-11-09, network):
- **Human ortholog:** Ensembl Compara 98, preferring one-to-one.
- **Drugs:** Open Targets drug and clinical candidates, giving the maximum clinical stage.
- **Literature:** Europe PMC co-mentions.
- **Disease genes:** PanelApp green tubulopathy and stone genes.
- **GWAS loci:** eGFR, urate and urine-metabolite loci, the gene sets of notebook 108.
- **Kidney literature count:** PubMed `(SYMBOL or alias)[tiab] AND (kidney or renal or proximal tubule)[tiab]`, for
  the shortlist only.

Coverage: 8,258 of 8,461 genes have a human ortholog. 605 are drug targets, 321 of them with an approved drug.

**Ranking:** `analysis/paper/fig56_rank.py`. Its tables are in `results/paper_figures/highlights/`.

**Drug annotations need judgment.** Several approved-drug hits are irrelevant to the kidney:
- Chrna4 (varenicline);
- Odc1 (eflornithine);
- Apoc3 (volanesorsen): the human PT does not express it, so it is an artefact.

## New read for Fig. 6: the injury response differs between S1 and S3 (post hoc, calibrated)

The response is the IR − control contrast of the segment means (label level, not the curve). A gene is called when
all of these hold:
- the S3 response minus the S1 response is at least 1 log2;
- all four IR × control pairs agree in sign, each by at least 0.5;
- both relabelings are at most 0.25 × the condition.

**Result:** 268 genes are called, against 0 and 4 under the relabelings (false-call share 0.007).

This read uses segment labels, not the continuous axis. It answers where injury hits, not whether the curve sees
more than clusters, so it belongs in Extended Data or in text unless the figure needs a clinical anchor.

## Fig. 5 — mouse vs human

**Reversed gradients** (category R calls; the external verdicts are from the full-text checks in the Fig. 5 doc):

| Gene | Protein / clinical hook | Mouse E | Human E | Kidney papers | External |
|---|---|---|---|---|---|
| MME | neprilysin; sacubitril (ARNI) target; degrades natriuretic peptides | +3.1 | −0.9 | 733 | both limbs reproduced |
| SLC17A3 | NPT4 urate exporter; urate and eGFR GWAS locus | +1.1 | −0.9 | 36 | both limbs reproduced |
| HPGD | 15-PGDH, inactivates PGE2; 15-PGDH inhibitors in development for regeneration | +2.4 | −1.5 | 32 | both limbs reproduced |
| SLC5A11 | SMIT2, myo-inositol uptake; urine-metabolite GWAS locus | −1.1 | +1.0 | 14 | both limbs reproduced (human-late half new) |
| DCXR | L-xylulose reductase; pentosuria | −1.0 | +1.2 | 9 | both limbs reproduced |
| SORBS1 | eGFR GWAS locus | −0.7 | +0.8 | 5 | both limbs reproduced |

**Recommended for panel a:** Mme, Slc17a3, Hpgd, Slc5a11 and Dcxr. Mme, Slc17a3 and Slc5a11 carry a drug, GWAS or
urate hook, and all five reproduce externally. This would replace Lpar6 and Pitpnc1, which have 2 and 0 kidney
papers and no clinical hook. Sorbs1 is the spare.

**Transporter handoffs that move**, each a drug or disease transporter. These are the subtle positional results.

| Pair | Δx\* (human − mouse) | Clinical hook | Human atlas |
|---|---|---|---|
| SGLT2 → SGLT1 | +0.30 | gliflozins (SGLT2); sotagliflozin (SGLT1/2); SGLT1 takes over glucose uptake under SGLT2 block | reproduces |
| SLC22A23 → URAT1 | +0.25 | uricosurics (dotinurad, verinurad, lesinurad); urate GWAS | reproduces |
| y+LAT1 → AGT1 | +0.24 | SLC7A7: lysinuric protein intolerance (tubulopathy panel); AGT1: cystine reabsorption | reproduces |
| SLC22A23 → OCT2 | +0.24 | OCT2: metformin and cisplatin uptake, a regulatory drug–drug interaction transporter | reproduces |
| SGLT2 → SGLT5 | +0.21 | gliflozin context (SGLT5: mannose and fructose) | reproduces |
| B0AT1 → B0AT3 | +0.18 | SLC6A19: Hartnup (tubulopathy panel); B0AT1 inhibitors pursued for PKU (*verify trial stage before citing*) | reproduces |
| MCT1 → MCT9 | −0.49 | MCT9 runs opposite in the two species (a reversal, not a moved switch) | reproduces |

Not recommended: OCTN2 → URAT1 and OCTN2 → OAT2. Both shift in our data, but the human atlas crosses later, the
opposite direction.

**Recommended for panel c:** SGLT2 → SGLT1, SLC22A23 → URAT1, SLC22A23 → OCT2, y+LAT1 → AGT1 and B0AT1 → B0AT3.
This is close to v2. It swaps MCT1 → MCT9 for OCT2, a drug–drug interaction transporter, and keeps MCT9 in the text
as a reversal.

**Optional anchor: drug targets that keep their position** (descriptive). Same sign and |E| ≥ 0.5 in all four
sections: SGLT2, URAT1, OAT3 (probenecid), ACE (ACE inhibitors, S3-high in both), PSMA/FOLH1 (¹⁷⁷Lu-PSMA
radioligands; renal uptake is dose-limiting), DPEP1 (cilastatin), PAH, GATM, GLUT2, NBCe1 (SLC4A4). These give the
contrast "most targets keep their place, these move". Amplitudes are not comparable (human flattening), so show
directions only.

"Mouse gradient, human flat" (52 clinical genes, e.g. PPARA, CPT2, GHR) is not shown: flattening of the human arm can
produce it.

## Fig. 6 — control vs AKI

**Handoffs that move upstream after IRI** (10 after the neighbour filter, 0 under the relabelings):

| Pair | Δx\* | Clinical hook |
|---|---|---|
| OAT1 → OCT2 | −0.14 | the two main PT drug secretors; OAT1: probenecid, tenofovir and adefovir uptake; OCT2: cisplatin and metformin uptake |
| 4F2hc → rBAT | −0.20 | SLC3A1: cystinuria (tubulopathy panel) |
| B0AT1 → GAT2 | −0.23 | SLC6A19: Hartnup |
| Slc46a3 → PCFT | −0.22 | PCFT: folate and antifolate (pemetrexed) uptake |
| MCT1 → MCT5 / MCT9 | −0.42 / −0.30 | MCT1 induced about 27× in injured S3 (MCT1 inhibitor AZD3965 in oncology trials) |
| ZIP5 → ZIP8 / ZIP14 | +0.17 / +0.20 | ZIP8 knockout protects in IRI; ZIP8 and ZIP14 take up cadmium |

**Recommended for panel d:** OAT1 → OCT2, 4F2hc → rBAT, B0AT1 → GAT2, PCFT and ZIP5 → ZIP8, with cystine
b0,+AT → AGT1 as the unshifted internal control. OAT1 → OCT2 is the clinically strongest: injury shifts the secretory
handoff that sets nephrotoxic-drug exposure.

**Gradients that injury builds or erases** (named, checked examples; Category A itself is not calibrated):

| Gene | Event | Clinical hook | Kidney / AKI papers |
|---|---|---|---|
| HAVCR1 (KIM-1) | built: climbs through S2 (+1.0 vs −0.3 in controls) | qualified urinary biomarker of tubular injury | 3,127 / 1,746 |
| CYP24A1 | erased in S2; induced far more in S3 than S1 (+7.7 vs +1.7) | vitamin D catabolism, FGF23 axis; tubulopathy panel (infantile hypercalcaemia) | 339 / 7 |
| AGT | erased in S2; falls in S3 (−2.5) | intrarenal RAS; urinary AGT studied as an AKI/CKD marker | 1,902 / 193 |
| CYP4A10/14 (human CYP4A11) | erased in S2 | 20-HETE synthesis (PPARα programme) | 68 / 16 |
| Respiratory chain set | built: falls through S3 | mitochondrial injury | — |
| DNA replication set | built: climbs through S2 | repair / proliferation | — |

**Recommended for panels a–b:** KIM-1, the replication set and oxidative phosphorylation as built gradients;
Cyp24a1, Agt and Cyp4a14 as erased ones. This is the current figure, minus Cyp51 (26 kidney papers; its direction
rests on our own reanalysis) and Anpep. KIM-1 is the anchor: it is the clinical injury biomarker, and here it is graded
within S2.

**From the new segment read** (text or Extended Data; drug targets whose injury response differs between S1 and S3):

| Gene | S1 / S2 / S3 response | Hook |
|---|---|---|
| XDH | +1.3 / +2.1 / +2.6 | xanthine oxidase: allopurinol, febuxostat; reperfusion ROS source |
| CDK4, CDK2, E2F1 | +0.9 / +1.4 / +1.9 (CDK4) | CDK4/6 inhibitors studied as AKI protectants (*verify sources*) |
| TYMS, GART | +1.1 / +1.8 / +2.6 (TYMS) | antifolate targets (pemetrexed, 5-FU) |
| HSD11B1 | −0.7 / −0.9 / −3.1 | 11β-HSD1 inhibitors; local glucocorticoid activation |
| SLC3A1 (rBAT) | −0.3 / −0.6 / −1.8 | cystinuria; matches the rBAT handoff |
| HNF1B | −0.2 / −0.8 / −1.5 | ADTKD / HNF1B nephropathy; master PT transcription factor |
| TTR | −1.5 / −1.2 / +1.4 | tafamidis, patisiran (amyloid) |

## What to drop from the candidate pools

- **Fig. 5:**
  - Apoc3: not expressed in human PT.
  - Chrna4 and Odc1: drug hits with no kidney relevance; Odc1's human limb is weak.
  - Hao2, Vnn1, Sh3gl2, Lrig1: androgen-dependent in mouse.
  - OCTN2 pairs: the human atlas disagrees.
- **Fig. 6:** Fgb and Serpina10 (control noise), Apob, Esr1 and Apoc3 (near the detection floor), and Slc5a9 → Slc5a3
  (neighbour signal).

## Open before the figures are final

- Read in full the sources for: B0AT1 inhibitors in PKU, CDK4/6 inhibition in AKI, 15-PGDH inhibitors, and urinary AGT
  in AKI (full-text rule).
- Every Fig. 5 statement is about one human donor.
