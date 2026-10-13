# Literature for non-1:1 orthology genes along the proximal tubule (pre-registration input)

Prepared from published literature only; no project data, results, notebooks or paper drafts were read.
Target setting: Visium HD probe-based tubule aggregates; continuous PT coordinate S1 -> S2 -> S3 (mouse late S3 grading
into thin-limb identity); mouse 2 healthy male C57BL/6 + 2 male IRI-AKI; human one healthy male cortex donor; separate
probe panel per species.

## How to read this file

- `verified` codes: `yes-full` = supporting sentence read in open full text; `yes-abs` = sentence is in the PubMed
  abstract; `yes-data` = value read from the paper's deposited data file (not article text); `pmid-only` = citation
  verified but claim not checkable in available text; `unverified` = could not verify.
- Every PMID in this file was checked against the PubMed record (title, year, journal). Quoted text is verbatim.
- Orthology classes come from the MGI/Alliance mouse-human homology report (`HOM_MouseHumanSequence.rpt`, downloaded
  2026-10-13). "class" = the MGI homology class of the gene itself; "FAMILY n:m" = the connected component of all mouse
  and human genes linked through any shared class. A gene can be 1:1 at class level and non-1:1 at family level
  (Havcr1, Ccl2, Slco1a1, Serpina3n). **Any analysis must state which of the two it filters on.**
- `strength`: strong = direct segment-resolved data; medium = PT-level or single study; weak = whole-kidney, indirect,
  or review.

## Orthology classes of the families in scope (MGI family-level components, mouse:human)

| family | class | mouse members (in report) | human members |
|---|---|---|---|
| Cyp4a | 8:2 | Cyp4a10, 4a12a, 4a12b, 4a14, 4a29, 4a30b, 4a31, 4a32 | CYP4A11, CYP4A22 |
| Cyp2d | 8:2 | Cyp2d9, 2d10, 2d11, 2d12, 2d22, 2d26, 2d34, 2d40 | CYP2D6, CYP2D7 |
| Cyp2j | 8:1 | Cyp2j5, 2j6, 2j7, 2j8, 2j9, 2j11, 2j12, 2j13 | CYP2J2 |
| Cyp4f | 5:4 | Cyp4f13, 4f14, 4f15, 4f16, 4f40 | CYP4F2, CYP4F8, CYP4F11, CYP4F12 |
| Ugt2b | 7:7 | Ugt2b1, 2b5, 2b34, 2b35, 2b36, 2b37, 2b38 | UGT2B4, 2B7, 2B10, 2B11, 2B15, 2B17, 2B28 |
| Ugt1a (1A7/8/9/10 arm) | 3:4 | Ugt1a7c, 1a8, 1a9 | UGT1A7, 1A8, 1A9, 1A10 |
| Ugt1a2/1a5 arm | 2:3 | Ugt1a2, Ugt1a5 | UGT1A3, 1A4, 1A5 |
| Ugt1a1 | 1:1 | Ugt1a1 | UGT1A1 |
| Akr1c | 9:4 | Akr1c6, 1c12, 1c13, 1c14, 1c18, 1c19, 1c20, 1c21, Akr1cl | AKR1C1-C4 |
| Ces1 | 8:2 | Ces1a-Ces1h | CES1, CES1P1 |
| Ces2 | 7:1 | Ces2a, 2b, 2c, 2e, 2f, 2g, 2h | CES2 |
| Nat8/Nat8f | 9:2 | Nat8, Nat8b-ps, Nat8f1-f7 | NAT8, NAT8B (NAT8B inactive in humans) |
| Slco1a | 6:1 | Slco1a1, 1a4, 1a5, 1a6, 1a7, 1a8 | SLCO1A2 |
| Slc22 Oat-like | 6:3 | Slc22a19, a26, a27, a28, a29, a30 | SLC22A9, SLC22A10, SLC22A25 |
| Slc22a22 | mouse-only | Slc22a22 | - |
| SLC22A11 (OAT4), SLC22A24 | human-only | - | no MGI class at all |
| Abcb1 | 2:1 | Abcb1a, Abcb1b | ABCB1 |
| Gsta | 5:4 | Gsta1, Gsta2, Gsta3, Gsta5, Gsta13 | GSTA1, A2, A3, A5 |
| Gsta4 | mouse-only in the report | Gsta4 | (human GSTA4 exists but is unclassed here - see caveat) |
| Gstm | 6:4 | Gstm1, 2, 3, 5, 6, 7 | GSTM1, M2, M3, M5 |
| Sult1a | 1:4 | Sult1a1 | SULT1A1, 1A2, 1A3, 1A4 |
| Sult1d1 | 1:1 (human is pseudogene SULT1D1P) | Sult1d1 | SULT1D1P |
| Sult2a | 8:1 | Sult2a1-Sult2a8 | SULT2A1 |
| Mt | 2:10 | Mt1, Mt2 | MT1A, B, E, F, G, H, L, M, X, MT2A |
| Serpina3 | 10:1 | Serpina3a, b, c, f, g, i, j, k, m, n | SERPINA3 |
| Serpina1 | 5:1 | Serpina1a-1e | SERPINA1 |
| Mup | mouse-only / unclassed | Mup3, Mup20 mouse-only; Mup1 not in report | - |
| Kap, Slc7a12 | mouse-only | Kap; Slc7a12 | - |
| Sprr2 | 6:6 | Sprr2b, d, e, f, g, i (Sprr2f class mouse-only) | SPRR2A, B, D, E, F, G |
| Saa | 3:3 | Saa1, Saa2, Saa3 | SAA1, SAA2, SAA3P (human SAA3P pseudogene) |
| Apol | 11:4 | Apol7a/b/c/e, Apol8, Apol9a/b, Apol10a/b, Apol11a/b | APOL1-APOL4 |
| H2/HLA class I | 15:7 | H2-K1, D1, Q*, T23, M* | HLA-A, B, C, E, F, G, H |
| Akr1b7/b8 | 3:2 | Akr1b7, Akr1b8, Akr1b10 | AKR1B10, AKR1B15 |
| Aox | Aox1 1:1, Aox3 1:1 (human AOX3P pseudogene), Aox4 mouse-only | Aox1, Aox3, Aox4 | AOX1, AOX3P |
| Havcr1 | class 1:1, FAMILY 4:1 | Havcr1, Timd2, Timd5, Timd6 | HAVCR1 |
| Rdh16 group | 6:1 | Rdh1, 7, 9, 16, 16f2, 19 | RDH16 |
| Acsm2 | 1:2 | Acsm2 | ACSM2A, ACSM2B |
| 1:1 (controls / landmarks) | 1:1 | Slc5a1, Slc5a2, Slc5a12, Slc34a1, Slc22a6, Slc22a7, Slc22a8, Slc22a12, Slc22a13, Slc7a13, Cyp2e1, Cyp7b1, Cyp4b1, Lrp2, Prlr, Acsm3, Abcc3, Abcc4, Fmo1, Fmo2, Fmo5, Lcn2, Krt20, Vcam1, Spp1, Sox9, Atf3, Hnf4a, Cdh6 | matching human symbol |

Caveats on the orthology source: (a) `Gsta4` has no mouse-human class in the MGI report although a human GSTA4 exists,
so "mouse-only" here means "unclassed in this report", not proven absence; (b) `Mup1` and some `Sprr2a*` genes are
likewise absent from the report; (c) HGNC lists mouse links for SLC22A25, APOL1, CYP4A22, UGT2B17, NAT8B but none for
GSTA4, SLC22A11 or SLC22A24, which agrees with the human-only calls for OAT4 and SLC22A24.

---

# 1. Gene-level positional claims

## 1A. Published positional claims (literature)

| id | gene(s) | species | sex | claim | direction | evidence type | PMID | verified | supporting text | strength |
|---|---|---|---|---|---|---|---|---|---|---|
| L01 | Kap | mouse | male (intact + castrate) | KAP mRNA in PT epithelium; in castrates only the juxtamedullary subpopulation keeps it, androgen recruits cortical PT | S3/juxtamedullary constitutive; S1/S2 androgen-inducible | ISH | 3153475 | yes-abs | "only a subpopulation of these cells, those located in the juxtamedullary region of the renal cortex, contain KAP mRNA in castrated males. After androgen treatment, cortical proximal tubules are also induced" | strong |
| L02 | Kap | mouse | male, castrate, female | two spatial compartments: whole-cortex PT under androgen, outer-stripe only after castration | S1/S2 androgen-dependent; S3 androgen-independent | ISH + Northern | 2739660 | yes-abs | "It is expressed throughout the cortex in the epithelial lining of the proximal tubules in response to androgen stimulation. After castration, only tubules in the outer stripe of the medulla express KAP mRNA" | strong |
| L03 | Kap | mouse (Tfm/Y, WT, ovx) | male + female | in AR-deficient mice Kap is S3-restricted; female S3 expression survives ovariectomy | S3 only without AR | ISH + hormone manipulation | 2293028 | yes-abs | "KAP mRNA was detected ... only in the S3 segment of the proximal tubules in the outer medulla"; "ovariectomy of normal female mice did not eliminate KAP gene expression in the S3 cells" | strong |
| L04 | Kap | mouse (hypox, dw/dw) | both | S3 component is pituitary/estrogen-controlled, cortical component androgenic; TSH gives the strongest S3 response | two regulatory domains along PT | ISH after hypophysectomy | 1338721 | yes-abs | "epithelial cells of proximal convoluted tubules of cortical nephrons express KAP mRNA in response to androgenic stimulation while similar cells in the juxtamedullary S3 segment ... express KAP mRNA under estrogenic and pituitary hormonal control" | strong |
| L05 | Kap | mouse (hyt/hyt) | not stated | thyroid axis sustains the outer-stripe (S3) component, androgen the cortical one | S3 != cortical regulation | ISH + TH replacement | 7956935 | yes-abs | "TSH proved to be necessary for the constitutive expression of the gene in the outer stripe of the outer medulla, androgens are responsible for expression in cortical segments of the proximal tubules" | strong |
| L06 | Kap | mouse | both | explicit pars recta vs pars convoluta split of Kap regulation | PR (S3) hormone-responsive; PCT (S1/S2) strictly androgenic | genetic + promoter | 16150866 | yes-abs | "Whereas the pars recta (PR cells) responds to thyroid and sexual hormones, the pars convoluta (PCT cells) represents a truly androgen-dependent compartment" | strong |
| L07 | Kap | mouse (WT, APRT-/-) | male + female | baseline Kap is an S3 gene, strongly male-biased | S3 > S1/S2; male >> female | ISH / RT-PCR | 12486089 | yes-abs | "In 1-month-old wild-type male mice, Kap was detected primarily in S3 proximal tubule segments, but expression was very low in female mice" | strong |
| L08 | Cyp4a10, Cyp4a12, Cyp4a14 | mouse (C57BL/6J) | both | microdissected PT: males express Cyp4a10 + Cyp4a12; females Cyp4a10 + Cyp4a14; CYP4A/4F in all segments, sex effect in PT and glomeruli | Cyp4a12 male-only, Cyp4a14 female-biased | microdissection RT-PCR | 12388424 | yes-abs | "In the proximal tubular segments from male mice, the 4A10 and -12 isoforms were present, whereas the 4A10 and -14 isoforms were detected in segments from female mice." | strong for sex; no S1/S2/S3 split |
| L09 | Cyp4a12a, Cyp4a10, Cyp4a14, Cyp4a12b | mouse (5 strains) | both | whole kidney: Cyp4a12a male-specific and the dominant 20-HETE synthase; females high Cyp4a10/14; Cyp4a12b ~undetectable; DHT induces Cyp4a12a >4-fold | Cyp4a12a male >> female | qPCR + immunoblot + activity | 17112342 | yes-full | "All female mice showed low AA hydroxylase activities ... and very low Cyp4a12a mRNA and protein levels, but high Cyp4a10 and Cyp4a14 expression."; "Treatment with 5alpha-dihydrotestosterone induced both 20-HETE production and Cyp4a12a expression more than 4-fold in male C57BL/6 mice." | medium (whole kidney) |
| L10 | Cyp4a14, Cyp2e1 | mouse (C57BL/6J) | 2M + 2F | RNAscope with Slc34a1: Cyp4a14 female-only and Cyp2e1 male-only inside the S2 subset; S2 and S3 sex domains abut sharply | both S2-restricted | multiplex ISH + scRNA | 31689386 | yes-full | "sex-restricted activity of Cyp2e1 (male only) and Cyp4a14 (female only) in the S2 subset of Slc34a1+ PT regions (Figure 2I)" | strong for S2; see conflict with L41 |
| L11 | Cyp4a2, Cyp4a3, Cyp4a8; CYP4A protein | rat | both | 4A2/4A3/4A8 mRNA in both PCT and PST; CYP4A protein strongest in PT; 4A2 4-fold male > female | PT-enriched | microdissection RT-PCR + IHC | 9486235 | yes-abs | "P-450 4A2 mRNA was readily detected in preglomerular arterioles, glomeruli, proximal convoluted tubule (PCT), proximal straight tubule (PST)"; "the expression of 4A2 mRNA was fourfold greater in the kidney of male vs. female rats" | medium |
| L12 | CYP4A11, CYP4F2 | human | not stated | immunocytochemistry: both proteins in **only S2 and S3**; CYP2C8/2C9/2E1 not in human kidney | S2+S3 > S1 (S1-absent) | IHC + immunoinhibition | 10660572 | yes-abs | "Immunocytochemistry showed CYP4F2 and CYP4A11 expression in only the S2 and S3 segments of proximal tubules in cortex and outer medulla." | strong (best human positional datum here) |
| L13 | CYP4A11 vs CYP4A22 | human | not stated | two 96%-identical CYP4A genes present in 100/100 individuals; CYP4A11 predominant in kidney | CYP4A11 >> CYP4A22 | Southern + RT-PCR | 12464262 | yes-abs | "These two genes share 96% sequence identity"; "RT-PCR amplification of human kidney RNA ... showed that CYP4A11 is the predominant isoform expressed in kidney" | medium |
| L14 | CYP4A22 | human | not stated | lower-expressed, lower-activity paralog; the reference sequence differs from the common allele | CYP4A22 low | qPCR (12464261); enzymology (31199497) | 12464261 | yes-abs | "CYP4A22 mRNA is expressed at significantly lower levels than CYP4A11 mRNA in human liver samples" | medium |
| L15 | CYP2D6 / mouse Cyp2d | human / humanised mouse | not stated | organ and cell-level only; **no PCT/PST or S1/S2/S3 assignment found anywhere** | no segment-resolved data found | immunoblot + immunocytochemistry | 16033950 | yes-abs | "The cell-specific expression of human CYP2D6 in liver, kidney, and intestine in humanized mice was similar to that reported in humans." (no renal segment named) | weak / negative |
| L16 | Cyp2j5 (CYP2J protein) | mouse | not stated | CYP2J protein primarily in proximal tubules and collecting ducts; Cyp2j5 mRNA in cortex and outer-medulla tubules | PT + CD, no subsegment | IHC + ISH + Northern | 10364221 | yes-abs | "CYP2J protein(s) are present primarily in the proximal tubules and collecting ducts"; "In situ hybridization confirms abundant CYP2J5 mRNA within tubules of the renal cortex and outer medulla." | medium |
| L17 | Cyp2j5 | mouse (+ castration/DHT/E2/Tfm/aERKO) | both | renal CYP2J5 higher in males after puberty; androgen up, estrogen down | male > female | Western + Northern | 14978252 | yes-abs | "after puberty, CYP2J5 protein is present at higher levels in male mice than in female mice"; "up-regulated by androgen and down-regulated by estrogen" | medium |
| L18 | Cyp2j5, Cyp2j13, Cyp2j11 | mouse | both | Cyp2j5 and Cyp2j13 high in male kidney; Cyp2j11 in both sexes | Cyp2j5/13 male-biased; Cyp2j11 flat | qPCR tissue panel | 25994032 | yes-abs | "Cyp2j5 and Cyp2j13 were highly expressed in male kidneys, and Cyp2j11 was highly expressed in both male and female kidneys." | medium |
| L19 | Cyp2j13 | mouse (C57BL/6) | M, F, nephron AR-KO | PT-S2, male, AR-dependent, colocalised with Cyp2e1; its promoter AR peak is the one most reduced after AR removal | S2, male | snRNA/snATAC + RNAscope | 37673062 | yes-full | "M-biased gene Cyp2j13 ... exhibited AR-dependent expression in PT-S2, colocalized with S2 marker Cyp2e1"; "out of the 6 M-biased DARs (5 intronic and 1 in the promoter) with AR binding annotated in the genomic region of Cyp2j13, the promoter region peak was identified as the most significantly down-regulated post AR removal (FDR<3.4E-5)" | strong |
| L20 | CYP2J2 | human | not stated | tissue-microarray IHC: wide distribution, **no renal segment assignment** | no segment-resolved data found | IHC microarray | 15033996 | yes-abs | "CYP2C8 and 2J2 were also widely distributed in human tissues but were less frequently associated with sEH." | weak / negative |
| L21 | Cyp4f13-4f18 | mouse (C57BL/6J) | both | CYP4F isoforms in every nephron segment; the CYP4F sex difference is in **glomeruli**, not PT | no PT sex difference | microdissection RT-PCR | 12388424 | yes-abs | "In glomeruli, sex differences in the expression pattern of CYP4F isoforms were also observed, with male mice expressing the 4F13, -14, and -15 isoforms, whereas female mice expressed the 4F13, -16, and -18 isoforms." | medium-strong |
| L22 | Ugt2b5/37/38, Ugt1a2, Ugt1a7c | mouse (C57BL/6) | both | kidney-expressed Ugt set; Ugt2b5/37/38 male-predominant, Ugt1a2 female-predominant (whole kidney) | as stated | bDNA mRNA, 14 tissues | 17050650 | yes-abs | "In kidney, Ugt1a2, Ugt1a7c, Ugt2b5/37/38, Ugt2b35, and Ugt3a1/2 were expressed."; "Male-predominant expression was observed for ... Ugt2b5/37/38 in kidney" | weak for position |
| L23 | Ugt2b38, Ugt1a2 | mouse (C57BL/6, gonadectomy/hypox/GH models) | both | the kidney Ugt sex difference is androgen-driven, unlike liver (GH-pattern) | Ugt2b38 male (T-induced); Ugt1a2 female (T-repressed) | bDNA + endocrine surgery | 19131521 | yes-abs | "In kidney, regulation of Ugt2b38 expression was attributed to inductive effects by testosterone. Conversely, Ugt1a2 expression in kidney was negatively regulated by testosterone." | medium |
| L24 | UGT1A (pan), UGT2B7 | human | not stated | IHC: UGT1A and UGT2B7 throughout proximal **and** distal tubules, loops, CD; glomerulus negative; **no S1/S2/S3 gradient reported** | PT-positive, flat | IHC + microsomal kinetics | 17698974 | yes-abs | "With the exception of the glomerulus, Bowman's capsule, and renal vasculature, UGT1A proteins and UGT2B7 were expressed throughout the proximal and distal convoluted tubules, the loops of Henle, and the collecting ducts" | medium |
| L25 | UGT1A locus | human | n/a | one locus, unique exon 1 per isoform and **4 shared 3' exons** | isoform assignment risk | genomic | 1339448 | yes-abs | "The 3' end of the locus contains 4 common exons which encode the identical carboxyl termini (246 amino acids)" | strong (methodological) |
| L26 | UGT1A9, UGT2B7, CYP4F2 (+CYP2B6, CYP3A5) | human cortex | not stated | targeted proteomics: UGT1A9 and UGT2B7 most abundant, CYP4F2 low; no subsegment | abundance only | quantitative proteomics | 34390491 | yes-abs | "CYP2B6, CYP3A5, and CYP4F2 showed comparable, but generally low expression, whereas UGT1A9 and UGT2B7 levels were the highest." | weak for position |
| L27 | Akr1c18 | mouse (C57BL/6J) | 2M + 2F | Akr1c18 among female-enriched **S3** PT genes | female > male within S3 | zone-dissected scRNA | 31689386 | yes-full | "Additional female S3 enriched genes are linked to hormone metabolism (including Dio1, Ttr, Akr1c18, Gm4450, Spp1 and Rdh16)" | medium |
| L28 | AKR1C9 (rat), AKR1A3, AKR1B4 | rat | not stated | IHC along nephron: AKR1C9 in the **pars recta (S3)**, AKR1A3 in PCT, AKR1B4 in loop of Henle | AKR1C9 S3-specific | IHC + Western | 19920056 | yes-abs | "the localization of AKR1A3 in proximal convoluted tubules, AKR1B4 in the loop of Henle, and AKR1C9 in the pars recta S3 segment of proximal tubules" | strong (rat protein) |
| L29 | Akr1c21 | mouse (CD1) | 3M, 3F, 3 pregnant F | Akr1c21 detected only in kidney among mouse tissues; human and mouse AKR1C enzymes diverge in phylogeny, substrate and tissue | kidney-restricted; no functional 1:1 map | RT-PCR + enzymology | 20003443 | yes-full | "Tissue specificity was also observed for Akr1c21, detected only in the kidney." | medium + key orthology caveat |
| L30 | Akr1a4, Akr1c18, Akr1c21, Akr7a5 | mouse | not stated | the kidney-predominant members of a 16-Akr panel | kidney-predominant set | qPCR tissue panel | 23660342 | yes-abs | "Akr1a4, 1c18, 1c21, and 7a5 in kidney" | weak |
| L31 | AKR1C1-AKR1C4 | human | not stated | 4-gene cluster >=84% identity; explicit contrast with the single rodent 3alpha-HSD; AKR1C4 liver-specific, AKR1C3 prostate/mammary, AKR1C1/C2 brain; **kidney not in the tissue list** | no segment-resolved data found | isoform-specific RT-PCR | 10998348 | yes-abs | "share at least 84% amino acid sequence identity"; "These studies are in stark contrast with those in rat where only a single AKR with positional- and stereo-specificity for 3alpha-hydroxysteroids exists" | weak / negative for kidney |
| L32 | AKR1C1 | human | not stated | Western across 10 tissues: AKR1C1 in only kidney, liver, testis; AKR1 members in proximal and distal convoluted tubules by IHC | PT-positive, no subsegment | Western + IHC | 10510318 | yes-abs | "AKR1C1 being observed in only kidney, liver and testis"; "AKR7A2 is expressed in a similar fashion to the AKR1 family members in proximal and distal convoluted renal tubules" | medium |
| L33 | CES2 | human | not stated | tissue-array IHC: highest renal CES2 in the **proximal convoluted tubule** | PCT (S1/S2) | IHC tissue array | 12607608 | yes-abs | "The highest carboxylesterase-2 expression occurs in hepatocyte, small intestine mucosa, kidney proximal convoluted tubule, and adrenal cortex cells." | medium-strong |
| L34 | Ces1 ("hydrolase B") | rat | not stated | the renal carboxylesterase identical to liver hydrolase B localises to proximal tubules | PT, no subsegment | cloning + immunocytochemistry | 7961958 | yes-abs | "hydrolase B is localized in the centrilobular region of the liver and in the proximal tubules of the kidney" | medium; no human CES1 segment data found |
| L35 | Fmo1-Fmo5 | mouse (129/SV) | both, several ages | ISH: every Fmo mRNA in distal **and** proximal tubules and CD (Fmo1 also glomerulus); kidney Fmo5 4-fold male > female at all ages | PT + DT + CD; Fmo5 male-biased | ISH + RNase protection | 15183119 | yes-abs | "In kidney, each FMO mRNA is localized to the distal and proximal tubules and collecting ducts; FMO1 mRNA is present also in the glomerulus."; "FMO5 mRNA abundance is fourfold greater in males than in females, at all stages of development." | medium-strong |
| L36 | Akr1b7 | mouse | not stated | renal AKR1B7 marks juxtaglomerular **renin cells**, not PT | PT-absent | IHC + lineage tracing | 26180185 | yes-abs | "throughout kidney development, AKR1B7 expression-together with renin-is progressively restricted in the kidney arteries toward the glomerulus" | strong for cell type |
| L37 | Akr1b8 / AKR1B10 | mouse / human | n/a | mouse Akr1b8 is the functional ortholog of human AKR1B10 (confusingly, a different mouse gene is named Akr1b10) | orthology caveat | phylogeny + enzymology | 21087085 | yes-abs | title and abstract: "identification of AKR1B8 as an ortholog of human AKR1B10" | medium |
| L38 | Gst isoforms (incl. Gsta/Gstm) | mouse | both | in whole kidney most Gst isoforms higher in females | female > male | bDNA, 19 Gsts x 14 tissues | 17890767 | yes-abs | "With a few exceptions, most of the Gst isoforms expressed in kidney showed higher expression in females than males" | weak for position |
| L39 | GST alpha class; GSTM3 | human | 8M + 8F (9622847) | class-alpha GST stains intensely in **kidney proximal tubular cells** (pi is distal/urinary tract); in cortex alpha is predominant, GSTM3 the main mu isoform | alpha PT-restricted | IHC (8372055) + activity/IHC (9622847) | 8372055 | yes-abs | "Intensive staining was visible in all hepatocytes, in kidney proximal tubular cells ..."; companion (9622847): "GSTM3 was the predominant isoenzyme of the mu class" | medium-strong |
| L40 | GST alpha vs pi along the human nephron | human | not stated | alpha solely in PT, pi restricted to DCT, thin loop and CD | segment contrast | IHC | 7843118 | yes-abs | "The alpha isoenzyme was present solely in the proximal tubular cells of the human kidney, while the pi form was restricted to the distal convoluted tubules, the thin loop of Henle, and the collecting ducts" | strong |
| L41 | Sult1d1, Sult1a1 | mouse | both | across 14 tissues Sult1d1 is kidney-highest; Sult1a1 highest in large intestine | Sult1d1 kidney-predominant | bDNA, 11 Sults | 16807285 | yes-abs | "Sult1b1, 2b1, and PAPSs2 in small intestine; Sult1a1 in large intestine; Sult1c2 in stomach; Sult1d1 in kidney" | medium for organ |
| L42 | SULT1A1, SULT1A3/4, SULT2A1 etc | human | not stated | quantitative immunoblot: "The kidney and lung contained low levels of SULT"; SULT1A3 and 1A4 quantifiable only jointly | kidney-low; no segment | quantitative immunoblot | 19679676 | yes-full | "The kidney and lung contained low levels of SULT." | weak / partial negative |
| L43 | NAT8 / NAT8B | human | n/a | NAT8B is inactive in humans (premature stop at codon 16); NAT8 almost exclusive to kidney + liver; no segment data | human has one functional gene | genomic + enzymology | 20392701 | yes-abs | "disclosed the systematic presence of a premature stop codon at codon 16"; "The almost exclusive expression of NAT8 in kidney and liver" | strong for dosage; none for segment |
| L44 | NAT8 (TSC501) | human | not stated | kidney- and liver-specific expression, no nephron resolution | tissue-level only | differential display | 9852678 | yes-abs | "we isolated a novel human gene, designated TSC501, that is expressed in kidney and liver" | weak |
| L45 | Nat8 (mouse) | mouse (C57BL/6J) | male | Nat8 is kidney-predominant among liver/kidney/intestine; no segment resolution | kidney > liver | RNA-seq + qPCR | 23382457 | yes-full | "Nat8 was predominant in kidney (Fig. 10A)" | medium for tissue |
| L46 | Nat8f1-f7 | mouse | - | **no segment-resolved data found for any Nat8f paralog in any species**; the family is a rodent tandem expansion | n/a | - | 20392701 | yes-abs | "NAT8 homologues ... are often encoded by multiple, tandemly repeated genes as many other genes encoding xenobiotic metabolism enzymes" | n/a |
| L47 | Slco1a1 (rat oatp1) | rat | sex not stated for IHC | protein on the **apical** membrane of S3 in the outer stripe; kidney mRNA concentrated in outer medulla | S3 > S1/S2 | immunoblot + immunomorphology | 8770038 | yes-abs | "apical plasma membrane localization in the S3 segment of the proximal tubule of the outer medulla"; "hybridization was greatest when RNA extracted from the outer medulla was used" | strong |
| L48 | Slco1a1 (rat oatp) | rat | both, gonadectomy | renal oatp strongly male-dominant and testosterone-driven (6 transcripts male, 1 female) | male > female | Northern + hormones | 8779895 | yes-abs | "In female rats, renal oatp mRNA expression was markedly less"; "Castration of male rats produced a dramatic drop in the steady-state level of all six renal oatp transcripts" | strong for sex |
| L49 | Slco1a1 (mouse Oatp1) | mouse | androgen-manipulated | strictly androgen-dependent renal expression; substrates steroid conjugates and PGE2 | male-restricted | Northern + oocyte transport | 11267661 | yes-abs | "we describe the strict androgen-dependent expression of mouse Oatp1 mRNA in kidney" | medium |
| L50 | Slco1a1, Slco1a6 | mouse | both | Oatp1a1 and 1a6 kidney-high; Oatp1a1 male > female | male > female | bDNA | 15843488 | yes-abs | "In kidney, expression of Oatp1a1, 3a1, and 4c1 was higher in males than in females" | medium |
| L51 | Slco1a1 | mouse (C57BL/6) | M, F, nephron AR-KO | **Slco1a1 is an S2 male-specific gene**, abolished by nephron AR removal; also has male-biased S3 accessible regions | S2 male-specific | RNAscope + multiome | 37673062 | yes-full | "Slco1a1 mRNA was highly male-specific among PT-S2, and entirely absent upon AR removal from male nephrons"; "there were 5 M-biased DARs ... in the genomic region of Slco1a1 in PT-S3" | strong |
| L52 | SLCO1A2 (OATP1A2) | human | not stated | the single human SLCO1A ortholog protein is in the **distal nephron**, not PT | PT-absent (distal) | IHC | 15632119 | yes-abs | "we demonstrate OATP1A2 protein expression in human brain capillary and renal distal nephron using immunohistochemistry" | medium (one study) |
| L53 | Slc22a19 (Oat5) | rat + mouse | both, gonadectomy | brush-border protein strongly in **S3** (outer stripe, medullary rays), weak in cortical S1/S2; female-dominant | S3 > S1/S2; female > male | IHC + Western + RT-PCR | 20865662 | yes-abs | "Oat5 protein was localized to the brush-border of S1/S2 in the cortex (CO) (weakly) and in S3 of the OS and medullary rays (strongly) with the F-dominant intensity"; "In adult mice, the sex-dependent pattern ... was comparable to that in rats" | strong |
| L54 | Slc22a19 (mOAT5) | mouse | not stated | apical membrane of PT of the **outer medulla** | S3 | immunolocalisation + transport | 16150593 | yes-abs | "mainly localized in the apical membrane of the proximal tubules of the outer medulla" | medium |
| L55 | Slc22a19 (Oat5) | mouse + rat | both | **negative/conflicting**: the identification paper found no sex difference | flat (sex) | Northern | 15068970 | yes-abs | "No sexual dimorphism in renal or hepatic expression of Oat5 was observed" | medium (conflicts with L53) |
| L56 | Slc22a19 (Oat5) | mouse (C57BL/6) | both, gonadectomy | female-predominant renal Oat5, estrogen-driven | female > male | bDNA + hormones | 19679677 | yes-full | "Oat5 mRNA expression was 100% higher in kidneys of female than male mice" | strong for sex |
| L57 | Slc22a22 (Oat-PG) | mouse (ICR) | not stated | kidney-specific; protein on the **basolateral** membrane of proximal convoluted tubules in cortex and outer stripe; colocalises with 15-PGDH | S1/S2 (PCT) + outer stripe | Northern + IHC + double IF | 20448048 | yes-full | "OAT-PG was detected at the basolateral membrane of proximal tubules ... in the cortex and the outer stripe of outer medulla"; "strong signals in the basolateral membrane of proximal convoluted tubules" | strong |
| L58 | Slc22a22 (rOAT-PG) | rat | both, castration | renal cortical expression much higher in males; castration lowers it, testosterone restores | male > female | Western/qPCR | 22031854 | yes-abs | "the expression level of rOAT-PG in the renal cortex was much higher in male rats than in female rats"; "Castration of male rat led to a remarkable reduction" | strong for sex |
| L59 | Slc22a19, Slc22a26-a30 | mouse/rat vs human | - | this subclade is a **rodent-specific expansion**; the human side is SLC22A9/A10/A24/A25 | no 1:1 human partner | phylogenetics | 32150922 | yes-full | "This subgroup also includes a large rodent-specific expansion, consisting of Slc22a19 and Slc22a26-30" | strong |
| L60 | SLC22A24 | human | donor sex not reported | **selectively expressed in human PT-S3**, absent elsewhere in the nephron; no rodent ortholog | S3 only | human snRNA re-analysis + transport | 31553721 | yes-full | "Overall, the plots indicate that SLC22A24 is selectively expressed in the proximal tubule segment 3 (PT(S3))"; "Phylogenetic analysis reveals the presence of SLC22A24 orthologs in apes, lemur, and certain other mammalian species, but none in rodents." | medium |
| L61 | SLC22A9 (OAT7), SLC22A10 | human | n/a | the human members of the rodent renal expansion are **liver** transporters | liver, not PT | siRNA/ChIP + 75 liver biopsies | 20829431 | yes-abs | "Three OATs are present in human liver: OAT2, OAT5, and OAT7" | medium |
| L62 | SLC22A10 | human | n/a | human **unitary pseudogene** (fixed P220); great-ape orthologs transport estradiol-glucuronide | non-functional in all humans | comparative functional | 38782905 | yes-abs | "Human SLC22A10 is a unitary pseudogene due to a fixed missense mutation, P220, while in great apes, its orthologs transport sex steroid conjugates" | strong |
| L63 | SLC22A11 (OAT4) | human | not stated | apical PT protein; **no rodent ortholog**, mouse Oat5 is a functional not genetic analog | PT apical; segment unresolved | IHC + transport (15037815); CHO kinetics (38627096) | 15037815 | yes-abs | "The immunoreactivities of OAT4 were observed in the apical membrane of proximal tubules"; (38627096) "Although there is no rodent ortholog of hOAT4, rodents express Oat5 (Slc22a19)" | medium |
| L64 | Slc22a20 (Oat6) | mouse | not stated | characterised in olfactory mucosa; **no kidney segment data found** | n/a | transport assays | 16478971 | yes-abs | "a novel Slc22 gene family member expressed in murine olfactory mucosa" | n/a |
| L65 | ABCB1 (P-gp) | human | not stated | P-gp concentrated on the **apical** surface of PT epithelium; no S1/S2/S3 | PT apical | IHC (MRK16) | 2444983 | yes-abs | "In kidney, P-glycoprotein was found concentrated on the apical surface of epithelial cells of the proximal tubules" | strong for PT |
| L66 | ABCB1 | human | not stated | MDR1 is **not PT-exclusive**: also thick limb, CD, mesangium | PT + TAL + CD | ISH + IHC | 9375821 | yes-abs | "The presence of MDR mRNA was confirmed in proximal tubules and mesangium, and demonstrated as well in thick limb of Henle's loops and in collecting ducts" | medium |
| L67 | Abcb1a / Abcb1b | mouse (C57BL/6) | both, gonadectomy/hypox | both renal Mdr1a and Mdr1b are **female-predominant** via androgen repression (opposite sign to Slco1a1/Slc22a22) | female > male | bDNA + endocrine surgery | 18854377 | yes-full | "Mdr1a mRNA was 83% higher in female than in male kidney"; "Both Mdr1a and Mdr1b in kidney show female-predominant expression patterns due to repression by androgens" | strong for sex |
| L68 | Abcb1a/1b (function) | mouse (WT vs KO) | not stated | functional P-gp efflux demonstrated in the isolated perfused **S2** segment | S2 (functional) | isolated perfused tubule | 11134265 | yes-abs | "studies were conducted in the isolated and perfused proximal tubule S2 segment" | medium |
| L69 | Mt1/Mt2 protein (MT-I/II) | rat | male | baseline/repleted renal MT localises primarily to **proximal convoluted tubules of the cortex** | S1/S2 (PCT) > S3 | MT immunostaining | 19486827 | yes-abs | "Dietary zinc repletion for 24h restored renal zinc and MT concentrations, the latter primarily in the proximal convoluted tubules of the cortex" | medium |
| L70 | MT1 cluster / MT2A | human | not stated | 10 functional + 7 non-functional human MT isoforms vs 2 mouse genes; MT immunoreactivity **confined exclusively to proximal tubules**; MT-1E/1F/1X/2A consistently expressed; MT-1A in only 2/6 kidneys; MT-1B/1G/1H/3/4 never detected | PT-restricted; per-gene presence donor-variable | IHC + isoform RT-PCR | 8619255 | yes-abs | "The organization of the metallothionein (MT) gene family has been demonstrated to be much more complex in humans than in the mouse"; "MT immunoreactivity to be confined exclusively to the proximal tubules"; "Messenger RNA for the MT-1A gene was detected in 2 of 6 renal samples" | strong for PT, no subsegment |
| L71 | MT-I/II | human | not stated | MT-I/II protein in proximal tubules, reduced in CKD | PT; segment unresolved | IHC biopsies | 34437439 | yes-abs | "MT-I/II is located on the proximal tubules and is notably reduced in patients with CKD" | medium |
| L72 | Mt (isoform not resolved) | rat | not stated | after ischemia-reperfusion MT predominates in the **outer medulla** - an injury pattern distinct from the baseline PCT pattern | S3/outer medulla under AKI | Northern + ISH + IHC | 12760485 | yes-abs | "MT was expressed exclusively in the renal tubular epithelial cells ... and that MT predominated in the outer medulla in the IARF rat kidney" | medium (AKI arm only) |
| L73 | Serpina1 | *Mus caroli* | male, T-treated | renal AAT confined to **a subset of PT cells**; renal expression is a species peculiarity | PT subset; species-specific | ISH + Southern | 3497926 | yes-abs | "The extrahepatic expression of alpha 1-antitrypsin is limited to the kidney, specifically within a subset of the proximal tubule cells" | medium (not *M. musculus*) |
| L74 | Serpina1 | mouse (CD-1) + human AKI | not stated | AAT is "renal silent" at baseline and induced 5-10x in cortex/PT by AKI | baseline low; injury-induced | qPCR/protein + human urine | 24848503 | yes-abs | "Does acute kidney injury (AKI) up-regulate the normally renal silent AAT gene?"; "Rapidly following injury, 5-10 fold renal cortical and isolated proximal tubule AAT mRNA and protein increases occurred" | strong for baseline-low |
| L75 | Serpina3c/k (cluster) | rat + human | not stated | SerpinA3 relocates cytoplasm -> apical tubular membrane under injury; **no healthy segment map found** | injury-state marker | transfection, IHC, urine (35147994); rat AKI-to-CKD IHC (31316093) | 35147994 | yes-abs | "serpinA3c/k relocates from the cytoplasm to the apical tubular membrane (ATM) in chronic kidney disease"; (31316093, yes-full) "serpinA3 is expressed mostly in the cytoplasm of tubular epithelium" in controls | weak for baseline position |
| L76 | Mup family | mouse | male > female | MUPs are synthesised in **liver**, not kidney; renal Mup signal is filtered/reabsorbed protein | kidney not a synthesis site | 14C incorporation, 9 tissues | 5827345 | yes-abs | "Mouse liver (but not eight other tissues examined) incorporated amino acids labeled with carbon-14 into the complex in vitro" | weak method, on-point |
| L77 | Mup / alpha2u-globulin uptake | rat (F344), mouse | male | the liver-made lipocalin accumulates in the **P2 (approximately S2)** segment; humans have no equivalent protein | S2 accumulation; human-absent | histopathology + IHC (review of primary data) | 7945892 | yes-abs | "this protein has been found to accumulate excessively in the P2 segment cells of renal proximal tubules, resulting in hyaline droplet formation" | medium (review) |
| L78 | Serpina6 (CBG) | mouse | not stated | renal CBG is a **developmental** PCT gene of juxtamedullary nephrons, undetectable in adult | adult PT-absent | ISH + IHC, postnatal series | 8425503 | yes-abs | "CBG mRNA in the kidney increased to a maximum during the third week of life, but was undetectable 3 weeks later" | medium |
| L79 | APOL1 | human | not stated | **contested PT localisation**: Madhavan reports proximal tubular epithelia; Blessing finds podocytes/endothelium only and blames antibody non-specificity | unresolved; not a safe PT gene | IHC (21997392); protein + in situ mRNA (34138902) | 21997392 | yes-abs | 21997392: "APOL1 localized to proximal tubular epithelia in normal kidney, FSGS, and HIVAN"; 34138902: "APOL1 protein and mRNA was detected in podocytes and endothelial cells, but not in tubular epithelia" | weak (conflict) |
| L80 | APOL1 | primates | n/a | primate-specific protein; mouse has only the many-to-many Apol7-11 family | human-only column | functional/evolutionary | 41894328 | yes-abs | "Some primates are immune to infection by most African trypanosome parasites due to apolipoprotein L-1 (APOL1), a primate-specific ion channel-forming protein" | strong |
| L81 | Slc7a13 / rBAT (Slc3a1), b(0,+)AT | mouse | not stated | rBAT highly expressed in **S3**, b(0,+)AT highest in **S1**; AGT1/Slc7a13 apical in S3 | S3 vs S1 split (1:1 genes) | biochemistry + localisation | 26739563 | yes-abs | "rBAT is highly expressed in the S3 segment, the late proximal tubules, whereas b(0,+)AT expression is highest in the S1 segment"; "AGT1 is localized in the apical membrane of the S3 segment" | strong (positive control) |
| L82 | SGLT2 / SGLT1 (human) | human | not stated | **hSGLT2 in S1/S2, SGLT1 in S3** by IHC - the 1:1 positive control for a paralog split | S1/S2 vs S3 | IHC | 25304002 | yes-abs | "hSGLT2 and SGLT1 were localized to the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively" | strong |
| L83 | OCT1 / OCT2 (rat) | rat | not stated | protein splits S1/S2 vs S2/S3 but **both mRNAs are present in all three segments** | mRNA flat, protein split | ISH + IHC | 10997918 | yes-abs | "rOCT1 protein was mainly observed in S1 and S2 segments ... rOCT2 protein was mainly expressed in the S2 and S3 segments"; "significant amounts of both rOCT1 and rOCT2 mRNA were detected in S1, S2, and S3" | strong (mRNA vs protein caveat) |
| L84 | human nephron OAT1/OAT2/OAT3 + NaDC3 | human | M + F | all three OATs more intense in S1/S2 than in S3 of medullary rays; S3 of the outer stripe unstained; **no sex difference detectable** | S1/S2 > S3 | validated antibodies | 27053689 | yes-abs | "all three OATs were stained more intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment in the outer stripe remained unstained"; "sex differences in their expression could not be detected" | strong (human) |
| L85 | human atlases: PT-segment resolution | human | mixed donors | Lake 2023 merged PT-S1/PT-S2 as unresolvable; Muto 2021 found S1/S2-vs-S3 unclear in snRNA | at most S1/S2 vs S3 in human references | snRNA/snATAC atlases | 37468583 | yes-full | 37468583: "Cell types that could not be accurately resolved (PT-S1/PT-S2) were kept merged."; 33850129: "The delineation between S1/S2 and S3 was less clear in the snRNA-seq dataset" | strong (methodological) |
| L86 | human atlases: these families | human | - | full-text search of Lake 2023, Muto 2021, Hansen 2022 found **no** mentions of MT1*, UGT*, GSTA1, CYP4A11, AKR1C*, APOL1, SLC22A11 or SLC22A24 | no atlas anchor | text search of open full text | 37468583, 33850129, 35675394 | yes-full (0 hits) | absence verified by search over open full text | strong negative |

**Explicit "no segment-resolved data found"** (useful as pre-registered gaps): human CYP2D6 and mouse Cyp2d in kidney;
human CYP2J2; human AKR1C1-C4 in kidney; human CES1; human UGT2B15/UGT2B17 by segment; human SULT1A/SULT2A by segment;
Aox anywhere along the nephron by ISH/IHC; Nat8/Nat8f in any species; Serpina3 cluster at healthy baseline; Fxyd2;
Slc22a13; Mt isoform-resolved position along mouse S1->S3.

## 1B. Deposited microdissection reference (mouse, male and female PT subsegments)

The strongest segment-resolved numbers for these families are in the **deposited data** of two Knepper-lab papers, not
in their article text. Values below are mean TPM read from those files (`yes-data`), and must be cited as deposited data:

- `GSE212213_TPM_Male_Female.txt.gz` (series "Multi-omics Analyses Reveal Sex Differences in Mouse Renal Proximal
  Tubules"; paper PMID **36758122**, J Am Soc Nephrol 2023) - columns `PTS1 PTS2 PTS3` (male) and `fPTS1 fPTS2 fPTS3`.
- `GSE150338_Mouse_TPM_Replicates.txt.gz` (paper PMID **33769951**, J Am Soc Nephrol 2021) - 14 segments, replicate level.

Sex check performed here: in GSE212213 the `PTS*` columns have Ddx3y 10.7-17.6 and Xist 0.0-0.1 TPM, while `fPTS*` have
Ddx3y 0 and Xist 11.7-27.9 - so the first three columns are male, matching this project's cohort. In the GSE150338
replicate file every PT replicate has Ddx3y 9.3-21.8 and Xist <= 0.28. (One caveat: a web page associated with
PMID 33769951 states that only female proximal tubules were included; the deposited PT data are male by these markers.
Flagged, not resolved.)

Article-text support for the summary claim is `yes-full` from the PMC version of 36758122: "A total of 1547 sex-specific
transcripts were detected ... 79 in PTS1, 1177 in PTS2, and 766 in PTS3."; "There were 11 transcripts that were
consistently expressed at higher levels in males in all three proximal tubule segments, including organic ion
transporters (Slc22a22, Slc22a30, Slc22a28) and cytochrome P450 family members (Cyp4b1, Cyp2j13)"; "Similarly, 11
transcripts were enriched in all female proximal tubule segments, including Prlr, Id1, Ugt1a2, Gsta4, Rdh16f2, Tmem37,
Mt1, Sgk1, Gm22133, Trf, and Tns1."; mice were "Six to 8-week-old male and female C57BL/6 mice from Taconic".

Mean TPM (M = male PTS1/2/3, F = female fPTS1/2/3, DTL1 = descending thin limb mean from GSE150338, included because
mouse late S3 grades into thin-limb identity). "shape" = which segments are within 2-fold of the male maximum.

| family | gene | M S1 | M S2 | M S3 | F S1 | F S2 | F S3 | DTL1 | male shape |
|---|---|---|---|---|---|---|---|---|---|
| Kap (mouse-only) | Kap | 304 | 36933 | 75375 | 56 | 64255 | 95736 | 2131 | peak S3 |
| Cyp4a 8:2 | Cyp4a10 | 28 | 1047 | 419 | 23 | 674 | 287 | 10 | peak S2 |
| Cyp4a 8:2 | Cyp4a12a | 9 | 26 | 55 | 1 | 0 | 0 | 18 | peak S3 |
| Cyp4a 8:2 | Cyp4a12b | 5 | 10 | 0 | 1 | 1 | 0 | 0 | peak S2 |
| Cyp4a 8:2 | Cyp4a14 | 1 | 123 | 3 | 1 | 315 | 2 | 2 | peak S2 |
| Cyp4a 8:2 | Cyp4a31 | 137 | 1001 | 556 | 136 | 621 | 569 | 17 | peak S2+S3 |
| Cyp4a 8:2 | Cyp4a29/30b/32 | 0 | 0-1 | 0-1 | 0-1 | 0-1 | 0 | 0 | PT-low |
| Cyp2d 8:2 | Cyp2d9 | 2 | 385 | 213 | 0 | 81 | 93 | 9 | peak S2+S3 |
| Cyp2d 8:2 | Cyp2d12 | 62 | 15 | 0 | 30 | 56 | 1 | 2 | peak S1 |
| Cyp2d 8:2 | Cyp2d22 | 25 | 3 | 2 | 26 | 17 | 5 | 51 | peak S1 |
| Cyp2d 8:2 | Cyp2d26 | 722 | 13 | 3 | 812 | 145 | 60 | 1 | peak S1 |
| Cyp2d 8:2 | Cyp2d10/34/40 | 0 | 0 | 0 | 0 | 0 | 0 | 0-1 | PT-low |
| Cyp2j 8:1 | Cyp2j5 | 816 | 1134 | 836 | 561 | 608 | 421 | 42 | flat S1-S3 |
| Cyp2j 8:1 | Cyp2j11 | 24 | 279 | 241 | 3 | 277 | 252 | 17 | peak S2+S3 |
| Cyp2j 8:1 | Cyp2j13 | 64 | 192 | 19 | 1 | 0 | 1 | 1 | peak S2 |
| Cyp2j 8:1 | Cyp2j6 | 6 | 14 | 16 | 7 | 5 | 4 | 12 | peak S2+S3 (low) |
| Cyp4f 5:4 | Cyp4f13 | 28 | 77 | 50 | 28 | 46 | 45 | 44 | peak S2+S3 |
| Cyp4f 5:4 | Cyp4f14/15/16 | 0-4 | 1-3 | 1 | 0-5 | 0-9 | 0-16 | 0-18 | PT-low |
| Ugt2b 7:7 | Ugt2b34 | 90 | 0 | 0 | 77 | 5 | 0 | 0 | peak S1 |
| Ugt2b 7:7 | Ugt2b37 | 16 | 407 | 1 | 2 | 2 | 0 | 1 | peak S2 |
| Ugt2b 7:7 | Ugt2b38 | 1240 | 1294 | 2 | 616 | 13 | 0 | 1 | peak S1+S2 |
| Ugt2b 7:7 | Ugt2b1/5/35/36 | 0-1 | 0-4 | 0 | 0-1 | 0 | 0 | 0-2 | PT-low |
| Ugt1a | Ugt1a2 | 22 | 0 | 0 | 340 | 954 | 27 | 1 | female S2 |
| Ugt1a | Ugt1a7c | 239 | 370 | 33 | 144 | 73 | 5 | 7 | peak S1+S2 |
| Ugt1a | Ugt1a1/1a6a/1a9 | 0-6 | 0-3 | 0-1 | 0-41 | 0-8 | 0-1 | 0-12 | PT-low |
| Akr1c 9:4 | Akr1c14 | 2 | 370 | 590 | 0 | 165 | 245 | 58 | peak S2+S3 |
| Akr1c 9:4 | Akr1c18 | 0 | 200 | 956 | 2 | 1192 | 4894 | 173 | peak S3 |
| Akr1c 9:4 | Akr1c21 | 1733 | 1564 | 1570 | 995 | 963 | 42 | 110 | flat S1-S3 |
| Akr1c 9:4 | Akr1c13 | 7 | 14 | 29 | 10 | 6 | 1 | 25 | peak S3 (low) |
| Akr1c 9:4 | Akr1c6/12/19/20 | 0-5 | 0-4 | 0-3 | 0-8 | 0-35 | 0-18 | 0-21 | PT-low |
| Ces1 8:2 | Ces1d | 497 | 327 | 53 | 307 | 146 | 17 | 4 | peak S1+S2 |
| Ces1 8:2 | Ces1e | 113 | 27 | 0 | 73 | 11 | 0 | 8 | peak S1 |
| Ces1 8:2 | Ces1f | 138 | 1218 | 710 | 8 | 198 | 423 | 27 | peak S2+S3 |
| Ces2 7:1 | Ces2c | 42 | 460 | 77 | 20 | 135 | 106 | 0 | peak S2 |
| Ces2 7:1 | Ces2b | 1 | 38 | 9 | 0 | 51 | 49 | 0 | peak S2 |
| Ces2 7:1 | Ces2g | 15 | 1 | 0 | 20 | 1 | 0 | 1 | peak S1 |
| Nat8 9:2 | Nat8 | 36 | 560 | 361 | 7 | 200 | 3 | 2 | peak S2+S3 |
| Nat8 9:2 | Nat8f1 | 968 | 1092 | 172 | 825 | 1519 | 378 | 38 | peak S1+S2 |
| Nat8 9:2 | Nat8f5 | 12 | 80 | 8 | 2 | 34 | 0 | 0 | peak S2 |
| Nat8 9:2 | Nat8f6 | 5 | 54 | 78 | 1 | 46 | 0 | 1 | peak S2+S3 |
| Nat8 9:2 | Nat8f2/f3/f4/f7 | 2-15 | 1-28 | 0-3 | 4-9 | 1-6 | 0-3 | 1-7 | PT-low |
| Slco1a 6:1 | Slco1a1 | 3 | 359 | 40 | 0 | 0 | 0 | 6 | peak S2 |
| Slco1a 6:1 | Slco1a6 | 67 | 189 | 54 | 13 | 304 | 422 | 3 | peak S2 (sex sign flips) |
| Slco1a 6:1 | Slco1a4/1a5 | 0-1 | 0-3 | 0 | 0 | 0-2 | 0-1 | 2-18 | PT-low |
| Slc22 Oat-like 6:3 | Slc22a19 | 1 | 534 | 524 | 1 | 670 | 957 | 40 | peak S2+S3 |
| Slc22 Oat-like 6:3 | Slc22a26 | 0 | 44 | 31 | 0 | 39 | 33 | 1 | peak S2+S3 |
| Slc22 Oat-like 6:3 | Slc22a28 | 131 | 108 | 166 | 16 | 29 | 13 | 6 | flat S1-S3 |
| Slc22 Oat-like 6:3 | Slc22a29 | 39 | 1 | 1 | 48 | 38 | 35 | 1 | peak S1 |
| Slc22 Oat-like 6:3 | Slc22a30 | 38 | 145 | 269 | 5 | 23 | 22 | 10 | peak S2+S3 |
| Slc22a22 (mouse-only) | Slc22a22 | 60 | 116 | 60 | 6 | 31 | 1 | 0 | flat S1-S3 |
| Abcb1 2:1 | Abcb1a, Abcb1b | 0-2 | 2 | 0-1 | 1-3 | 3-13 | 0-2 | 1-12 | PT-low (both) |
| Gsta 5:4 | Gsta2 | 1006 | 22 | 8 | 1300 | 588 | 132 | 4 | peak S1 |
| Gsta 5:4 | Gsta3 | 69 | 47 | 19 | 87 | 54 | 2 | 21 | peak S1+S2 |
| Gsta 5:4 | Gsta1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | PT-absent |
| Gsta4 (unclassed) | Gsta4 | 75 | 1 | 1 | 263 | 140 | 81 | 33 | peak S1 |
| Gstm 6:4 | Gstm1 | 233 | 2064 | 5794 | 239 | 1994 | 3418 | 7257 | peak S3 |
| Gstm 6:4 | Gstm5 | 363 | 252 | 432 | 282 | 494 | 622 | 156 | flat S1-S3 |
| Gstm 6:4 | Gstm7 | 7 | 13 | 20 | 7 | 18 | 19 | 10 | peak S2+S3 (low) |
| Gstm 6:4 | Gstm2/3/6 | 0-6 | 0-2 | 0-1 | 0-7 | 1-7 | 1-3 | 1-1668 | PT-low (Gstm2 is DTL) |
| Sult | Sult1d1 | 342 | 181 | 293 | 189 | 333 | 318 | 87 | flat S1-S3 |
| Sult | Sult1a1 | 1 | 1 | 4 | 2 | 11 | 10 | 93 | PT-absent (DTL/DCT gene) |
| Sult2a 8:1 | Sult2a1, Sult2a8 | 0 | 0 | 0 | 0 | 0-1 | 0 | 0-2 | PT-absent |
| Mt 2:10 | Mt1 | 248 | 3 | 1 | 596 | 189 | 146 | 42 | peak S1 |
| Mt 2:10 | Mt2 | 10 | 0 | 0 | 48 | 4 | 3 | 11 | peak S1 (low) |
| Serpina3 10:1 | Serpina3a/c/k/m/n | 0 | 0 | 0 | 0 | 0 | 0 | 0-1 | PT-absent at baseline |
| Serpina1 5:1 | Serpina1a | 0 | 15 | 70 | 0 | 5 | 2 | 1 | peak S3 (low) |
| Mup | Mup1, Mup3, Mup20 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | PT-absent |
| Aox | Aox1, Aox3, Aox4 | 0-1 | 0 | 0 | 0 | 0-1 | 0 | 4-75 | PT-absent (Aox3 is DTL1) |
| Fmo 1:1 | Fmo1 | 946 | 524 | 250 | 1062 | 554 | 232 | 21 | peak S1+S2 |
| Fmo 1:1 | Fmo2 | 169 | 136 | 57 | 150 | 59 | 2 | 245 | peak S1+S2 |
| Fmo 1:1 | Fmo5 | 10 | 250 | 179 | 2 | 10 | 5 | 65 | peak S2+S3, male |
| Akr1b 3:2 | Akr1b7, Akr1b8 | 0-2 | 0 | 0 | 0-2 | 0 | 0 | 0-13 | PT-absent |
| Rdh16 6:1 | Rdh16 | 52 | 89 | 32 | 55 | 133 | 102 | 5 | peak S1+S2 |
| Rdh16 6:1 | Rdh16f2 | 12 | 132 | 8 | 46 | 888 | 488 | 0 | peak S2 |
| Slc7a12 (mouse-only) | Slc7a12 | 1 | 11 | 3 | 2 | 1613 | 3048 | 4 | female S2+S3 |
| Acsm2 1:2 | Acsm2 | 4478 | 8281 | 3310 | 3156 | 561 | 9 | 50 | peak S1+S2 |
| Apol 11:4 | Apol7a/9a/9b/10a | 0-1 | 0-2 | 0-1 | 0 | 0-2 | 0-1 | 1-14 | PT-absent |
| H2 class I 15:7 | H2-D1 | 168 | 250 | 282 | 184 | 246 | 159 | 1678 | flat S1-S3 |
| H2 class I 15:7 | H2-T23 | 72 | 202 | 331 | 69 | 253 | 460 | 277 | peak S2+S3 |
| H2 class I 15:7 | H2-K1 | 137 | 55 | 87 | 163 | 165 | 126 | 944 | S1 and S3 |
| Saa 3:3 | Saa1, Saa2, Saa3 | 0 | 0 | 0 | 0 | 0 | 0 | 1-2 | PT-absent at baseline |
| 1:1 landmarks | Slc5a2 | 2722 | 1 | 0 | 2808 | 1 | 0 | 13 | peak S1 |
| 1:1 landmarks | Slc5a12 | 332 | 0 | 0 | 376 | 0 | 0 | 1 | peak S1 |
| 1:1 landmarks | Slc34a1 | 9487 | 5646 | 44 | 10241 | 2227 | 17 | 2 | peak S1+S2 |
| 1:1 landmarks | Slc7a13 | 6 | 2042 | 1502 | 0 | 1 | 3 | 43 | peak S2+S3, male |
| 1:1 landmarks | Slc22a7 | 2 | 2463 | 1460 | 1 | 877 | 1155 | 21 | peak S2+S3 |
| 1:1 landmarks | Cyp7b1 | 0 | 397 | 1068 | 0 | 0 | 3 | 68 | peak S3, male |
| 1:1 landmarks | Cyp2e1 | 63 | 991 | 0 | 3 | 94 | 0 | 7 | peak S2, male |

Pre-registrable patterns visible in the male columns (to be treated as the expected profile, not as our result):
S1-restricted = Cyp2d26, Cyp2d12, Ces1e, Ces2g, Ugt2b34, Gsta2, Gsta4, Mt1, Slc22a29 (and Fmo1/Fmo2 declining);
S2-peaked = Cyp4a10, Cyp4a31, Cyp4a14, Cyp2d9, Cyp2j13, Cyp2j11, Cyp4f13, Ces1f, Ces2c, Ugt2b37, Ugt1a7c, Nat8f5,
Slco1a1, Slco1a6, Rdh16f2; S3-rising = Kap, Gstm1, Akr1c18, Akr1c14, Cyp4a12a, Slc22a30, Nat8f6, H2-T23;
flat = Cyp2j5, Akr1c21, Gstm5, Sult1d1, Slc22a28, Slc22a22, H2-D1; PT-absent = all Aox, all Sult2a, Sult1a1, Akr1b7/b8,
Gsta1, all Serpina3, all Saa, all Mup, Apol, Abcb1a/b, Ugt1a9.

Two literature-vs-data conflicts to pre-register rather than smooth over:
1. Cyp4a14 is reported female-**only** in S2 by RNAscope (31689386) but reads 123 TPM in male PTS2 vs 315 in female
   (about 2.6-fold), and Akr1c18 is called "female S3 enriched" yet reads 956 TPM in male S3. In an all-male cohort,
   treat ISH "female-only" statements as relative, not absolute.
2. Slco1a6 flips sex sign between segments (male-biased in PTS1, female-biased in PTS3). A constant male/female offset
   model cannot represent that gene.

---

# 2. Sex dimorphism along the mouse PT

## 2.1 Key papers

| paper | PMID | verified | species / strain / sex / method | what it reports (quote) |
|---|---|---|---|---|
| Ransick 2019 Dev Cell | 31689386 | yes-full | mouse C57BL/6J 8-10 wk, 2M + 2F, zone-dissected scRNA + RNAscope | "all proximal tubule clusters separated by sex"; "The highest disparity in gene activity between the sexes was in the S3 region."; "S2 and S3 sex-specific expression domains form a tight, non-overlapping boundary at the S2-S3 junction"; "Slc7a12 and Cyp7b1 transcripts localized exclusively to female (cluster 9) and male (cluster 9m) subsets, respectively, of Slc22a7+ S3 cells" |
| Xiong 2023 Dev Cell | 37673062 | yes-full | mouse C57BL/6; castration, ovariectomy, testosterone, nephron-specific Six2-Ar-KO / Six2-Esr1-KO, systemic Sox2-Ar-KO; snRNA+snATAC multiome; RNAscope | "sexually dimorphic gene activity predominantly maps to PT segments 2 and 3 (PT-S2 and PT-S3)"; "PT segment 1 (PT-S1) cluster comprised S1 cells from all genotypes and sexes, indicative of a low level of AR-dependent variability"; "Of note, PT-S3 showed the highest number of sex-biased genes"; "After castration, 85.5% of the core M-biased genes (288 of 337) were down-regulated and 77.6% of core F-biased genes (166 of 214) were up-regulated"; "nephron removal of either Ar or Esr1 had little effect on gene expression in the female kidney"; "F-biased DARs showed a strong motif enrichment for Stat5a/b and Bcl6 which lie downstream of prolactin and growth hormone signaling" |
| Xiong 2023 bioRxiv (same study) | 37205355 | yes-full | as above | "M-biased genes such as Slco1a1, Cyp2j13, Cyp7b1, Slc7a13, Acsm3, and Tbx10, were significantly down-regulated in both nephron-specific and systemic removal of AR; F-biased genes Ugt1a2, Scl7a12, Prlr, and Cyp4a14 were most up-regulated" |
| Chen 2023 JASN | 36758122 | yes-full | mouse C57BL/6 (Taconic) 6-8 wk, both sexes; microdissected PTS1/2/3 RNA-seq + ATAC-seq + whole-kidney proteomics | "Sex-specific transcripts include: 79 in PTS1, 1177 in PTS2, and 766 in PTS3."; male in all three segments: "Slc22a22, Slc22a30, Slc22a28 ... Cyp4b1, Cyp2j13"; female in all three: "Prlr, Id1, Ugt1a2, Gsta4, Rdh16f2, Tmem37, Mt1, Sgk1, Gm22133, Trf, and Tns1"; "neither androgen nor estrogen response elements (ARE and ERE) were over-represented in the promoter regions" |
| Rodriguez-Montes 2023 Science | 37917687 | yes-full | human, mouse, rat, rabbit, opossum, chicken; 5 organs incl. kidney; developmental bulk RNA-seq + scRNA re-analysis | "in mouse, only 15 genes are sex-biased in the brain, whereas 1891 genes are Sex-biased in the kidney"; "Male-biased genes are expressed specifically in the proximal tubule cells, while female-biased genes are not cell-type specific"; "in the kidney, most mouse sex-biased genes have 1:1 orthologs in the rat (87%), in the liver, 28% of mouse sex-biased genes do not have a 1:1 ortholog"; "For most genes in the kidney, the sex differences are recent, but the genes themselves are old"; newly emerged sex-biased genes "derived from the expansion of three gene families, the cytochrome P450 family (20/77 genes), the major urinary protein family (16/77 genes) and the Slc22 transporter family (5/77)"; "The rodent-specific expansion of the Slc22 transporter family has been associated with the transport of conjugated sex hormones" |
| Pihlajamaa 2014 EMBO J | 24451200 | yes-abs | mouse prostate, **kidney**, epididymis; in vivo AR ChIP-seq | "in vivo AR cistromes and their respective androgen-dependent transcription programs are highly tissue specific"; "two novel collaborating factors for AR signaling in vivo--Hnf4a ... in mouse kidney" |
| Harris 2021 AJP Renal (KS-AR-KO) | 34605272 | yes-abs | mouse, kidney-specific AR deletion, both sexes | "In male but not female mice, KS-AR-KO increased ammonia excretion, which eliminated sex differences"; "KS-AR-KO decreased kidney size, cortical proximal tubule volume density, and cortical proximal tubule cell height in males" |
| Harris 2021 AJP Renal (review) | 33308019 | yes-abs | review, mouse | "the canonical testosterone target receptor, androgen receptor (AR), is present exclusively in the proximal tubule" |
| Cheng 2006 Mol Pharmacol | 16807376 | yes-abs | mouse: gonadectomy, hypophysectomy, lit/lit, GH patterns | "Androgens increased Oatp1a1 mRNA in liver and kidney, whereas male-pattern GH administration increased Oatp1a1 mRNA in livers but not in kidneys."; "In kidney, gender-divergent Oatp expression is exclusively caused by stimulation by androgens." |
| Cheng 2009 DMD | 19679677 | yes-full | mouse, both sexes, hormone manipulation, 14 renal transporters | "Gender-divergent expression of male-predominant (Urat1 and Oatp4c1) and female-predominant (Oat5) transporters in mouse kidneys is primarily due to stimulatory effects of androgens and estrogens, respectively" |
| Maher 2006 Biochem Pharmacol | 16529719 | yes-abs | mouse, gonadectomy/hypophysectomy, DHT/E2/GH patterns | "Mrp3 and Mrp4 are expressed at much higher levels in female than male kidney"; "Mrp4 mRNA is higher in females due to repression by both DHT and MP-GH secretion in males" |
| Buckley 2009 DMD | 19131521 | yes-abs | mouse, hypophysectomy + GH patterns, gonadectomy | "gender differences in mouse Ugt mRNA expression were influenced by male-pattern GH secretion in liver, whereas gender differences were regulated by the effects of androgens in kidney" |
| Rinn 2004 Dev Cell | 15177028 | yes-abs | mouse, both sexes; microarray of hypothalamus, kidney, liver, reproductive tissues | "Extensive differential gene expression was observed not only in the reproductive tissues, but also in the kidney and liver. The differentially expressed genes are involved in drug and steroid metabolism, osmotic regulation" |
| McEvoy 2022 Nat Commun | 36496458 | yes-full | **human**, 19 living donors, sex-balanced, pre-implantation biopsies; scRNA + primary PT cells | "Sex-based differences in gene expression within proximal tubular cells were observed, specifically, increased anti-oxidant metallothionein genes in females and aerobic metabolism-related genes in males."; "we identified genes with significant differential expression between males and females (n = 75 genes...)"; "Several sex-biased genes are consistent with previous reports of genes upregulated in murine male (NAT8, FKBP5, KDM5D, DGKG) and female (MGST3, SLC3A1, CYP4A11, RPS29) PT cells" |
| Oliva 2020 Science (GTEx v8) | 32913072 | yes-full | human, 838 donors, 44 tissues incl. kidney cortex | "A total of 37% of all genes exhibit sex-biased expression in at least one tissue."; "Only 30 genes (0.23%) ... exhibited consistent sex bias across all 44 tissue sources" |
| Idda 2021 Biomolecules | 34439872 | yes-full | human GTEx v8, 289 pharmacogenes, 41 tissues incl. renal cortex | "In the kidney only two SBDR genes with a single transcript were identified: the phospholipase A2 group IIA membrane enzyme (PLA2G2A) and the solute carrier family 2 member 9 (SLC2A9)" |
| Huang 2020 Gene | 32439376 | yes-abs | human kidney scRNA, both sexes | "the SLC22 family of anion transporters, including SLC22A6 and SLC22A8, had different expression profiles between male and female proximal tubule clusters but not sex-dependent abundance at the protein level" |
| Breljak 2016 AJP Renal | 27053689 | yes-abs | **human** kidneys, M and F, validated antibodies | "sex differences in their expression could not be detected" (OAT1/OAT2/OAT3, NaDC3) |
| Clotet-Freixas 2024 Sci Transl Med | 38446901 | yes-abs | human primary PTECs (M/F), mouse FCG, metabolic flux | "Primary human proximal tubular epithelial cells (PTECs) from healthy males displayed increased mitochondrial respiration, oxidative stress, apoptosis"; "The male human PTEC phenotype was enhanced by dihydrotestosterone and mediated by the transcription factor HNF4A and histone demethylase KDM6A" |
| Ljubojevic 2004 AJP Renal | 15010355 | yes-abs | rat, adult + prepubertal, gonadectomy + hormones; immunocytochemistry | "OAT1 was strong in proximal tubule S2 and weak in the S3 segments, whereas OAT3 was stained in proximal tubule S1 and S2 segments"; "Gender differences ... (M > F) were observed only in the cortical tubules"; "both a stimulatory effect of androgens ... and an inhibitory effect of estrogens" |
| Ljubojevic 2007 AJP Renal (OAT2) | 16885152 | yes-abs | rat + mouse, gonadectomy + replacement | "the antibody exclusively stained brush border (BB) of the proximal tubule S3 segment (S3) in the OS and medullary rays (F > M)"; "caused by strong androgen inhibition"; "In mice, the expression pattern largely resembled that in rats" |
| Breljak 2013 AJP Renal | 23389457 | yes-abs | mouse (Oat3-KO + WT), castration/testosterone | "mOat3 and mOat1 exhibit sex-dependent expression with an opposite pattern; mOat3 is female dominant due to androgen inhibition, while mOat1 is male dominant due to androgen stimulation" |
| Sabolic 2006 AJP Renal (SGLT1) | 16204409 | yes-abs | rat, gonadectomy + hormones | "strongly stained BBM of the S3 PT segments in the OS and medullary rays (F > M)"; "gender differences (F > M) ... caused by androgen inhibition" |
| Sabolic 2012 AJP Cell (SGLT2) | 22262063 | yes-full | rat and mouse, both sexes, castration + hormones | "rSGLT2 protein was localized to the brush-border membrane (BBM) of the renal proximal tubule S1 and S2 segments (S1 > S2) with female-dominant expression in adult rats"; "At variance to the findings in rats, the abundance of mSGLT2 protein in the mouse kidneys was male dominant, whereas the expression of mSglt2 mRNA was female dominant" |
| Buist 2004 DMD | 15155553 | yes-abs | mouse 129J and C57BL/6, both sexes | "Mouse Oat1 mRNA was primarily expressed in kidney of both strains, with male predominance."; "only 129J mice had higher Oat3 mRNA levels in female kidney than in male kidney"; "The most notable species differences ... were a lack of Oat2 female predominance in mouse kidney" |
| Veiras 2017 JASN | 28774999 | yes-abs | rat + C57BL/6 mice, both sexes; quantitative immunoblotting | "the proximal tubule of females had greater phosphorylation of Na+/H+ exchanger isoform 3 (NHE3) ... and less abundant expression of Na+/Pi cotransporter 2, claudin-2, and aquaporin 1"; "we observed similar abundance profiles in female versus male C57BL/6 mice" |
| Breljak 2026 Int J Mol Sci (Klotho) | 41596364 | yes-abs | rat vs mouse | "Female-dominant expression of renal rKl, affected by androgen's inhibitory effect"; "Renal mKl expression was male-dominant" |
| Landau 2026 bioRxiv (not kidney) | 41889965 | yes-abs | mouse vs human salivary gland, liver, pancreas | "These genes cluster in regions shaped by recent gene duplication, such as the kallikrein gene cluster, a mouse-specific expansion that accounts for ~16.4% of male-biased submandibular expression"; "68% and 73% of expression from genes encoding secreted proteins derives from lineage-specific genes that lack one-to-one human orthologs" |

Discarded for kidney claims after checking: **Yang 2006 Genome Res (16825664)** - tissues are liver, adipose, muscle,
brain only, no kidney; **Naqvi 2019 Science (31320509)** - its 12 tissues verifiably exclude kidney and it is restricted
to 12,939 1:1 orthologs, so it is structurally blind to expanded families. Both are usable for general framing only.

## 2.2 Answers

**(a) Which segments carry the sex differences.** S2 and S3; S1 is near-flat. Three independent designs agree:
microdissection counts 79 / 1177 / 766 sex-specific transcripts in PTS1 / PTS2 / PTS3 (36758122); the snRNA multiome
finds PT-S1 pooling across sexes and genotypes while "PT-S3 showed the highest number of sex-biased genes" (37673062);
zone-dissected scRNA reports the greatest disparity in S3 (31689386). The two resources disagree on whether S2 or S3 is
the peak, and no paper reconciles that. Importantly for a continuous coordinate, the sex boundary is reported as sharp:
"S2 and S3 sex-specific expression domains form a tight, non-overlapping boundary at the S2-S3 junction" (31689386) -
i.e. sex-restricted genes behave as step functions at a segment edge, not as smooth gradients.
Canonical male-biased: Slco1a1 (S2), Cyp2j13 (S2), Cyp2e1 (S2), Cyp7b1 (S3), Slc7a13 (S3), Acsm3 (S3->S2), Cyp4a12a,
Cyp4b1, Fmo5, Slc22a22/a28/a30 (all three segments), Kap (cortical component), Odc1, Gusb, Tbx10.
Canonical female-biased: Abcc3 (S3), Slc7a12 (S3), Prlr (S3->S2), Cyp4a14 (S2), Cyp4a10, Gsta4, Ugt1a2, Mt1, Akr1c18
(S3), Dio1, Ttr, Gm4450, Rdh16/Rdh16f2, Slc22a7 (S3), Slc22a19 (S3). Most sex-biased transcripts are autosomal
(36758122), and the bias appears at puberty (37673062).

**(b) Androgen vs estrogen/prolactin/GH.** Androgen acting directly through AR in PT cells is the primary mouse driver:
castration reverses 85.5% of core male-biased and 77.6% of core female-biased genes, nephron-specific Six2-Ar-KO
reproduces it, and Six2-Esr1-KO or ovariectomy do little (37673062). AR-dependent male genes: Slco1a1, Cyp2j13, Cyp7b1,
Slc7a13, Acsm3, Tbx10; androgen-repressed female genes: Ugt1a2, Slc7a12, Prlr, Cyp4a14 (37205355), plus Abcc3 and Gsta4
(37673062). The AR binding data are Pihlajamaa's in vivo mouse-kidney AR cistrome with Hnf4a as the collaborating
pioneer factor (24451200). The female program lacks AR motifs and instead shows Stat5a/b + Bcl6 (prolactin/GH-STAT5)
motifs - a hypothesis in that paper, not a perturbation (37673062). Dissenting measurement: Chen found no ARE/ERE
over-representation **in promoters** and proposed indirect control via Tead1, Nfia/b, Pou3f3 with Hnf1b/Hnf4a
(36758122); Xiong's AR-motif enrichment is in distal accessible regions. Promoter-only vs distal scanning is the likely
reconciliation but no paper states it, so it stays an open conflict.
Single-gene pharmacology agrees and separates kidney from liver: kidney Oatp sex bias is "exclusively caused by
stimulation by androgens" while liver uses male-pattern GH (16807376); kidney Ugt2b38 is testosterone-induced and
Ugt1a2 testosterone-repressed, again with the GH-pattern mechanism confined to liver (19131521); Urat1 and Oatp4c1 are
androgen-stimulated and Oat5 estrogen-stimulated (19679677); Mrp3/Mrp4 are female-biased, Mrp4 repressed by both DHT and
male-pattern GH (16529719); Cyp4a12a is DHT-inducible >4-fold (17112342). Not everything is AR: the **S3 component of
Kap** is AR-independent, estrogen-inducible and pituitary (TSH/thyroid)-dependent (2293028, 1338721, 7956935), and Gusb
needs pulsatile GH on top of androgen (10459853, 8363554). Phenotype-level: kidney-specific AR deletion eliminates the
sex difference in ammonia excretion and reduces cortical PT volume density and cell height in males only (34605272),
with AR reported as present exclusively in the PT (33308019).

**(c) Are androgen-dependent / sex-biased kidney genes enriched in rodent-expanded families?** **No published paper
tests family expansion as an enrichment hypothesis in kidney.** State this plainly; it is an open, testable question.
The closest is Rodriguez-Montes 2023 (37917687), the only paper that analyses gene age and duplication status of
sex-biased genes. In mouse **kidney** it finds the opposite of the naive expectation: "the sex differences are recent,
but the genes themselves are old", and 87% of mouse kidney sex-biased genes have 1:1 rat orthologs. Its expanded-family
annotation (P450 20/77, Mup 16/77, Slc22 5/77) refers to *recently emerged* sex-biased genes and is driven mainly by
liver; it does separately note that "The rodent-specific expansion of the Slc22 transporter family has been associated
with the transport of conjugated sex hormones". What is firmly established is only the weaker, function-level claim:
rodent kidney sex differences concentrate in drug/xenobiotic and organic-anion handling (15177028; 36758122; 37917687) -
exactly the functions carried by the expanded families (Cyp4a 8:2, Cyp2j 8:1, Slco1a 6:1, Slc22 Oat-like 6:3, Nat8f 9:2,
Ugt2b 7:7, Akr1c 9:4), and several top sex-biased PT genes have no human ortholog at all (Kap, Slc7a12, Slc22a22,
Gsta4-unclassed, Mup, Gm4450). The step from "sex-biased genes are in these functions" to "sex-biased genes are enriched
in expanded families against an expression-matched background" has not been published for kidney. Precedent from another
organ: in mouse submandibular gland, sex-biased expression does concentrate in a mouse-specific expansion (kallikrein
cluster, about 16.4% of male-biased expression) (41889965, unreviewed preprint).

**(d) Is human PT sexually dimorphic in these families?** Human PT **is** the dimorphic cell type, but the gene set is
different and much smaller. McEvoy's 19-donor sex-balanced scRNA finds PT separation with about 75 DE genes, female
metallothioneins (MT1F/MT1G/MT1H) and male aerobic metabolism (36496458). Direct conservation testing is negative for
autosomes: GTEx v8 kidney re-analysis gives 1 female-biased (KDM5C) and 73 male-biased genes (63 Y-linked), and the only
mouse-human conserved sex-biased genes are Y-linked UTY, DDX3Y, KDM5D; human snRNA (2F + 3M) showed PT co-clustering by
sex (37673062). Ortholog matching does recover 23 female- and 15 male-biased shared genes, including Cyp4a12a. Family
hints exist but with flipped direction: CYP4A11 is female-biased in human PT while mouse Cyp4a12a is male-specific, both
inside the same 8:2 family; NAT8 (family 9:2) is male-biased in human PT (36496458). Human SLC22A6/A8 differ by sex at
mRNA but not protein (32439376) and human OAT1/2/3 show no detectable sex difference by antibody (27053689); human kidney
cortex yields only two sex-biased pharmacogenes, PLA2G2A and SLC2A9 (34439872). Functional human PT dimorphism is
nonetheless real and androgen/HNF4A-dependent (38446901) - the same HNF4A node that licenses AR binding in mouse kidney
(24451200). Net: cell type conserved, mechanism node plausibly conserved, gene identity largely not. All human negatives
are small-n, so read them as "not detected at current power".

**(e) Rat vs mouse.** Rat supplies the segment-resolved protein map: OAT1 strong in S2, OAT3 in S1/S2, both M>F in cortex
only, androgen-stimulated and estrogen-inhibited (15010355); OAT2 and Oat5 are F>M on the **S3** brush border and
androgen-repressed (16885152, 20865662); SGLT1 is F>M in S3 (16204409); SGLT2 is S1>S2 and female-dominant in rat
(22262063). Mouse differs from rat in specific documented cases: mouse Oat3 is female-dominant via androgen inhibition
while Oat1 is male-dominant via androgen stimulation, "an opposite pattern" (23389457); SGLT2 protein is female-dominant
in rat but male-dominant in mouse, with mouse mRNA going the other way again (22262063); Klotho flips the same way
(41596364); mouse reportedly lacks rat's Oat2 female predominance (15155553, in conflict with the immunochemistry of
16885152). Strain matters inside mouse too (Oat3 female bias only in 129J; renal Cyp4a12a/20-HETE varies about 4-fold
across strains with C57BL/6 at the low end, 17112342). What is conserved rat->mouse is the Na+-handling PT profile
(28774999). Practical implication: rat is not a safe proxy for mouse at the level of individual transporter sex
direction, and mRNA is not a safe proxy for protein in either species.

---

# 3. Hypotheses from duplicate-gene evolution that a spatial expression axis can test

Notation to freeze with the protocol: `f_g(t)` = fitted mean expression of gene g along its own species' PT coordinate;
`p_g(t)` = shape = `f_g(t) / mean_t f_g(t)`; `D(a,b)` = shape distance (Euclidean on p, as in Assis 2015; 1 - Pearson as
secondary). `B` = baseline distribution of `D(mouse g, human g)` over PT-variable 1:1 orthologs matched on mean
expression, with cutoff `c` = median + semi-interquartile range of B (the Assis 2015 choice). `Sum_m(t)` = sum over all
mouse family members. Levels are compared only as **within-species percentiles**, never raw counts across panels
(the Qian 2010 design). `tau_axis` = tissue-specificity index tau computed over K equal-width t bins.

| id | hypothesis | key papers | PMID | verified | quoted key finding | testable prediction in PT data | what would count against it |
|---|---|---|---|---|---|---|---|
| D01 | Expression divergence between paralogs increases with duplicate age | Gu 2002 Trends Genet (yeast); Makova 2003 Genome Res (human tissues); Guschanski 2017 Genome Res (9 amniotes, 8 organs); Assis 2015 BMC Evol Biol (8 mammals, M+F tissues) | 12446139; 12840042; 28743766; 26173681 | yes-abs; yes-abs; yes-full; yes-abs | Gu: "expression divergence increases with evolutionary time". Makova: "a significant negative correlation between R and K(S)". Guschanski: "We observed an arch-shaped relationship between expression divergence ... and duplication age in all studied species". Assis: "young mammalian duplicates are often functionally conserved, and that expression divergence rapidly increases over evolutionary time" | Across mouse paralog pairs, `D(m_i, m_j)` of PT shapes rises with pair age (dS or node age); most rodent-specific pairs fall inside c; the most divergent within-species pairs are ancient paralogs (e.g. Slc5a1/Slc5a2) | Spearman(age, D) <= 0 after expression matching, or young pairs as divergent as ancient ones |
| D02 | A species-specific duplication lowers cross-species agreement of the orthologs' profiles | Huminiecki 2004 Genome Res (human vs mouse, 16 tissues) | 15466287 | yes-abs | "orthologs that have undergone recent duplication are less likely to have strongly correlated expression profiles than those that remain in a one-to-one relationship between human and mouse" | For each 1:n family, `D(best-matching mouse member, human gene)` exceeds B on average (permutation over expression-matched 1:1 sets) | The 1:n best-member D distribution is indistinguishable from B |
| D03 | Dosage sharing: duplicates are down-regulated so the summed family approximates the singleton | Qian 2010 Trends Genet; Lan 2016 Science (GTEx 46 human / 26 mouse tissues) | 20708291; 27199432 | yes-full; yes-full | Qian: "when there are more human genes than mouse genes in an orthologous set, the mean expression level per gene tends to be lower in human than in mouse, and vice versa". Lan: "the median summed expression of the human duplicates is close to the expression of the singleton orthologs in macaque (median expression ratio 1.11; this is significantly less than the 2:1 expression ratio expected on the basis of copy number" | Level: each mouse member's within-species PT percentile sits below the human singleton's within-species percentile, while `Sum_m` sits near it. Shape: `D(Sum_m, human) <= min_i D(m_i, human)` and `Sum_m` within c. Best test families: Slco1a 6:1, Cyp2j 8:1, Ces2 7:1, Serpina3 10:1, Sult2a 8:1 (mouse expanded); and in reverse, Sult1a 1:4 and Mt 2:10 (human expanded) | Individual mouse members at or above the singleton percentile, or `Sum_m` further from the human shape than the best single member |
| D04 | Subfunctionalisation of the expression domain is rare; conservation is commonest, neofunctionalisation second | Force 1999 Genetics; Lan 2016 Science; Assis 2015; Huminiecki 2004; Farre 2010 Mol Biol Evol | 10101175; 27199432; 26173681; 15466287; 19822635 | yes-abs; yes-full; yes-full; yes-abs; yes-abs | Force: "the usual mechanism of duplicate gene preservation is the partitioning of ancestral functions rather than the evolution of new functions". Lan: "For duplicates with dS < 0.7, just 15.2% of duplicates are classified as potentially sub- or neofunctionalized in expression". Assis: "382 cases of conservation, 213 cases of neofunctionalization ..., 9 cases of subfunctionalization, and 50 cases of specialization". Huminiecki: "only one case of subfunctionalization". Farre: "around 25% of the families, shows patterns compatible with tissue-expression partitioning" | Apply the Assis rules to PT shapes with the human gene as outgroup O and cutoff c from B. Conserved: both `D(m_i,O) <= c`. Neo: exactly one > c. Sub: both > c but `D(m_i+m_j, O) <= c` (a complementary S1/S3 split recreating the human shape). Specialisation: all three > c. Pre-registered expectation: conserved is the plurality and sub is <= 25% of classifiable pairs | Subfunctionalisation is the commonest class, or exceeds 25% of pairs |
| D05 | Asymmetric divergence: the major paralog keeps the ancestral profile, the minor one becomes more specific | Kryuchkova-Mostacci 2016 PLoS Comput Biol; Farre 2010; Lan 2016 | 28030541; 19822635; 27199432 | yes-full; yes-abs; yes-full | Kryuchkova: "the higher expressed paralog has a stronger correlation with the outgroup, thus appears to keep more the ancestral tissue-specificity, while the lower expressed paralog has a lower correlation and appears to become more tissue-specific". Farre: "broad expression of one of the gene copies and reduced, overlapping, expression of the other copy or copies" | Within each mouse family the major member (highest PT mean) has smaller `D` to the human gene and lower `tau_axis` than the minor members, in most families (sign test across families) | No relation between expression rank and D or tau_axis, or the same effect appears in expression-matched 1:1 genes (i.e. it is a low-count artefact) |
| D06 | Tandem duplicates are co-regulated, so they share the axis profile | Lan 2016 Science; Ibn-Salem 2017 Nucleic Acids Res (TADs); Eraly 2003 BBRC (SLC22 pairs); Ghanbarian 2015 Mol Biol Evol; Almeida-Silva 2025 Plant Cell | 27199432; 27634932; 12504088; 25743543; 41092097 | yes-full; yes-full; yes-abs; yes-abs; yes-abs | Lan: "Gene expression is more correlated for tandem duplicates than for singleton neighbors (P = 10-19; t test)"; "Separated duplicates have roughly 50% lower correlation of expression across tissues". Ibn-Salem: "paralog gene pairs are enriched for co-localization in the same TAD, share more often common enhancer elements than expected". Eraly: "pair-members have similar tissue distributions". Almeida-Silva: "Duplication mechanisms that preserve cis-regulatory landscapes typically generate paralogs with more preserved expression profiles" | Shape correlation among mouse paralogs in the same genomic cluster exceeds (a) expression-matched random PT gene pairs, (b) distance-matched non-paralog neighbours in the same cluster, (c) separated paralogs of similar dS. Should hold in both species for tandem pairs such as SLC22A6/A8 and Slc22a6/a8 | r(tandem paralogs) not above r(non-paralog neighbours) or r(random PT genes) |
| D07 | Competing hypothesis to D06: neighbourhood, not paralogy, explains co-expression | Shoja 2007 Comp Funct Genomics (361 human / 212 mouse tandem arrays); Ibn-Salem 2017 | 18273387; 27634932 | yes-abs; yes-full | Shoja: "both chromosomal proximity and expression correlation in TAGs do not differ significantly from their neighboring non-TAG gene pairs". Ibn-Salem: "the positive correlation of gene expression over different tissues is not significantly higher for paralogs whether they are in the same TAD or not" | r(tandem paralogs) approximately equals r(non-paralog neighbours) and both exceed r(random). Pre-register the (paralog minus neighbour) effect size as the D06-vs-D07 contrast | A clear paralog-specific excess over neighbours (supports D06 instead) |
| D08 | Ortholog conjecture for expression: 1:1 orthologs are more similar than within-species paralogs | Chen 2012 PLoS Comput Biol; Kryuchkova-Mostacci 2016; Altenhoff 2012; Nehrt 2011 (contra); Studer 2009 | 23209392; 28030541; 22615551; 21695233; 19368988 | yes-abs (all) | Chen: "the expression similarity between orthologs is significantly higher than that between within-species paralogs, supporting the ortholog conjecture and refuting the cellular context hypothesis for gene expression". Kryuchkova: "Paralogs in human show much lower conservation, even for recent Primate-specific paralogs". Altenhoff: "orthologs have generally more similar functional annotations than paralogs". Nehrt: "paralogs are often a much better predictor of function than are orthologs" | For ancient paralog pairs that are 1:1 in both species (SLC5A2/SLC5A1, SLC22A1/SLC22A2): `D(A_mouse, A_human) < D(A_mouse, B_mouse)`; more generally B sits below the within-species paralog-pair D distribution | Within-species paralogs systematically closer than cross-species orthologs - in this design that most likely indicates a panel/coordinate artefact and must be resolved before any biological reading |
| D09 | Young or lineage-specific duplicates are more axis-specific, and specificity rises with family size | Huminiecki 2004; Farre 2010; Huerta-Cepas 2011 Brief Bioinform; Guschanski 2017 | 15466287; 19822635; 21515902; 28743766 | yes-abs; yes-abs; yes-abs; yes-full | Huminiecki: "decreased breadth and increased specificity of expression as gene family size increases". Farre: "gene duplication is significantly associated with a reduction in expression breadth and intensity". Huerta-Cepas: "gene duplication leads to increased levels of tissue specificity and that this tends to occur promptly after the duplication event". Guschanski: "Lineage-specific and young duplications also contained higher than expected proportions of genes specifically expressed in the liver" | `tau_axis` of mouse lineage-specific paralogs exceeds that of their human co-ortholog(s) and of expression-matched 1:1 genes; tau_axis and the fraction of members absent from PT both increase with mouse family size | No tau_axis difference after expression matching, or the human singleton is more axis-specific than the mouse members |
| D10 | When two young small-scale paralogs are both tissue-specific they are specific to the **same** tissue; ohnologs differ | Kryuchkova-Mostacci 2016 | 28030541 | yes-full | "In the rare case where both paralogs are tissue-specific, small-scale young paralogs are expressed in the same tissue, while genome-wide old paralogs (ohnologs) are expressed in different tissues" | When two or more young mouse paralogs are both segment-restricted they peak in the **same** segment more often than chance (chance set by the marginal segment distribution of segment-restricted PT genes); ancient pairs (SLC5A1/2) split | Young segment-restricted pairs peak in different segments at or above the chance rate |
| D11 | Family-level profiles are conserved across species despite member turnover; lineage-specific expansions are the exceptions | Guschanski 2017 | 28743766 | yes-abs + full | "the overall spatial expression profiles of gene families tend to be conserved, with several exceptions of pronounced tissue specificity shifts among lineage-specific gene family expansions"; "lineage-specific expansions of the same gene family can evolve unique expression patterns in different lineages" | `D(Sum_mouse family, Sum_human family)` lies within c for most families; the minority outside c concentrate in the most expanded families (Serpina3 10:1, Cyp2j 8:1, Slco1a 6:1, Sult2a 8:1) | Family sums mostly outside c, or deviations that do not track the expansion ratio |
| D12 | Xenobiotic families evolve by birth-and-death in unstable clusters; rodent clusters are expanded and members diverge in tissue distribution | Nelson 2004 Pharmacogenetics; Thomas 2007 PLoS Genet; Emes 2003 Hum Mol Genet; MGSC 2002 Nature; Holmes 2010 Mamm Genome; Mackenzie 2005; Velica 2009 Mol Cancer | 15128046; 17500592; 12651866; 12466850; 20931200; 16141793; 20003443 | yes-abs; yes-abs + full; yes-abs; yes-abs; yes-abs; yes-abs (count only); yes-full | Nelson: "The seven CYP clusters are greatly expanded in the mouse with 72 functional genes versus only 27 in the human". Thomas: "Nearly all unstable CYP450 genes in the mouse and human genomes reside in a few dense gene clusters, forming unstable gene islands that arose by recurrent local gene duplication"; CYP2J "includes multiple genes in most species but only a single gene in humans and macaques". Holmes: "Tandem duplications have added more genes for some families, particularly for mouse and rat genomes". Emes: mouse-specific clusters "appear to be under diversifying positive selective pressure". Velica: "Tissue specificity was also observed for Akr1c21, detected only in the kidney" | Cross-species agreement of family sums (the D11 metric) is lower for the unstable xenobiotic families (Cyp2d, Cyp2j, Cyp4a, Ugt2b, Akr1c, Ces1, Slco1a, Serpina3, Nat8) than for phylogenetically stable 1:1 PT genes; these families also have a larger fraction of members absent from PT | Unstable families as conserved across species as stable ones, or all members expressed uniformly along the PT |
| D13 | Neofunctionalised mammalian duplicates are over-represented in **male kidney** | Assis 2015 BMC Evol Biol | 26173681 | yes-full | "Relative to single-copy genes, there was an underrepresentation of outgroup genes and an overrepresentation of neofunctionalizad genes with highest expression in male kidney" | Mouse paralogs classified neo/specialised in D04 are enriched, relative to D04-conserved members, for an external list of male-biased kidney genes (Section 2 supplies the list). Sex cannot be tested directly because all specimens are male | No enrichment. Note a positive result is also what androgen-driven expression alone would produce, so it is not evidence for neofunctionalisation per se |
| D14 | Positive control: repeated-unit axes split work spatially, and ancient PT paralogs split S1 vs S3 | Halpern 2017 Nature (liver); Moor 2018 Cell (villus); Ben-Moshe 2019 Nat Rev (review); Karbach 2000 AJP Renal (rat OCT1/2); Vrhovac 2015 Pflugers Arch (human SGLT1/2) | 28166538; 30270040; 30936469; 10997918; 25304002 | yes-full; yes-abs; yes-abs; yes-abs; yes-abs | Halpern: "more than 50% of liver genes are significantly zonated"; also "zonated liver enzyme pairs that are significantly spatially discordant and either operate on the same substrate, produce the same product, or act sequentially". Vrhovac: "hSGLT2 and SGLT1 were localized to the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively". Karbach: "rOCT1 protein was mainly observed in S1 and S2 segments ... rOCT2 protein was mainly expressed in the S2 and S3 segments", but "significant amounts of both rOCT1 and rOCT2 mRNA were detected in S1, S2, and S3" | Must pass before any young-paralog result is read: SLC5A2 peaks on the S1 side and SLC5A1 on the S3 side in both species with `D(SGLT2, SGLT1) > c`. For Slc22a1 vs Slc22a2 expect at most a modest mRNA shift. No liver or villus paper analyses paralogs within one family, so there is no mammalian precedent to borrow | The SGLT pair is not separated in either species - then the axis or the probe data cannot resolve a paralog split, and all "shared profile" results (D04-conserved, D06, D10) are uninterpretable |
| D15 | Calibration: even 1:1 orthologs conserve spatial profiles only partly across species | Brosch 2018 Nat Commun (human liver, 3 microdissected zones, vs mouse) | 30297808 | yes-full | "Our zonation comparison was restricted on 313 genes zonated in humans that featured one to one matching orthologues with high confidence expression values in the mouse dataset. We detected 111 genes with conserved zonation signatures between mouse and human" | Report the fraction of PT-variable 1:1 orthologs within c and the same-direction (S1-high vs S3-high) agreement rate; make every paralog divergence call relative to B, not to zero divergence | If 1:1 agreement is near chance the cross-species paralog tests have no power; then report only within-species tests (D01, D05, D06/D07, D10) |
| D16 | A continuous axis of cell states reveals more reciprocal partitioning of paralogs than discrete type contrasts (plant evidence only) | Coate 2020 Front Genet (Arabidopsis root scRNA, 11,470 pairs); Almeida-Silva 2025 Plant Cell | 33240334; 41092097 | yes-abs; yes-abs | Coate: "More than twice as many pairs showed reciprocal expression bias between cell states than between cell types or between roots and leaves". Almeida-Silva: "most small-scale duplicates diverge asymmetrically, consistent with neofunctionalization" | On the same mouse data, the number of paralog pairs with reciprocal bias along continuous t is at least the number found by discrete S1/S2/S3 label contrasts | Discrete segment contrasts find as many or more reciprocal pairs |

Corrections to commonly mis-cited sources: Emes 2003 is *Hum Mol Genet* (12651866), not Genome Res; Gu 2002 uses yeast
data only; Assis & Bachtrog 2013 PNAS (24101476) is *Drosophila* - the mammalian application with the outgroup
classification rules is Assis 2015 BMC Evol Biol (26173681); Freilich 2005 (15998445) is about phyletic age, not
duplicates, and does not apply.

## Caveats on section 3

- **The evidence comes from across tissues, not within one tissue axis.** Every mammalian duplicate study above compares
  organs in bulk. The zonation papers (Halpern, Moor, Brosch) report division of labour but do not analyse paralogs. The
  only within-tissue duplicate studies found are in plants (D16). Effect sizes along a single epithelial continuum may be
  much smaller than between organs. **No mammalian paper analyses paralog partitioning within one tissue axis** - that
  is the gap this project can fill.
- **Metric and expression-level confounding.** Chen 2012: comparing profile similarity with "either Pearson's
  correlation or Euclidian distances is problematic, because the variation of tissue-specificity of expression among
  genes interferes with such a comparison". Guschanski: expression divergence correlates negatively with expression
  level (rho -0.35 and -0.47). Low-count minor paralogs will look divergent and specific in D05 and D09 unless every test
  has an expression-matched null.
- **Separate panels, separate coordinates.** Cross-species distances carry platform and coordinate effects, which is why
  B (D15) is required and why a "cellular context" outcome in D08 is the expected artefact rather than a finding.
- **The human gene is a proxy, not the ancestor.** Huminiecki: "Often, the expression of both copies of a duplicated
  gene appears to have changed relative to the ancestral state". Several "singletons" are themselves 1:2 (CYP4A11/A22,
  CYP2D6/D7, NAT8/NAT8B), and Ugt2b is 7:7 with no clean singleton at all.
- **Non-independence.** Members of one family share a history, so a 6-10 member family yields few independent contrasts;
  gene conversion in tandem clusters makes dS underestimate age, weakening D01; Assis 2015 excluded tandem duplicates,
  yet nearly all target families here are tandem.
- **Sex and design.** All specimens are male and many rodent xenobiotic family members are androgen-regulated in kidney,
  so male-specific mouse expression can mimic neofunctionalisation against a male human (D13) and cannot be separated
  from species. n = 2 healthy mice and 1 human donor; use healthy controls only for these baseline tests.
- **mRNA vs protein.** Karbach 2000 found OCT1/OCT2 split the PT at protein level while both mRNAs were present in S1-S3.
  An mRNA axis can understate splits that are real in protein.
- **Axis end and aggregation.** Mouse late S3 grades into thin-limb identity, so a "late" peak can reflect composition
  rather than a within-PT program (the DTL1 column in section 1B is included for exactly this check). Tubule aggregates
  average many cells, so aggregate-level co-variation (D06) can arise from mixing rather than co-regulation.
- **The literature itself conflicts on co-regulation** (D06 vs D07), so register them as competing hypotheses with a
  pre-specified effect-size contrast, not as a confirmatory test.

---

# 4. AKI and non-1:1 genes (mouse ischemia-reperfusion)

## 4.1 Gene claims

| id | gene | family-level orthology | IRI direction | segment / location reported | timing | paper | PMID | verified | quote | strength |
|---|---|---|---|---|---|---|---|---|---|---|
| E01 | Havcr1 (Kim1) | class 1:1; FAMILY 4:1 (Havcr1, Timd2/5/6 : HAVCR1) | up | PT; separate "injured S1/2" and "injured S3" clusters | 4 h, 12 h | Kirita 2020 PNAS | 32571916 | yes-full | "\"Injured S1/2\" and \"injured S3\" were primarily composed of cells from the 4- and 12-h time points and expressed Myc ... and Havcr1" | strong |
| E02 | Havcr1 | as E01 | up | outer stripe of outer medulla (S3 territory) | 6 h | Melo Ferreira 2021 JCI Insight (Visium) | 34003797 | yes-full | "The distribution of Havcr1 expression was localized to the outer stripe in the IRI model and not upregulated in the CLP model." | medium |
| E03 | Havcr1 | as E01 | up | **throughout cortex**, not confined to S3 | 4 h - 2 d | Wu 2024 Nat Commun (HybISS) | 38360882 | yes-full | "the spatial distribution of the injPT cell state and the injury marker Havcr1 expression extends throughout the cortex at early timepoints" | medium |
| E04 | Havcr1, Lcn2 | 1:1 (class) | up, transient | whole kidney | peak day 2, near basal by 1 wk | Liu 2017 JCI Insight (bulk) | 28931758 | yes-full | "This group included two ischemic/nephrotic AKI marker genes that dramatically increased at day 2: Havcr1 and Lcn2 (also known as NGAL)" | medium |
| E05 | Lcn2 (Ngal) | 1:1 | up | distal nephron: TAL, macula densa, CD intercalated cells; **absent from PT** | 3-6 h peak | Paragas 2011 Nat Med (reporter, M+F) | 21240264 | yes-full | "Ngal RNA and mC fluorescence was expressed by the Thick Ascending Limbs of Henle (TAL), the macula densa (MD), and the intercalated cells (IC) of the collecting ducts (CD) after ischemia."; "no expression was found in proximal tubules (PT)" | strong |
| E06 | Lcn2 vs Kim1 | 1:1 | up | Lcn2 in loop of Henle/DT/CD; Kim1 in injured PT | UIR day 1 | Rudman-Melnick 2020 JASN | 33115917 | yes-full | "Kim1 is predominantly expressed in the injured prox, whereas Lcn2 labels the loop of Henle, distal tubule, and collecting ducts of UIR day 1." | strong |
| E07 | Lcn2 | 1:1 | up | medulla and papilla > cortex | 12 h - 2 d; baseline by 6 wk | Dixon 2022 JASN (**female**, Visium) | 34853151 | yes-full | "the injury marker Lcn2 was rapidly upregulated in the medulla and papilla and to a lesser degree in the cortex" | medium |
| E08 | Krt20 | 1:1 | up, de novo | LTL+ PT at the cortex/outer-medulla junction (day 3); failed-repair PT (1 wk); atrophic cortical PT (4 wk) | day 3 -> 4 wk | Liu 2017 | 28931758 | yes-full | "In uninjured and sham surgery control kidneys, no Krt20 was detected in any renal tubule segment"; "by day 3 after IRI, Krt20 was readily detected in proximal tubules" | strong |
| E09 | Krt20 | 1:1 | up | early at the cortico-medullary boundary, later the outer cortex | day 7 -> day 28 | Gerhardt 2021 PNAS | 34183416 | yes-full | "the predominant location of these cells shifted from the cortico-medullary boundary (CMB) ... in early injured samples to the outer cortex at 28 d" | strong |
| E10 | Krt8, Krt18 (with Krt20, Havcr1) | 1:1 | up | PT at the cortex/outer-medulla junction | day 3 | Liu 2017 | 28931758 | yes-full | "a sparse thin layer of elongated cells lining the LTL+ segment of the proximal tubule were colabeled by antibodies for Krt20, Krt8/18, and the kidney injury marker Havcr1" | strong |
| E11 | Vcam1 | 1:1 | up, persistent | failed-repair PT (FR-PTC) inside Kim1+ tubules, then atrophic tubules | from 2 d; about 30% of PT at 14 d | Kirita 2020 | 32571916 | yes-full | "Vcam1-positive FR-PTCs emerged within Kim1-positive injured tubule in a scattered manner at 2 d after IRI and then expanded and remained within atrophic tubules" | strong |
| E12 | Slc22a30 (+ Slc5a12, Slc7a13) | **Slc22a30 FAMILY 6:3**; others 1:1 | **down** | FR-PTC lose terminal differentiation markers | persists late | Kirita 2020 | 32571916 | yes-full | "this cluster additionally down-regulated expression of terminal differentiation markers such as Slc5a12, Slc22a30, and Slc7a13 even at late time points" | strong |
| E13 | Slc22a30 | FAMILY 6:3 | (identity marker) | S3 proximal straight tubule | 1-14 d | Balzer 2022 Nat Commun | 35821371 | yes-full | "Two subgroups expressed S3 proximal straight tubule (PST) markers (Slc22a30, Atp11a, Inmt)." | strong |
| E14 | Ccl2, Ccl7, Ccl8 (+ Il34, Cxcl10, Csf1) | Ccl2 class 1:1 but **FAMILY 5:5** | up, late | FR-PTC signalling to leukocytes | rises to 6 wk | Kirita 2020 | 32571916 | yes-full | "FR-PTC signaling to leukocytes included a variety of proinflammatory and profibrotic cytokines, including Csf1, Il34, Ccl5, Tnf, Ccl2, Ccl7, Ccl8, and Cxcl10." | medium (ligand-receptor inference) |
| E15 | Vcam1, Ccl2 | Vcam1 1:1; Ccl2 FAMILY 5:5 | up, late | Vcam1+/Ccl2+ PT at the CMB from early-replicating cells | day 28 | Gerhardt 2021 | 34183416 | yes-full | "At the cortico-medullary boundary (CMB), where injury initiates, the majority of Vcam1+/Ccl2+ PTCs arose from early replicating PTCs." | strong |
| E16 | **Slc7a12** | **mouse-only** | **up (ectopic) in male PT** | male S3 cluster; Vcam1+ injured PT | day 7 / 28 | Gerhardt 2021 | 34183416 | yes-full | "The S3 cluster 8, for example, showed ectopic activation of Slc7a12, a female-specific S3-marker gene, in IRI replicates (23), while also showing expected expression of characteristic male S3 segment-specific genes like Cyp7b1"; "The predicted up-regulation of Slc7a12 in male PTCs after IRI was confirmed by RNAscope in situ hybridization studies" | strong |
| E17 | Cdh6 | 1:1 | up, de novo | early injured PT (normally a DTL marker) | day 7 | Gerhardt 2021 | 34183416 | yes-full | "The early injured clusters showed enriched expression of Havcr1 and the cadherin encoding genes Cdh6 and Cdh13." | strong |
| E18 | Cxcl1, Cxcl2 (+ Spp1, Ccl2, Vcam1, Sox9, Cdh6) | **Cxcl1/2 FAMILY 3:3**; Ccl2 5:5; rest 1:1 | up | damage-associated PT (DA-PT) | 6 h - 21 d | Ide 2021 eLife | 34279220 | yes-full | "The cells in the DA-PT cluster also show enrichment of developmental genes (Sox4, Sox9, Cited2, and Cdh6), damage-induced genes (Havcr1 and Vcam1) and inflammatory cytokines and chemokines (Spp1, Ccl2, Cxcl1, and Cxcl2)." | strong |
| E19 | **Acsm2** (+ Slc34a1, Lrp2, Hnf4a, Miox, Gpx3) | **Acsm2 FAMILY 1:2**; others 1:1 | **down** | DA-PT | 6 h - 21 d | Ide 2021 | 34279220 | yes-full | "cells in the DA-PT cluster show reduced expression of homeostatic marker genes (Acsm2, Slc34a1, and Lrp2), oxidative stress defense genes (Hmox, Miox, and Gpx3), and a gene encodes transcription factor for PT maturation in renal development (Hnf4a)" | strong |
| E20 | Cxcl2, Il1b, S100a8, S100a9 | **Cxcl2 3:3; S100a9 2:1**; others 1:1 | up | PT lineage-2 (maladaptive) endpoint | 14 d after long IRI only | Balzer 2022 | 35821371 | yes-full | "The lineage 2 endpoint was exclusively from samples 14 d after long IRI and showed increased inflammatory (e.g., Il1b, Cxcl2, S100a8, and S100a9) and stress marker (e.g., Junb, Fos) expression." | medium (possible myeloid admixture) |
| E21 | Ccl3, C3, Il1b, Cxcl2 | **Ccl3 FAMILY 1:4**; Cxcl2 3:3 | up | proinflammatory PT subcluster | 1-14 d | Balzer 2022 | 35821371 | yes-full | "we identified a proinflammatory subcluster expressing a variety of chemokines (Il1b, Cxcl2, Ccl3, Tyrobp, and C3)" | medium |
| E22 | Cxcl1, Cxcl2 (+ Lif, Il6) | FAMILY 3:3 | up | **interstitial** (Foxd1-TRAP), not nephron | 24 h | Liu 2017 | 28931758 | yes-full | "A number of cytokine-encoding genes that showed extensive upregulation (Lif, Cxcl1, Cxcl2, Il6) ... were enriched in Foxd1-TRAP analysis of interstitial cell types at 24 hours after IRI" | medium |
| E23 | Saa1, Saa3 | **FAMILY 3:3** (human SAA3P is a pseudogene) | up, transient | whole kidney, cell source unassigned | peak about day 2 | Liu 2017 | 28931758 | yes-full | "15 genes related to wound repair (Kng1, Arg1, Fgg, Dysf, Procr, Fga, Saa1, Fgb, Ccr1, F13a1, Saa3, Pf4, Entpd1, Orm2, and F2r)" | weak (bulk) |
| E24 | Saa3 | FAMILY 3:3 | up | macrophage-derived in sepsis-AKI (not PT) | sepsis model | Peng 2024 Int Immunopharmacol | 38134592 | yes-abs | "we identify Saa3, a gene primarily expressed in macrophages, as a potent pro-inflammatory cytokine in SAKI"; "further exploration of human homologs is warranted" | medium (not IRI) |
| E25 | Saa3 promoter | FAMILY 3:3 | up | injured kidney (reporter bioluminescence), nephrotoxic models | early | Kudo 2024 Int J Mol Sci | 38791162 | yes-abs | "the Saa3-mediated bioluminescent signal was specifically detected in the injured kidney" | medium (reporter, not IRI) |
| E26 | **Kap** | **mouse-only** | **down** | listed among genes with endogenous PT expression | downregulated module | Liu 2017 | 28931758 | yes-full | "identified 34 genes with endogenous expression in the proximal tubules (Table 1 and Supplemental Table 6), including Kap (kidney androgen-regulated protein) and Lrp2" | medium (bulk) |
| E27 | Mfsd2a | 1:1 | **down** | within its normal site, the **S3** segment of LTL+ PT | downregulated module | Liu 2017 | 28931758 | yes-full | "Mfsd2a ... was markedly downregulated within its normal site of expression in the S3 segment of LTL+ proximal tubule" | medium |
| E28 | Slc7a13, Slc22a6, Slc34a1, Lrp2, Cyp2e1, Miox | 1:1 | **down** | PT | UIR day 1 | Rudman-Melnick 2020 | 33115917 | yes-full | "transmembrane transport ( Slc34a1 , Aqp1 , Lrp2 , Slc7a13 , Slc22a6 ), oxidation-reduction ( Slc37a4 , Slc25a2 , Aco1 , Aco2 , Cyp2e1 , Miox )" | medium (juvenile mice) |
| E29 | **Gsta1** | **FAMILY 5:4** (and Gsta1 is PT-absent at healthy baseline, section 1B) | up | PT, **restricted to cortex** | 12 h and 2 d | Wu 2024 | 38360882 | yes-full | "Glutathione S-Transferase Alpha 1 (Gsta1) is a newly identified disease gene whose expression was activated in PT during the injury phase of AKI"; "activated at 12 h and 2d of IRI and restricted to cortex" | medium |
| E30 | Ifi27, Ifitm3 | **Ifi27 FAMILY 3:3; Ifitm3 FAMILY 5:3** | up | PT and mixed-identity injured clusters | UIR day 1 and 4 | Rudman-Melnick 2020 | 33115917 | yes-full | "GO analysis revealed elevation of cytokine response ( Ifi27 , Ifitm3 , Anxa2 , Lcn2 )" | medium |
| E31 | Spp1, Krt7/8/18 | 1:1 | up | widespread tubular, stromal, immune; medulla and cortex | day 1-4 | Rudman-Melnick 2020 | 33115917 | yes-full | "UIR caused widespread Spp1 elevation in both medulla and cortex." | medium |
| E32 | Atf3 | 1:1 | up | **S3 in the outer stripe**; neutrophil-colocalised | 6 h | Melo Ferreira 2021 | 34003797 | yes-full | "The expression of Atf3 specifically colocalized to the outer stripe of the medulla (Figure 7F)." | medium |
| E33 | Atf3, Klf6 (+ JUN, FOS, EGR1 ...) | 1:1 | up, immediate early | whole kidney; **conserved in human transplant biopsies** | 2-4 h | Liu 2017 | 28931758 | yes-full | "15 of the 21 transcription factors in the mouse early response module I (71.4%) were identified in human kidney graft biopsies: JUN, FOS, BTG2, ZFP36, FOSB, EGR1, KLF4, KLF6, JUNB, CSRNP1, ATF3, JUND, NFIL3, DDIT3, and MAFF." | strong (for conservation) |
| E34 | Sox9 | 1:1 | up | PT | rapid, persistent where damage unresolved | Kumar 2015 Cell Rep | 26279573 | yes-abs | "Translational profiling after AKI revealed a rapid upregulation of Sox9 within proximal tubule (PT) cells, the nephron cell type most vulnerable to AKI." | strong |
| E35 | Cryab | 1:1 | up in S3, down in papilla | S3 vs papilla, opposite directions | injury vs repair | Dixon 2022 (female) | 34853151 | yes-full | "in the S3 segment, it was upregulated with injury and downregulated with repair" | medium |
| E36 | **Sprr2f** | class **mouse-only**; FAMILY 6:6 | up several hundredfold (UUO, not IRI); KO worsens IRI | not segment-resolved | UUO day 5 | Huynh 2020 AJP Renal (**female**) | 33017192 | yes-abs | "expression of small proline-rich region 2f (Sprr2f) ... is increased several hundredfold after renal injury using a unilateral ureteral obstruction (UUO) mouse model" | weak for IRI/segment |
| E37 | Metallothionein (isoform not resolved) | **Mt 2:10** | up | tubular epithelium, predominant in **outer medulla** | 3 h -> 24-48 h max | Takahashi 2001 (**rat**) | 12760485 | yes-abs | "MT predominated in the outer medulla in the IARF rat kidney at transcriptional and translational levels" | medium (rat) |
| E38 | **Serpina3n** | class 1:1; **FAMILY 10:1** | no kidney IRI segment data found | - | - | (see note) | 31316093 | yes-full | the closest primary kidney data are rat AKI-to-CKD: "An early and progressive elevation of urinary SerpinA3 (uSerpinA3) was observed during the AKI to CKD transition together with SerpinA3 relocation from the cytoplasm to the apical tubular membrane in the rat kidney" - cortex-level, no S1/S2/S3 | weak |
| E39 | (S3 vulnerability, classic) | - | S3 cell death; S1/S2 recover | S3 selective | 5 min - 48 h | Venkatachalam 1978 Kidney Int (**rat**) | 682423 | yes-abs | "S1 and S2 cells alike are only reversibly injured and recover completely to normalcy within 4 hr, whereas S3 cells selectively undergo progressive cell injury and death" | strong (rat) |
| E40 | (S3 loss, spatial) | - | PTS3 cell loss at the CMB | PTS3 strongly reduced | 4 wk after mild IRI | Polonsky 2024 Nat Commun (seqFISH) | 39237549 | yes-full | "At the cortico-medullary boundary ... the proximal tubule segment 3 (PTS3) and vascular endothelium (Vasc_2) were strongly reduced" | strong |
| E41 | (S3 fibrosis) | - | - | fibrosis restricted to outer medulla next to damaged S3 | 2 wk - 12 mo | Liu 2017 | 28931758 | yes-full | "fibrosis was restricted to the outer medullary region adjacent to damaged S3 segments" | medium |
| E42 | HUMAN: PT/S3 depletion, injured-S3-like state | - | - | human AKI PT (especially S3) depleted; PT-New 2 (hypoxia) resembles mouse injured S3 | post-mortem AKI | Hinze 2022 Genome Med | 36085050 | yes-full | "Profound depletion in AKI was observed within cells of the PT (in particular those representing the S3 segment)"; "PT-New 1 (oxidative stress) showed similarity to injured mouse S1/2 cells, whereas PT-New 2 (hypoxia) resembled injured S3 cells" | medium (critical illness) |
| E43 | HUMAN: SERPINA1, VCAM1, LCN2, HAVCR1 | mouse Serpina1a-e : SERPINA1 **FAMILY 5:1**; others 1:1 | up | SERPINA1 marks PT-New 3 (inflammation), VCAM1 PT-New 4; LCN2 mainly CNT/CD-PC; HAVCR1 mainly PT | - | Hinze 2022 | 36085050 | yes-full | "LCN2 was primarily upregulated in CNT and CD-PC [12, 47], while HAVCR1 was primarily upregulated in PT"; "inflammation response (PT-New 3, marker gene: serpin family A member 1, SERPINA1), and EMT (PT-New 4, marker gene: vascular cell adhesion molecule 1, VCAM1)" | medium |
| E44 | Cyp7b1, **Cyp4a14**, **Slco1a6**, Slc22a7 (sham, sex context) | Cyp7b1 1:1; Cyp4a14 8:2; Slco1a6 6:1 | n/a (healthy) | Cyp7b1 male-enriched; Cyp4a14 female-enriched; Slco1a6 female-enriched in PTs3 | sham | Dixon 2022 | 34853151 | yes-full | "DEGs ( Cyp7b1 , male enriched; Cyp4a14 , female enriched; Slco1a6 , female enriched in PTs3) in male and female sham kidneys" | medium |
| E45 | **Akr1c18**, **Cyp4a12**, Cyp2e1, Abcc3 (sham, sex context) | Akr1c18 9:4; Cyp4a12 8:2 | n/a (healthy) | Akr1c18 and Abcc3 female S3 only; Cyp2e1 and Cyp4a12 male PT only | sham | Wu 2024 | 38360882 | yes-full | "Abcc3 ... and Akr1c18 ... are expressed in the S3 segment of the female PT only, whereas two of the cytochrome P450 superfamily enzymes, Cyp2e1 and Cyp4a12, are expressed in the male PT only" | medium |

## 4.2 The IRI studies

| paper (PMID) | model | sex / strain / age | method | timepoints |
|---|---|---|---|---|
| Kirita 2020 (32571916) | bilateral IRI 18 min | male C57BL/6J, 8-10 wk | snRNA-seq + IF; human allograft deconvolution | 4 h, 12 h, 2 d, 14 d, 6 wk |
| Gerhardt 2021 (34183416) | IRI 18 min | male, 11-19 wk; Krt20-CreERT2, Ki67-CreERT2 | lineage tracing + snRNA + RNAscope | 7 d, 28 d |
| Ide 2021 (34279220) | unilateral IRI 20 min (mild) / 30 min (severe) | male C57BL/6J, 8-16 wk | Drop-seq | 6 h, 1, 7, 21 d |
| Rudman-Melnick 2020 (33115917) | unilateral IRI 30 min | male Swiss-Webster, **4 wk juvenile** | scRNA + RNAscope | 1, 2, 4, 7, 11, 14 d |
| Liu 2017 (28931758) | bilateral IRI 21 min | male C57BL/6NCrl, 9-11 wk | whole-kidney RNA-seq + IF/ISH | 2 h - 12 mo |
| Dixon 2022 (34853151) | bilateral IRI 34 min (female-equivalent) | **female** C57BL/6J | Visium | 4 h, 12 h, 2 d, 6 wk |
| Melo Ferreira 2021 (34003797) | IRI and CLP | 129/SvEv, 8-10 wk; sex not stated in the searched text | Visium + scRNA mapping | 6 h |
| Polonsky 2024 (39237549) | mild bilateral IRI | male C57BL/6J | seqFISH | 4 wk |
| Wu 2024 (38360882) | bilateral IRI 23 min | male C57BL/6J (+ one female sham) | HybISS + Visium + snRNA imputation | 4 h, 12 h, 2 d, 6 wk |
| Balzer 2022 (35821371) | bilateral IRI 23 / 30 min | male C57BL/6 | scRNA-seq | 1, 3, 14 d |
| Paragas 2011 (21240264) | unilateral 15 / 30 min | male and female Ngal reporter mice | reporter imaging, ISH | hours |
| Hinze 2022 (36085050) | **human** AKI, critical illness (4/8 COVID-19) | human, 8 AKI + controls | snRNA + RNAscope | post-mortem |

## 4.3 Takeaways

- **A strict 1:1 filter silently drops injury genes.** Mouse-only: Kap (lost after IRI, E26), Slc7a12 (ectopically
  induced in male S3 after IRI, E16), Sprr2f (class-level mouse-only, E36). Paralog families: Saa 3:3 (E23-E25),
  Cxcl1/2/3 3:3 (E18, E20, E22), Ccl2/7/8 5:5 (E14, E15), Ccl3 1:4 (E21), Gsta 5:4 (E29), Acsm2 1:2 (E19, lost),
  Slc22a30 6:3 (E12, E13 - an S3 identity marker lost in failed repair), S100a9 2:1, Ifi27 3:3, Ifitm 5:3, Mt 2:10 (E37),
  Serpina3 10:1 (E38), Serpina1 5:1 on the human side (E43).
- **Class vs family matters most for Havcr1.** Its own MGI class is 1:1 with HAVCR1, but the family is 4:1 because
  Timd2/5/6 also map to HAVCR1. Ccl2, Slco1a1 and Serpina3n behave the same way. Pre-register which rule is used.
- **The core PT injury and failed-repair markers are 1:1** - Krt20, Krt8/18, Vcam1, Spp1, Sox9, Cdh6, Atf3, Hnf4a, Lrp2,
  Slc34a1, Slc7a13 - so they are safe for cross-species comparison. Lcn2 is 1:1 but distal (TAL/CD) in both species, so
  Lcn2 signal near the S3 end of a PT coordinate more likely reflects neighbouring TAL than PT.
- **"Injury = S3" is time-dependent, not universal.** Early (hours to day 3) the injury is S3/outer-stripe/CMB biased
  (E39, E01, E02, E32, E08, E27). Late (about 4 wk) S3 loss persists at the CMB (E40) but Krt20+ and Vcam1+ cells shift
  to the outer cortex, i.e. S1/S2 (E09, E15). Counter-evidence at early times: Havcr1 runs throughout the cortex (E03)
  and Gsta1 induction is cortex-restricted (E29).
- **Sex.** Nearly all IRI single-cell studies used males (Dixon needed 34 min in females to match 22 min in males). The
  sex-biased PT genes that change context after IRI sit largely in non-1:1 families (Cyp4a14, Cyp4a12, Slco1a6,
  Akr1c18; E44, E45), and the clearest "feminisation" signal in male IRI (Slc7a12, E16) is a mouse-only gene that cannot
  be mapped to human at all.
- **Human counterparts.** Immediate-early transcription factors are conserved mouse IRI -> human transplant biopsies
  (E33); human AKI depletes S3 most and its hypoxia state resembles mouse injured S3 (E42); human failed-repair-like
  states are marked by VCAM1 and SERPINA1, whose mouse orthologs are a 5:1 family (E43). **No paper was found that
  discusses how mouse-only or expanded-family injury genes (Kap, Slc7a12, Saa3, Serpina3n, Mt, Sprr2) translate to
  human.** That is a stateable gap.
- **Cell-source caveats.** Saa, Cxcl5 and Ccl7 come from whole-kidney bulk (Liu 2017); Cxcl1/2 at 24 h were interstitial
  (E22); the S100a8/9 and Il1b endpoint may include myeloid admixture (E20); Saa3 is macrophage-derived in sepsis-AKI
  (E24). PT attribution of these families is a hypothesis, not an established fact.
- **Baseline check from section 1B:** at healthy male baseline Serpina3*, Saa*, Mup*, Gsta1 and Apol* are all
  approximately 0 TPM in PT microdissection, so any PT signal for them in the AKI arm is induction from a true zero,
  not modulation of a resident program.

---

# 5. Human-specific or human-expanded families in PT, and single-donor CNV risk

| id | gene/family | human-specific feature | allele frequency if given | kidney/PT expression and segment | PMID | verified | quote | implication for a single male donor |
|---|---|---|---|---|---|---|---|---|
| H01 | UGT2B17 | whole-gene deletion polymorphism (0-2 copies); 117 kb deletion by NAHR; selection signals | "varies in copy number from zero to two per individual and also differs in mean number between populations from Africa, Europe, and East Asia" | not segment-resolved | 18760392 | yes-abs | "The human UGT2B17 gene varies in copy number from zero to two per individual"; "a 117 kb deletion arising by nonallelic homologous recombination between approximately 4.9 kb segmental duplications" | a zero may be a germline null, not PT biology |
| H02 | UGT2B17 (discovery) | the deletion spans the entire gene, linked to UGT2B15 | differs between African Americans and Caucasians | n/a | 15475248 | yes-abs | "the presence or absence of a 150 kb genomic interval spanning the entire UGT2B17 gene, revealing that UGT2B17 is present in the human genome as a deletion polymorphism linked to UGT2B15" | a whole-gene null cannot be rescued by another probe position |
| H03 | UGT2B17 (phenotype) | del/del abolishes urinary testosterone glucuronide | del/del 66.7% Korean vs 9.3% Swedish men (all male) | urine phenotype | 16332934 | yes-abs | "The deletion/deletion genotype was seven times more common in the Korean (66.7%) than the Swedish population (9.3%)" | prior of a complete null in one donor: about 9% (European) to 67% (East Asian) |
| H04 | UGT2B17 + UGT2B28 | the only UGT genes with common whole-gene deletions | 0/0 measured: UGT2B17 15%, UGT2B28 4% in US controls | colorectal study, no kidney data | 23575887 | yes-full | "UGT2B17 and UGT2B28 are the only UGT genes that exhibit common whole gene deletion polymorphisms"; "The frequencies of the null (0/0) variants in the control population for UGT2B17 and UGT2B28 were 15% and 4%" | about 15-19% chance of a null at one of these two loci |
| H05 | common deletions genome-wide | the UGT2B17-type deletion is a general class, enriched for sex-steroid and drug-response genes | 541 deletion variants; 120 seen homozygous | n/a | 16468122 | yes-abs | "The coding exons of ten expressed genes were found to be commonly deleted, including multiple genes with roles in sex steroid metabolism, olfaction and drug response" | the one-donor problem is broader than the listed candidates |
| H06 | UGT2B15 vs UGT2B17 | >95% sequence identity | n/a | n/a | 15475248 | yes-abs | "Two members of this family, UGT2B15 and UGT2B17, share over 95% sequence identity" | not individually resolvable by sequence-based assays (see section 6) |
| H07 | UGT2B28 | 97% cDNA identity to UGT2B11 | n/a | detected in breast and adipose, **not** kidney; UGT2B11 in liver/breast/kidney | 21679149 | yes-abs | "Human UGT2B28 was clearly detected in the breast and adipose tissue"; "the UGT2B28 has 97% identity with UGT2B11 at the cDNA sequence level" | weak PT candidate; check any kidney signal against UGT2B11 |
| H08 | UGT1A9, UGT2B7, UGT2B17, CES1 | - | n/a | targeted proteomics: "UGT1A9 was expressed in the liver and kidney"; kidney the most variable tissue | 32955894 | yes-abs | "high interindividual variability (>15-fold) was observed for UGT2B17, CES2 (intestine), SULT1A1 (liver), UGT1A9, UGT2B7, and CES1 (kidney)" | >15-fold between-donor spread for these genes; one donor cannot represent human PT level |
| H09 | SULT1A1 | deletion and duplication, 1 to about 5 copies; copy number dominates activity | 5% of Caucasians single copy, 26% >=3 copies; 63% of African Americans >=3 copies | kidney low (H11) | 17189289 | yes-abs | "5% of Caucasian subjects contained a single copy of the gene and 26% had three or more copies, while 63% of African-American subjects had three or more copies"; "the variability in the level of enzyme activity ... was best explained by gene copy number differences" | human SULT1A1 level is partly a copy-number draw |
| H10 | SULT1A3 / SULT1A4 | human-specific tandem duplication; both copies in every sample, both transcriptionally active; mouse has one Sult1a1 | present in all samples tested | not kidney-specific here | 15358107 | yes-abs | "2 SULT1A3 genes-SULT1A3 and SULT1A4-were present on chromosome 16 in all human DNA samples studied ... both SULT1A3 and SULT1A4 can be transcriptionally active" | a human "SULT1A3" measurement sums two identical genes against one mouse gene |
| H11 | SULT1A1, SULT1A3/4 in kidney | - | n/a | "The kidney and lung contained low levels of SULT"; SULT1A3 and 1A4 quantifiable only jointly | 19679676 | yes-full | as quoted | expect weak PT signal; even protein assays merge 1A3/1A4 |
| H12 | CYP4A11 / CYP4A22 | 96% identical pair, both present in all individuals; CYP4A22 lower expressed and lower activity | not a presence/absence CNV | "CYP4A11 is the predominant isoform expressed in kidney" | 12464262 | yes-abs | as quoted in L13 | treat CYP4A11/A22 as one measurement |
| H13 | CYP4A11 + CYP4F2 PT segment | - | n/a | "in only the S2 and S3 segments of proximal tubules in cortex and outer medulla" | 10660572 | yes-abs | as quoted in L12 | a genuine human S2/S3 > S1 prediction from one antibody study |
| H14 | AKR1C1-C4 | human 4-gene cluster, >=84% identity; explicit contrast with a single rodent enzyme | n/a | AKR1C4 liver, AKR1C3 prostate/mammary, AKR1C1/C2 brain; kidney not listed | 10998348 | yes-abs | as quoted in L31 | a 4-vs-9 family problem, not gene-to-gene |
| H15 | AKR1C1 in kidney | - | n/a | Western: "AKR1C1 being observed in only kidney, liver and testis"; IHC at proximal and distal convoluted tubule resolution | 10510318 | yes-abs | as quoted in L32 | a real human PT gene, and the one whose Visium probe is excluded by default (section 6) |
| H16 | MT1 cluster + MT2A | human-expanded: 10 functional + 7 non-functional isoforms vs 2 mouse genes | MT-1A mRNA in only 2 of 6 kidneys | MT immunoreactivity "confined exclusively to the proximal tubules"; MT-1E/1F/1X/2A consistent; MT-1B/1G/1H/3/4 never detected | 8619255 | yes-abs | as quoted in L70 | per-gene human counts are split over 10 genes vs mouse 2 and per-gene presence varies by donor - compare family sums only |
| H17 | GSTM1 | common whole-gene deletion | "28.8% of blacks and 52.1% of whites had 0 copies of GSTM1" (n = 46,983) | kidney-failure cohort, no segment data | 31555322 | yes-abs | as quoted | about 50% prior that one European-ancestry donor is GSTM1-null |
| H18 | GSTT1 | common whole-gene deletion | "lower among Caucasians (10-20%) than among Asians (50-60%)" | no kidney tubule data found | 17000715 | yes-abs | as quoted | 10-60% prior for a null, **and** GSTT1 has no probe in the human Visium v2.0 panel - a double null |
| H19 | GST alpha vs pi in human nephron | mouse/human family sizes differ (GSTA 5:4; GSTP 3:1) | n/a | "The alpha isoenzyme was present solely in the proximal tubular cells ... while the pi form was restricted to the distal convoluted tubules, the thin loop of Henle, and the collecting ducts" | 7843118 | yes-abs | as quoted | alpha-GST is the canonical human PT marker - and mouse Gsta1's probe is excluded by default, so this is a species-artefact hot spot |
| H20 | APOL1 | primate-specific; no 1:1 mouse ortholog (mouse side is the 11:4 Apol family) | n/a | **contested**: proximal tubular epithelia (21997392) vs podocytes/endothelium only, with antibody non-specificity blamed (34138902) | 21997392 / 34138902 | yes-abs (both) | as quoted in L79 | do not pre-register APOL1 as a PT gene; a PT signal is unresolved |
| H21 | SLC22A11 (OAT4) | human gene with **no rodent ortholog**; mouse Oat5 is a functional not genetic analog | n/a | "expressed in the apical membrane of renal proximal tubules" | 38627096 | yes-abs | as quoted in L63 | a strong human-only PT transporter; a species difference here carries no biology |
| H22 | SLC22A24 | orphan human transporter, absent from the MGI homology report | two functional isoforms in human kidney | "selectively expressed in the proximal tubule segment 3 (PT(S3))" | 31553721 | yes-full | as quoted in L60; also "none in rodents" | the best human-only S3 marker candidate, but low-expressed, so one donor may give few counts |
| H23 | SLC22A10 | human **unitary pseudogene** (fixed P220) | fixed, not polymorphic | great-ape orthologs transport estradiol-17beta-glucuronide | 38782905 | yes-abs | as quoted in L62 | any signal is pseudogene transcript; never report it as a human-mouse functional loss |
| H24 | NAT8 / NAT8B | NAT8B arose by primate duplication and is **inactive in humans**; ongoing gene conversion into NAT8 | fixed stop codon | "The almost exclusive expression of NAT8 in kidney and liver" | 20392701 | yes-abs | as quoted in L43 | NAT8 is a good human kidney gene whose paralog is a dead near-identical sequence; mouse side is 1 of 9 |
| H25 | CES1 / CES1P1 / CES1A2 | segmental-duplication region with pseudogene-hybrid genes; genotyping assays often not specific | n/a | kidney CES1 >15-fold interindividual variability (H08) | 29427553 | yes-abs | "The CES1-related pseudogene ... has been implicated in gene exchange with CES1 and in the formation of hybrid genes including the carboxylesterase 1A2 gene (CES1A2)"; "many procedures for CES1, CES1A2 and CES1P1 genotyping appear to lack specificity" | human CES1 counts may be a donor-variable composite |
| H26 | CYP2D6 / CYP2D7 | hybrid-gene conversions plus whole-gene deletion/duplication | CYP2D7::CYP2D6 2.44%; CYP2D6::CYP2D7 5.84%; deletion 0.15%; duplication 5.98% (n = 106,474) | renal expression low | 38372405 | yes-abs | as quoted | low CNV prior and low renal expression: deprioritise CYP2D6 |
| H27 | UGT1A locus | one locus, 4 shared 3' exons | n/a | UGT1A9 a major renal isoform (H08) | 1339448 | yes-abs | as quoted in L25 | any 3'-biased or shared-exon measurement of UGT1A* is a locus sum, not an isoform |
| H28 | human atlases | - | n/a | Lake 2023 merged PT-S1/PT-S2; Muto 2021 found S1/S2-vs-S3 unclear; none of Lake 2023, Muto 2021, Hansen 2022 mentions MT1*, UGT*, GSTA1, CYP4A11, AKR1C*, APOL1, SLC22A11 or SLC22A24 | 37468583, 33850129, 35675394 | yes-full | as quoted in L85/L86 | human atlases cannot arbitrate the segment for any Section-5 candidate except SLC22A24 and CYP4A11/CYP4F2 |

Citation corrections made while verifying: SULT1A1 copy number is Hebbring 2007 **Hum Mol Genet** (17189289), not Am J
Hum Genet; the CYP4A11/A22 pair of papers is Bellamine 2003 (12464262) with companion Savas 2003 (12464261); "Ohno &
Nakajin 2009 DMD" (18838504) contains no kidney-specific UGT2B17 statement; no paper was found whose primary claim is
that the mouse genome lacks an APOL1 ortholog (the primate-specificity statement is 41894328); no primary paper was
found for GSTT1 tubule-level kidney expression.

---

# 6. Probe and mapping caveats for near-identical paralogs

| id | source | PMID or URL | verified | what it says | implication |
|---|---|---|---|---|---|
| P01 | Robert & Watson 2015 Genome Biol - 12 quantification methods vs simulated truth | 26335491 | yes-full | "Having many members with identical or close-to-identical sequences, gene families are often enriched for multi-mapped reads; therefore, the results of RNA-Seq quantification depend on the choice of aligner, the choice of the reference, a huge array of parameters, and algorithmic details relating to how multi-mapped reads are handled"; "We test a wide range of techniques and find systematic biases within each of them, resulting in severe underestimation or overestimation of gene expression in hundreds of genes"; proposes that "multi-mapped or ambiguous reads can instead be uniquely assigned to groups of genes" | per-gene counts inside a near-identical family are method artefacts in both directions; the defensible unit is the gene group / family sum |
| P02 | Deschamps-Francoeur 2020 Comput Struct Biotechnol J - review of multi-mapped reads | 32637053 | yes-abs | "Such repeated sequences complicate gene/transcript quantification during RNA-seq analysis due to reads mapping to more than one locus"; "Genes of different biotypes have dissimilar levels of sequence duplication" | duplication level differs by gene class, so global normalisation does not correct a family-wide bias |
| P03 | Srivastava 2020 Genome Biol - alignment vs quantification | 32894187 | yes-abs | "the effect of choosing a different alignment methodology, or aligning reads using different parameters, on quantification estimates can sometimes be large and can affect downstream differential expression analyses as well" | both species' matrices must come from an identical pipeline and version, or a species difference can be a mapping difference |
| P04 | Pool 2023 Nat Methods - missing genes in droplet scRNA | 37697162 | yes-abs | "the observed sensitivity deficits stem from three sources: (1) poor annotation of 3' gene ends; (2) issues with intronic read incorporation; and (3) gene overlap-derived read loss" | an all-zero gene can be a reference artefact; zeros for closely spaced or overlapping paralogs need a reference check |
| P05 | Hallinan 2026 eLife - off-target probe binding in 10x panels (OPT tool) | 42065925 | yes-full | "OPT identified 121 probe target sequences across 37 genes as having off-target binding based on perfect sequence homology"; "while we previously found 37 affected genes using GENCODE, we found 14 when using RefSeq and 23 when using CHESS"; "for some genes, the expression patterns detected by Xenium demonstrably reflect the aggregate expression of the target and predicted off-target genes ... rather than the target gene alone" | probe-based signal for a near-identical paralog can be the **family aggregate**, and the count of affected genes is itself annotation-dependent |
| P06 | Du 2025 Genome Biol - SpatialBenchVisium (poly-A vs probe-based Visium, mouse) | 40156041 | yes-full | "Reads mapped with more than 85% confidence to the probe set in our probe-based experiments ... probe sets are highly specific to known sequences in the reference transcriptome"; "For these probe-based samples, only genes that are both expressed and within the probe set can be detected."; "In total, 38% of 19,465 genes targeted in Visium Mouse Transcriptome Probe Set v1.0 were detected" | probe-based data are bounded by panel membership: a gene absent from the panel is a structural zero, exactly like an unfed ortholog column |
| P07 | 10x Genomics Technical Note CG000621 Rev A, "Custom Probe Design" (PDF) | cdn.10xgenomics.com/image/upload/v1674514995/support-documents/CG000621_CustomProbeDesign_RevA.pdf - docs (no PMID) | yes-full | "If probes can bind to sequences other than the target mRNA sequence, off-target signal may be observed ... Matches to off-target genes should have at least five mismatches in at least one of the LHS or RHS probes to prevent efficient hybridization."; "10x Genomics probe panels consist of three pairs of probes for most target mRNAs, with each probe containing 25 bps complimentary to the target mRNA sequence" | the vendor's specificity criterion is >=5 mismatches in one 25 bp arm; genes >95% identical to a paralog have few qualifying windows, which is why such genes carry 1 probe instead of 3 |
| P08 | Visium **Human** Transcriptome Probe Set v2.0 reference CSV (GRCh38-2020-A) | cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv - docs (no PMID) | yes-full | file audit: 18,536 non-deprecated gene ids, 54,018 probes, **514 probes (0.95%) with `included=FALSE` covering 451 genes (2.43%)**; every one of those genes has **all** its probes excluded; 562 further rows carry `DEPRECATED_` ids and are always FALSE | a gene flagged for off-target activity becomes an all-zero row under default filtering, indistinguishable from "not expressed in PT" |
| P09 | Visium **Mouse** Transcriptome Probe Set v2.0 reference CSV (mm10-2020-A) | cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv - docs (no PMID) | yes-full | file audit: 19,405 gene ids, 55,538 probes, **351 probes (0.63%) with `included=FALSE` covering 346 genes (1.78%)**; no DEPRECATED rows | the human panel loses about 1.4x more genes to off-target filtering than the mouse panel (2.43% vs 1.78%): the two species' matrices are **not** censored symmetrically |
| P10 | human v2.0 off-target activity CSV | cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.offtarget.csv - docs (no PMID) | yes-full | 477 probe rows; rows for families in scope: UGT2B15 -> UGT2B17, UGT2B10 and pseudogenes; AKR1C1 -> AKR1C2; GSTM2 and GSTM4 -> GSTM1 (and GSTM5/pseudogenes) | 10x itself declares UGT2B15/UGT2B17/UGT2B10, AKR1C1/AKR1C2 and GSTM1/M2/M4 non-separable - three of the Section-5 candidate families |
| P11 | mouse v2.0 off-target activity CSV | cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.offtarget.csv - docs (no PMID) | yes-full | 346 probe rows; rows include Gsta1 -> Gsta2 (+ Gm10639, Gm3776); Apol9b -> Apol9a; Nat8f6 -> Nat8f3, Nat8f7; Nat8f7 -> ...; Ugt2a2 -> Ugt2a1 | **mouse Gsta1 is excluded by default while human GSTA1 keeps 3 probes** - the canonical alpha-GST PT marker is censored on the mouse side only, so a "human > mouse" alpha-GST result would be a panel artefact |
| P12 | per-gene probe audit of the two v2.0 panels for the candidate families | same URLs as P08/P09 - docs (no PMID), own tabulation of the `included` column | yes-full | **Human, no probe at all:** SULT1A3, SULT1A4, GSTT1, UGT1A6, UGT1A7, NAT8B, CYP2D7, CES1P1, MT1L. **Human, single probe and that probe excluded (all-zero by default):** UGT2B15, AKR1C1, GSTM2. **Human, one included probe:** SULT1A1, SULT1A2, MT1E, MT1G, MT1H, MT1M, MT2A, UGT1A3/1A4/1A5/1A8, UGT2B10, CYP4A22, AKR1C2, GSTM1, CYP2D6. **Human, 2 probes:** UGT2B17, UGT2B7, UGT2B28, UGT2B11, CYP4A11, NAT8, MT1A, MT1B, MT1F, UGT1A10. **Human, 3-4 probes:** UGT2B4, UGT1A1, UGT1A9, CYP4F2, CYP4F11, AKR1C3, AKR1C4, MT1X, GSTA1, GSTA2, GSTP1, GSTM3, APOL1-4, SLC22A6/8/11/12/24/25/10/9, CES1, CES2. **Mouse, single probe excluded:** Gsta1, Apol9b, Nat8f6, Nat8f7; otherwise 2-3 probes for every family member audited | sensitivity per gene is uneven inside a family and asymmetric between species (MT2A 1 probe vs mouse Mt2 2 probes; human MT family 14 included probes over 9 genes with MT1L absent, vs mouse 5 probes over 2 genes). Any per-gene human-vs-mouse level claim in these families confounds probe count with expression |
| P13 | Visium Human Transcriptome Probe Set **v2.1** (the Visium HD human panel) | 10xgenomics.com/support/spatial-gene-expression-hd/.../visium-human-transcriptome-probe-set-v2-1 - docs (no PMID) | **unverified** - the page returned HTTP 429 and support.10xgenomics.com 403 on every attempt; figures come from a web-search summary, not read text | reported there: v2.1 targets 18,530 genes with 54,018 probes, of which 398 genes (2.1%) targeted by 499 probes are excluded by default for predicted off-target activity | the v2.1 probe total equals the audited v2.0 non-deprecated count, so the P08/P10/P12 audit should transfer - but confirm the panel version actually used before quoting these numbers |
| P14 | Space Ranger probe filtering behaviour | 10xgenomics.com/support/... probe-sets pages - docs (no PMID) | **unverified** (HTTP 429) | reported: probes with predicted off-target activity are excluded by default and marked FALSE in `included`; `--filter-probes=false` retains all probe UMIs; `DEPRECATED`-prefixed probes are always excluded | the `included=FALSE` / `DEPRECATED` behaviour is directly observable in the CSVs (P08-P11); only the flag wording and version-specific counts remain unverified |

## Practical caveats

1. **Analyse these families as family sums, not per gene.** Three independent lines agree: the vendor's own off-target
   tables (P10/P11), the >=5-mismatch design rule (P07), and the short-read literature recommending group-level
   assignment (P01, P05). Pre-register the MGI family component as the unit and treat any single-member claim as
   exploratory.
2. **Check panel membership and the `included` flag before calling anything absent.** In human v2.0, SULT1A3, SULT1A4,
   GSTT1, UGT1A6/1A7, NAT8B and CYP2D7 have no probe at all, and UGT2B15, AKR1C1 and GSTM2 have one probe that is
   excluded by default. A zero for these is a structural zero, exactly like an unfed cross-species ortholog column.
3. **Censoring is asymmetric and hits a canonical PT marker.** Human 2.43% vs mouse 1.78% of genes excluded by default,
   and mouse Gsta1 - the alpha-class GST that human IHC places "solely in the proximal tubular cells" - is one of the
   four mouse genes zeroed while human GSTA1 keeps three probes. Any human-vs-mouse difference in these families must
   survive a probe-count and panel-membership control before being called biology.
4. **One human donor cannot support a presence/absence or level claim for a CNV/null gene.** Priors for a single
   European-ancestry male: GSTM1 zero copies about 52%, GSTT1 null 10-20%, UGT2B17 del/del about 9-15%, UGT2B28 null
   about 4%, SULT1A1 single copy about 5% (>=3 copies about 26%), plus donor-variable MT1A detection (2/6 kidneys).
   SLC22A10 is a fixed pseudogene and NAT8B a fixed inactive gene in all humans. Gate these behind a
   "genotype unknown" flag and never report a human-mouse loss for them.
5. **Do not expect human atlases to arbitrate the segment** (H28). Only CYP4A11/CYP4F2 (S2+S3, immunocytochemistry) and
   SLC22A24 (S3, snRNA) have published human segment evidence, and SLC22A24 has no mouse gene at all.

---

# 7. Ranked testable questions for this dataset

Ranked by (i) strength of the prior in the literature, (ii) how much a continuous S1->S3 coordinate adds over S1/S2/S3
bins, and (iii) how interpretable the result is given two male mice, one male human donor and two probe panels.
Each entry: prediction / genes / species and sex / expected direction / PMID support / what counts against it.

**Q1 (control that must pass first). The continuous coordinate resolves the two ancient 1:1 PT paralog splits.**
Genes: Slc5a2 vs Slc5a1 (SGLT2/SGLT1), and Slc3a1/Slc7a9-type pair if present. Species/sex: both species, male.
Direction: SLC5A2 peaks on the S1 side, SLC5A1 on the S3 side, in each species separately, with shape distance above the
1:1 baseline. Support: 25304002 (human IHC, "S1/S2 and S3 segments, respectively"), 26739563 (mouse rBAT S3 vs
b(0,+)AT S1), plus deposited mouse values (Slc5a2 2722/1/0 vs Slc5a1 0/93/39 TPM). Against: the pair is not separated in
one or both species - then the axis or the probe data cannot resolve a paralog split, and every "shared profile" result
below (Q5, Q7, Q12) is uninterpretable.

**Q2 (calibration). What fraction of PT-variable 1:1 orthologs have the same shape in mouse and human?**
Genes: all 1:1 PT-variable orthologs; landmarks Slc34a1, Slc22a7, Lrp2, Slc5a2, Slc5a1. Species/sex: both, male.
Direction: report the distribution of mouse-human shape distance and the same-direction (S1-high vs S3-high) agreement
rate; expect partial conservation, roughly a third by the liver precedent. Support: 30297808 ("We detected 111 genes
with conserved zonation signatures" of 313 human-zonated 1:1 orthologs). Against: agreement near chance - then all
cross-species paralog tests (Q3, Q4, Q6, Q11) have no power and only within-species tests should be reported.

**Q3. Dosage sharing: the summed mouse family matches the human singleton's shape better than any single member.**
Genes: Slco1a 6:1 (vs SLCO1A2), Cyp2j 8:1 (vs CYP2J2), Ces2 7:1 (vs CES2), Sult2a 8:1 (vs SULT2A1), Serpina3 10:1
(vs SERPINA3); reverse direction for Sult1a 1:4 and Mt 2:10. Species/sex: both, male. Direction: `D(Sum_mouse, human)`
<= the minimum over single members, and within the Q2 baseline; levels compared only as within-species percentiles.
Support: 20708291 ("the mean expression level per gene tends to be lower in human than in mouse, and vice versa"),
27199432 ("median expression ratio 1.11"). Against: the family sum is further from the human shape than the best single
member, or individual mouse members sit at or above the human singleton's within-species percentile.
Note the strong negative controls here: SLCO1A2 protein is **distal nephron**, not PT (15632119), and the human members
of the rodent Slc22a19/a26-a30 expansion are **liver** genes (20829431) - so for those two families the expected result
is explicitly "no human positional counterpart", not a matched profile.

**Q4. Mouse paralog families expanded relative to human diverge from the human profile more than stable 1:1 genes.**
Genes: unstable families (Cyp2d, Cyp2j, Cyp4a, Ugt2b, Akr1c, Ces1, Slco1a, Slc22 Oat-like, Nat8, Serpina3) vs
expression-matched 1:1 PT genes. Species/sex: both, male. Direction: lower cross-species family-sum agreement for the
unstable families, with deviation tracking the mouse:human expansion ratio. Support: 15128046 ("72 functional genes
versus only 27 in the human"), 17500592 ("unstable gene islands that arose by recurrent local gene duplication"),
28743766 ("exceptions of pronounced tissue specificity shifts among lineage-specific gene family expansions").
Against: unstable families as conserved as stable ones, or deviation uncorrelated with expansion ratio.

**Q5. Within a mouse family, paralogs co-localise on the axis rather than splitting it.**
Genes: Cyp4a (Cyp4a10, 4a31, 4a14, 4a12a), Nat8f (Nat8f1, f5, f6 + Nat8), Slc22 Oat-like (Slc22a19, a26, a28, a29, a30),
Akr1c (Akr1c14, c18, c21), Ces1 (Ces1d, e, f), Ugt2b (Ugt2b34, b37, b38). Species/sex: mouse, male.
Direction: segment-restricted members of one family peak in the **same** part of the axis more often than the chance
rate set by the marginal distribution of segment-restricted PT genes; sub-functionalising (complementary) pairs are
<= 25% of classifiable pairs. Support: 28030541 ("small-scale young paralogs are expressed in the same tissue, while
genome-wide old paralogs ... are expressed in different tissues"), 27199432 ("just 15.2% of duplicates are classified as
potentially sub- or neofunctionalized"), 26173681 (9 subfunctionalised of 654 pairs). Against: complementary S1-vs-S3
splits are the commonest class, or exceed 25% of pairs. **Known partial counterexamples to pre-register:** mouse Cyp4a
already contains both an S2-peaked member (Cyp4a10/4a14) and an S3-rising member (Cyp4a12a), and Cyp2d contains an
S1-restricted member (Cyp2d26) and an S2/S3 member (Cyp2d9) in the deposited microdissection data - so the expectation
is a bias toward co-localisation, not an absolute rule.

**Q6. A continuous coordinate finds more reciprocal paralog partitioning than the S1/S2/S3 labels do.**
Genes: all mouse paralog pairs in the families above; all human paralog pairs with >=2 PT-expressed members.
Species/sex: both, male. Direction: the number of pairs with reciprocal bias (each member significantly higher somewhere
along t) is at least as large on the continuous coordinate as from discrete segment contrasts; genes whose peak sits at
a segment boundary should be the gain. Support: 33240334 ("More than twice as many pairs showed reciprocal expression
bias between cell states than between cell types"); mechanistically motivated by 31689386's sharp S2-S3 boundary.
Against: discrete contrasts recover as many or more reciprocal pairs - then the continuous coordinate adds nothing here.

**Q7. Tandem-cluster paralogs are more correlated than distance-matched non-paralog neighbours (competing hypotheses).**
Genes: tandem clusters Nat8f1-f7, Slc22a19/a26-a30, Ugt2b34-38, Akr1c12-c21, Ces1a-h, Slco1a1-1a8, Serpina3*; human
counterparts UGT2B*, AKR1C1-C4, GSTM1-M5, MT1 cluster. Species/sex: both, male. Direction: register both D06 and D07 -
report the (paralog minus neighbour) effect size against an expression-matched null. Support for co-regulation: 27199432
("more correlated for tandem duplicates than for singleton neighbors (P = 10-19)"), 27634932 ("enriched for
co-localization in the same TAD, share more often common enhancer elements"). Support for the null: 18273387
("expression correlation in TAGs do not differ significantly from their neighboring non-TAG gene pairs"), 27634932
("not significantly higher for paralogs whether they are in the same TAD or not"). Against either: whichever side the
effect size falls on. **Mandatory artefact control:** probe cross-reactivity (section 6) and tubule-aggregate mixing can
both manufacture co-expression, so the test must exclude genes flagged in the off-target CSVs and compare against
non-paralog neighbours, not against random genes.

**Q8. Human single-copy genes are flatter along the PT than their mouse expanded families.**
Genes: SLCO1A2 vs mouse Slco1a1/1a6; CYP2J2 vs Cyp2j5/11/13; CES2 vs Ces2b/c/g; SULT2A1 vs Sult2a*; SERPINA3 vs
Serpina3*; NAT8 vs Nat8/Nat8f*. Species/sex: both, male. Direction: higher axis-specificity (tau_axis) for mouse family
members than for the human co-ortholog, after expression matching; the mouse family members, not the human gene, carry
the zonation. Support: 15466287 ("decreased breadth and increased specificity of expression as gene family size
increases"), 21515902 ("gene duplication leads to increased levels of tissue specificity"), 19822635 ("a reduction in
expression breadth and intensity"); and the deposited mouse data show steep profiles for Slco1a1, Cyp2j13, Ces2c,
Nat8f5/f6 while human evidence for the corresponding singletons is either flat (UGT1A/UGT2B7 IHC, 17698974) or
non-PT (SLCO1A2, 15632119). Against: no tau_axis difference after expression matching, or the human singleton more
axis-specific than the mouse members. **Confound that must be addressed:** probe count per gene differs between the
panels (section 6, P12), so a flat human profile can be a low-sensitivity artefact.

**Q9. Kap has a bimodal profile in intact males, not a single S3 peak.**
Genes: Kap (mouse-only). Species/sex: mouse, male (both controls). Direction: Kap should be very high toward the S3 end
**and** detectably present in the cortical S1/S2 stretch, because the cortical component is androgen-induced and our
mice are intact males; a pure S3-restricted profile is the castrate/female pattern. Support: 3153475 ("After androgen
treatment, cortical proximal tubules are also induced"), 2739660 ("expressed throughout the cortex ... After castration,
only tubules in the outer stripe of the medulla"), 16150866 (pars recta vs pars convoluta compartments), 12486089
("detected primarily in S3 proximal tubule segments"); deposited male data: 304 / 36,933 / 75,375 TPM across PTS1/2/3.
Against: a strictly S3-restricted profile with no S1/S2 shoulder. This is a good positional test precisely because the
two components are predicted to sit at opposite ends of the same coordinate, which bins cannot express. Kap has **no
human ortholog**, so it is a within-mouse test only.

**Q10. Genes that the literature calls "sex-restricted" by ISH are only quantitatively biased in a male cohort.**
Genes: Cyp4a14 (female-only by RNAscope but 123 TPM in male PTS2), Akr1c18 (female S3 but 956 TPM in male S3),
Slc7a12 (female S3, 11 TPM in male PTS2). Species/sex: mouse, male. Direction: Cyp4a14 should be detectable with an
S2 peak and Akr1c18 with an S3-rising profile in healthy male mice; Slc7a12 should be near-absent at baseline.
Support: 31689386 ("Cyp4a14 (female only) in the S2 subset"), 38360882 ("Akr1c18 ... expressed in the S3 segment of the
female PT only"), vs the deposited male microdissection values in section 1B. Against: Cyp4a14 and Akr1c18 undetectable
in male PT - which would instead indicate a sensitivity or annotation problem, given the microdissection values.

**Q11. The human-only S3 transporter SLC22A24 marks the late human PT, with no mouse counterpart to compare.**
Genes: SLC22A24 (human-only), SLC22A11/OAT4 (human-only), against mouse Slc22a19/a26-a30. Species/sex: human, male.
Direction: SLC22A24 rises toward the human S3 end; SLC22A11 is PT-wide with no published gradient; the mouse expansion
members peak S2/S3 but cannot be aligned gene-to-gene. Support: 31553721 ("selectively expressed in the proximal tubule
segment 3 (PT(S3))" and "none in rodents"), 38627096 ("there is no rodent ortholog of hOAT4"), 32150922 ("a large
rodent-specific expansion, consisting of Slc22a19 and Slc22a26-30"). Against: SLC22A24 flat or absent - plausible,
since the source paper notes "low overall transcriptomic expression in the kidney", so a negative is weak evidence.

**Q12. The Nat8/Nat8f expansion is the cleanest untouched case: no segment data exist for any member in any species.**
Genes: Nat8, Nat8f1-f7 (mouse 9 : human 2, NAT8B inactive); NAT8 (human). Species/sex: both, male. Direction: purely
exploratory - the deposited mouse data suggest Nat8f1 is S1/S2, Nat8 and Nat8f6 S2/S3, Nat8f5 S2 - i.e. partial
splitting inside one tandem cluster. Human NAT8 is male-biased in PT by scRNA (36496458). Support: 20392701 (NAT8B
inactive; "often encoded by multiple, tandemly repeated genes"), 23382457 ("Nat8 was predominant in kidney"), 36496458
(NAT8 male-biased human PT). Against: all Nat8f members with indistinguishable profiles (also the expected signature of
probe cross-reactivity: mouse Nat8f6 and Nat8f7 probes are excluded by default, P11).
This is the highest-novelty question in the list because there is no prior segment-resolved claim to confirm.

**Q13. Mouse Mt1/Mt2 and the human MT1 cluster: compare family sums, never single genes.**
Genes: Mt1, Mt2 (mouse) vs MT1A/B/E/F/G/H/M/X and MT2A (human). Species/sex: both, male. Direction: mouse Mt1 peaks at
the **S1** end (deposited: 248/3/1 TPM) and human MT immunoreactivity is PT-restricted but without a published
subsegment, so the pre-registered test is whether the **human family sum** reproduces the mouse Mt1+Mt2 shape.
Support: 8619255 ("confined exclusively to the proximal tubules"; MT-1A in only 2/6 kidneys; 10 vs 2 genes), 19486827
(rat MT "primarily in the proximal convoluted tubules of the cortex"), 34437439 (human PT MT-I/II), 36758122 (Mt1
female-biased in all three segments). Against: a human per-gene profile that contradicts the family sum, which would
most likely reflect the one-probe sensitivity of MT1E/G/H/M and MT2A and the missing MT1L probe (P12).

**Q14. Injury genes that a 1:1-only analysis would miss are induced along specific parts of the mouse PT axis.**
Genes: mouse-only Kap (expect loss), Slc7a12 (expect ectopic gain in male S3), Sprr2f; families Saa1/2/3, Cxcl1/2,
Ccl2/7/8, Gsta1, Acsm2 (loss), Slc22a30 (loss), Mt1/Mt2, Serpina3n. Species/sex: mouse, male, IRI vs control.
Direction: at the late timepoint expect S3-end loss of Slc22a30, Slc7a13 and Kap with Krt20/Vcam1-type gain extending
into the cortical (S1/S2) stretch; Gsta1 induction cortex-restricted; Serpina3*/Saa*/Mup* rising from a true zero.
Support: 32571916 (FR-PTC lose "Slc5a12, Slc22a30, and Slc7a13"), 34183416 (Krt20+ cells shift "to the outer cortex at
28 d"; "ectopic activation of Slc7a12"), 28931758 (Kap among downregulated PT genes; Mfsd2a lost "in the S3 segment"),
38360882 (Gsta1 "restricted to cortex"), 34279220 (Acsm2 lost in damage-associated PT). Against: injury-associated
change uniform along the axis, or S3-end loss absent. Note Lcn2 must **not** be read as a PT injury gene (21240264,
33115917, 36085050).

**Q15. Sex-biased mouse PT genes are enriched in rodent-expanded families relative to an expression-matched background.**
Genes: the canonical male- and female-biased PT lists from section 2 against family size from the MGI components.
Species/sex: mouse (needs an external sex-bias list because our cohort is all male). Direction: enrichment of expanded
families among sex-biased genes; the honest prior is uncertain, because the one paper that tested gene age found kidney
sex-biased genes are **old** genes with recent sex differences (87% with 1:1 rat orthologs). Support: 37917687 (both the
expansion annotation "the cytochrome P450 family (20/77 genes), the major urinary protein family (16/77 genes) and the
Slc22 transporter family (5/77)" and the contrary age finding), 41889965 (kallikrein-cluster precedent in salivary
gland, preprint), 15177028 and 36758122 (function-level enrichment in drug and organic-anion handling). Against: no
enrichment after matching on expression and function - which is a publishable negative, since **no paper has tested
this in kidney**. Caveat: with all-male specimens, androgen-driven expression and species cannot be separated, so a
positive result cannot be attributed to duplication per se (see D13).

**Q16. The thin-limb end of the mouse axis is a composition change, not a within-PT program.**
Genes: diagnostic DTL markers from the 14-segment data - Aox3 (PT 0, DTL1 75 TPM), Sult1a1 (PT 1-4, DTL1 93), Gstm2
(PT 1-3, DTL1 1668), Fst, Gdf15; versus genuinely S3-rising PT genes Kap, Gstm1, Akr1c18, Cyp4a12a. Species/sex: mouse,
male. Direction: Aox3, Sult1a1 and Gstm2 should appear only in the final stretch of the coordinate and should correlate
with each other, flagging thin-limb admixture; Kap/Gstm1/Akr1c18 should rise before that point. Support: 33769951 (the
14-segment map), 31689386 (Fst marks "thin descending limb just past S3"); deposited values in section 1B. Against:
DTL markers distributed along the whole axis (indicating aggregate-level mixing rather than an axis end) or absent
entirely. This question matters because any "late S3" family claim is otherwise confounded.

**Q17. Are the human PT families flat because of biology or because of the probe panel?**
Genes: the human genes with one or zero included probes (SULT1A1/1A2, MT1E/G/H/M, MT2A, CYP4A22, AKR1C1, AKR1C2, GSTM1,
GSTM2, UGT2B15, GSTT1, SULT1A3/1A4, UGT1A6/1A7, NAT8B, CYP2D7) versus 3-4-probe genes in the same families (GSTA1,
GSTA2, AKR1C3, AKR1C4, MT1X, UGT1A9, CYP4F2, CES1, CES2). Species/sex: human, male. Direction: detection rate and
apparent flatness should track probe count; after restricting to >=2-probe genes in both panels, the human-vs-mouse
family comparison should change materially if the panel is the cause. Support: section 6 P07-P12 (vendor files),
42065925 (probe signal can be "the aggregate expression of the target and predicted off-target genes"), 40156041 ("only
genes that are both expressed and within the probe set can be detected"). Against: no relation between probe count and
detection/flatness - which strengthens every human-vs-mouse claim above.

**Q18. Does the alpha-GST "human PT marker" survive the mouse-side censoring?**
Genes: GSTA1/GSTA2 (human, 3-4 probes) vs mouse Gsta1 (probe excluded by default), Gsta2, Gsta3, Gsta4.
Species/sex: both, male. Direction: mouse alpha-GST activity should appear as **Gsta2** peaking at the S1 end
(deposited: 1006/22/8 TPM) with Gsta1 at zero; the human GSTA1+GSTA2 sum should be PT-wide or S1-weighted. A naive
"human GSTA1 >> mouse Gsta1" result is a panel artefact, not biology. Support: 7843118 and 8372055 (alpha solely in
human PT cells), section 6 P11 (mouse Gsta1 probe excluded, off-target to Gsta2). Against: mouse Gsta1 nonzero (then
the panel version in use does not censor it) or a human alpha-GST profile incompatible with the mouse Gsta2 profile
after family summation.

**Deprioritised deliberately** (recorded so they are not revisited): APOL1 as a PT gene (contested localisation,
21997392 vs 34138902); CYP2D6 (low renal expression, no segment data, CNV); Aox, Sult2a, Akr1b7/b8, Mup, Abcb1a/b
(PT-absent in microdissection, so a PT positional question is not well posed); SLC22A10 (fixed human pseudogene);
NAT8B (fixed inactive, and no probe).
