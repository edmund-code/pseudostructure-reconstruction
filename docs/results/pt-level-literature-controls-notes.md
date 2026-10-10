# Pre-registered positive controls: human vs mouse expression LEVEL in healthy kidney PT

Built **blind to the project's own data**. Nothing under `results/`, `docs/results/`, `docs/paper/`
or any notebook output was opened; no git command was run; no repository file was edited. All
evidence comes from Europe PMC / NCBI E-utilities / Unpaywall / OpenAlex / Semantic Scholar API
responses and from open-access full text. Every PMID was resolved through an API call and checked
against the returned author / year / journal / title.

**Revision 2 (2026-10-10).** Rebuilt to put mouse-specific directions on a firm footing, after the
coordinator asked for mouse rather than "rodent" or "preclinical" evidence.

- **34 controls** (16 high, 18 medium; 24 higher-in-mouse, 10 higher-in-human; 32 gene, 2 pathway).
- **30 of 34 have `direction_verified_for_mouse = yes`.**
- New column: `direction_verified_for_mouse`.

---

## 1. What I could and could not open

| Target | Result |
|---|---|
| **Basit 2019 per-species table** (PMID 31123036, DOI 10.1124/dmd.119.086579) | **UNREACHABLE.** Unpaywall: `is_oa false`, `oa_status closed`, `has_repository_copy false`, zero OA locations. OpenAlex and Semantic Scholar both report closed with no OA PDF. No PMC record. The DOI now redirects to Elsevier (`linkinghub.elsevier.com/retrieve/pii/S0090955624076736`); ScienceDirect returns HTTP 403 behind a Cloudflare challenge, and the old `dmd.aspetjournals.org` URLs return 403. No author manuscript or supplementary file exists in any index I checked. I did not attempt to bypass the access control. |
| **Bleasby 2006 tables** (PMID 17118916, DOI 10.1080/00498250600861751) | **UNREACHABLE.** Unpaywall: `is_oa false`, `oa_status closed`, `has_repository_copy false`, zero OA locations. |
| **Thakur 2024** (PMID 38711199, PMC11218045, *Clin Pharmacol Ther*, DOI 10.1002/cpt.3277) | **OPENED** — NIH author manuscript, full text including **Table 1**. This is a *better* source than either target and is what the revision is built on. |
| Motohashi & Inui 2013 (PMC3675737) | Abstract only; its statement is "rodent", not mouse, so it cannot verify anything under the new rule. |
| Papers citing Basit 2019 (82 found) | Screened for reproduced per-species values; this is how Thakur 2024 was found. |

### The replacement source

Rather than Basit's pooled "preclinical species" statement, **Thakur et al. 2024** quantified **78
membrane drug-metabolising enzymes and transporters** by LC-MS/MS (DIA, DIA-NN, total protein
approach) in kidney membrane fractions of **human (n = 15), rat (n = 10) and mouse (n = 14)** in one
study, and Table 1 reports absolute abundances in **pmol/mg membrane protein with a separate mouse
column**. Every `direction_verified_for_mouse = yes` entry is set from the human and mouse values in
that table, quoted in the `evidence` field, and where the paper also states the mouse comparison in
words (OAT1, URAT1, MRP2, EPHX1, POR, CYP2E1, UGT1A6) the quotation is included.

---

## 2. Directions that CHANGED

Three substantive changes, all because a mouse column replaced a pooled "preclinical" claim.

1. **SLC22A4 / OCTN1 (L03) REVERSED: higher in human -> higher in mouse.** Basit 2019 named OCTN1 as
   one of three transporters *not* higher in preclinical species. Thakur Table 1 gives human 0.6 vs
   **mouse 1.3** pmol/mg, 2.2-fold higher in mouse. The two are partly reconcilable: Thakur's text
   says human OCTN1 is comparable to **rat**, so the pooled "preclinical" average can sit near human
   while mouse specifically sits above it. Two same-method studies from the same laboratory
   nonetheless disagree on the mouse direction, so this is **medium** and is the single entry I would
   treat as least safe.
2. **MATE aggregate (old L18) REVERSED and DROPPED.** v1 predicted higher-in-human apical MATE
   capacity from a gene-count argument: human kidney has MATE1 *and* the kidney-specific MATE2-K,
   mouse only Mate1. Quantitatively this is wrong. Summed MATE protein is 1.2 pmol/mg in human
   (MATE1 1.1 + MATE2 0.1) versus **9.3 in mouse**, because mouse Mate1 alone is 8.5-fold higher than
   human MATE1. The aggregate was deleted and replaced by **L24 (SLC47A1 / MATE1, higher in mouse,
   8.5-fold, high)**. L13 (SLC47A2, higher in human) survives on its own but only as a low-abundance
   call. This is the clearest lesson of the revision: **gene-count arguments do not substitute for
   measured abundance.**
3. **FMO1 (L31) ADDED after being wrongly excluded in v1.** v1 excluded FMO1 because the known
   human/mouse FMO1 difference (Shephard 2007) is a *liver* difference, with both species expressing
   it in kidney. A direct kidney measurement now exists: human 17.0 vs mouse 8.7 pmol/mg, and the
   paper states "~50% lower FMO1 abundance in both the rodent species". Direction: **higher in
   human**.

No other direction flipped; the remaining carried-over entries kept their v1 direction and simply
gained mouse-specific support.

---

## 3. Dropped since v1, and why

| Dropped | Reason |
|---|---|
| **L11 ATP1A1** | Per coordinator instruction. Housekeeping-like in PT, so a library-size or normalisation difference could create or erase it. Also removed from the L16 aggregate gene list for the same reason. |
| **L15 GGT1** | Per coordinator instruction. Rested on a kidney/liver activity *ratio* with human values imported from other papers, not a kidney-to-kidney measurement. |
| **L04 SLC22A5 / OCTN2** | **Fails the mouse-specific bar.** Table 1: human 1.4 vs mouse 1.7 pmol/mg, only 1.2-fold, and the paper states explicitly that "the abundance of human OCTN2 and MCT5 was comparable to the mouse proteins". Basit's 7.4-fold figure was a range across all five species, not human-vs-mouse. Keeping it would have been a predicted-positive that the best data says is null. |
| **L10 SLC5A2 / SGLT2** | **Fails the mouse-specific bar.** Table 1: human 11.5 vs mouse 14.3 pmol/mg, 1.24-fold. Basit's 5.2-fold was again a five-species range. |
| **L08 ABCC3 / MRP3** | Mouse value is 0.3 pmol/mg, below the paper's own 0.5 pmol/mg confirmation threshold, with human ND. Too weak to pre-register. |
| **L18 MATE pathway** | Direction reversed by the data; see section 2. |

Also worth recording as an explicit **null**: **ABCB1 / P-gp**. Thakur Table 1 gives human Abcb1a-
equivalent 1.2 vs mouse 1.0 pmol/mg and states "the average abundance of P-gp in humans is
comparable with rats and mice (Abcb1a)". Basit's "higher in human" for MDR1 is therefore **not**
reproduced in mouse-specific terms. ABCB1 was already excluded from v1 for lack of a 1:1 ortholog;
it now also fails on magnitude, so it should not be resurrected even if the ortholog map collapsed
`ABCB1 -> Abcb1a`.

---

## 4. Caveats specific to the new source

These apply to all 30 mouse-verified entries and matter for how a miss should be read.

1. **Tissue mismatch, the most important one.** Rat and mouse samples are **whole-kidney**
   homogenates ("obtained by homogenizing the entire kidney, ensuring a consistent representation of
   the cortex-to-medulla ratio"). The human samples are "a small section of kidney tissue" whose
   cortex-to-medulla composition was *estimated* from marker proteins (AQP1, claudin 2, KCC3, OAT1
   for PT; AQP2, claudin 7, claudin 10, NCC for distal). Our data is **cortex**, PT-enriched. A
   PT-enriched gene is diluted by medulla in the mouse whole-kidney numbers, so **mouse-higher
   effects for PT genes are if anything underestimated**, while human values carry an unmodelled
   composition uncertainty. `tissue` is recorded as `whole kidney` for these rows rather than
   `cortex`, deliberately.
2. **Protein, not mRNA.** Our data is mRNA. Renal transporter protein/mRNA correlation is imperfect,
   so a failed recovery is weaker evidence against the pipeline than a success is for it.
3. **Table 1 is sex-pooled.** The study found significant sex dependence for **46 of 78 proteins in
   mice** (and 28 in rats) but none in humans. All our specimens are male. Where the paper gives the
   male direction I used it: **Cyp2e1 is 7-fold higher in male than female mice** (L39), and mouse
   `Oct2` (L01) and `Urat1` (L20) are androgen-driven and male-predominant from independent sources
   (Alnouti 2006 PMID 16381671; Cheng 2009 PMID 19679677). For those three the sex-pooled average
   **understates** what male mice should show. For the rest, sex pooling is unmodelled noise in the
   expected magnitude, not in the direction.
4. **ND / BLQ-based directions.** Ten entries rest on "not detected" or "below limit of
   quantification" in one species (L12, L13, L22, L28, L30, L34, L36, L38, L39, L40). A cross-species
   ND can reflect surrogate-peptide selection or sequence divergence rather than true absence. All
   are marked medium except L12 and L34, where an independent source or a very large value in the
   other species supports the call. **L22 (ABCG2 / BCRP)** is the one to watch: human renal BCRP is
   detectable in other work, so the human ND is probably partly methodological even though mouse
   46.1 pmol/mg makes the direction safe.
5. **Low-abundance calls.** The paper flags proteins below 0.5 pmol/mg as needing confirmation. Rows
   where one species falls below that line are marked medium and say so: L06, L09, L13, L23, L26,
   L27.
6. **Wide dispersion.** L40 (CYP4B1, mouse 46.2 +/- 39.6) and L30, L33 have SDs approaching the mean.
   Direction only; do not pre-register a magnitude.
7. **Strain not stated in the main text** (it is in the supplementary file, which I did not retrieve).
   Our mice are C57BL/6; if Thakur used a different strain, strain and species are confounded.
8. **One human donor in our data.** Every control is a species-level claim but our human side is
   n = 1, so donor identity and species are confounded by design. Recovering a control cannot
   distinguish "species difference" from "this donor". A limit of the data, not of the list.

---

## 5. Scope limit on the aggregate control (L16)

L16 is marked `direction_verified_for_mouse = no` because both source papers state the aggregate for
"rodents" / "preclinical species", never for mouse alone. A tally **I derived myself** from Thakur
Table 1 — flagged as my derivation in the `evidence` field, not a claim in the paper — supports it
for transporters but **not** for enzymes:

- **Transporters**, numeric in both human and mouse: **13 of 21 are >= 1.5-fold higher in mouse**, 3
  higher in human, 5 comparable. A further 14 were quantified only in mouse versus only 3 only in
  human.
- **Membrane enzymes**, numeric in both: **4 higher in mouse, 6 higher in human.**

So the expectation is a transporter-specific mouse-high shift, **not** a global one. L16's gene list
was narrowed to transporters accordingly, and ATP1A1 removed. This also means the ten
higher-in-human controls are not a rag-bag: nine of them (EPHX1, EPHX3, FMO1, FMO4, MGST1, MGST2,
MGST3, GUSB, plus SLC22A3) are enzymes, which is exactly where the enzyme tally points. **The
enzyme-versus-transporter asymmetry is itself a pre-registered prediction** and a more
falsifiable one than any single gene.

---

## 6. Family-level candidates with no 1:1 ortholog

Not in the CSV. Any apparent signal from these in a 1:1 ortholog analysis is a mapping artefact, so
they double as a **negative** check on the ortholog map. Thakur Table 1 values given where available.

- **ABCB1** vs **Abcb1a / Abcb1b** (1:2). Now also a measured null, see section 3.
- **UGT1A locus** — human `UGT1A1`/`UGT1A6`/`UGT1A7`/`UGT1A9` vs the mouse `Ugt1a` cluster
  (`Ugt1a1`, `Ugt1a2`, `Ugt1a6a/b`, `Ugt1a7c`, `Ugt1a9`, `Ugt1a10`), a shared-exon cluster with no
  clean 1:1 mapping. Excluded despite strong measured differences that would otherwise qualify:
  **UGT1A9** human 43.7, mouse ND; **UGT1A6** 5-fold lower in mice than humans (p < 0.05);
  **UGT1A7** human 1.8 vs mouse 49.7. These are the highest-value items in this section — if the
  project's ortholog map does assign a 1:1 partner for `UGT1A9` or `UGT1A6`, they become strong
  controls.
- **CES family** — human `CES1`/`CES2` vs mouse `Ces1c`, `Ces1d`, `Ces1f`, `Ces2c`, `Ces2g`.
- **CYP4A** — human `CYP4A11` (12.9 pmol/mg, mouse ND), `CYP4A22` vs mouse `Cyp4a10`, `Cyp4a12a/b`,
  `Cyp4a14`. `Cyp4a12a/b` is the male-specific mouse isoform.
- **CYP4F** — human `CYP4F2` (1.0), `CYP4F3` (3.5), both mouse ND, vs the mouse `Cyp4f` cluster.
- **FMO2** — mouse 3.3 pmol/mg, human ND; the common human allele `FMO2*2` is truncated and
  nonfunctional, so this is a human loss-of-function polymorphism, not a clean expression control.
- **SLC22A11 / OAT4** — human 0.5 pmol/mg, mouse **ND**, confirming no functional mouse counterpart.
  Treat any mouse `Slc22a11` call as suspect.
- **Slco1a1 / Slco1a6 (Oatp1a)** — mouse 3.1, rat 1.9, human **ND**; `Oatp1a3` rat-only. Rodent-
  specific, no human 1:1 partner.
- **Slc22a19 / Oat5** — rodent, and female-predominant in mouse (Cheng 2009), doubly unusable.
- **Slc22a6-like Oat-Pg** — rodent only (mouse 0.8, human ND).
- **Ostalpha / Ostbeta (SLC51A/B)** — mouse only (0.2 and 3.0, human ND); orthology to human
  `SLC51A`/`SLC51B` is 1:1 but the human ND is at the detection floor, so held back as medium-weak.
- **Kap** (kidney androgen-regulated protein) — rodent only, abundant and androgen-driven in male
  mouse PT, no human ortholog.
- **Mup** family / rat alpha-2u-globulin — rodent only.
- **UOX** — functional in mouse, unitary pseudogene in human. **GULO** — functional in mouse,
  `GULOP` pseudogene in human; mouse site is liver.
- **AOX** — human `AOX1` only; mouse `Aox1`, `Aox3`, `Aox4`, `Aox3l1`.
- **NAT8 / NAT8B** — human `NAT8B` is a pseudogene; mouse has `Nat8` plus `Nat8f1-f6`.
- **AKR1B** — mouse `Akr1b7`, `Akr1b8` have no human partner.
- **Ren2** — present in some mouse strains but **not** C57BL/6, which carries `Ren1c` only.

---

## 7. Other candidates still excluded

| Candidate | Why |
|---|---|
| FMO3 | Thakur: human 0.1, mouse **ND**, rat 11.2. Both human and mouse at the floor, so no usable direction. The real difference is human-vs-rat. |
| FMO5, MARC2, CES2C, SLC16A4 (MCT5), SLCO4C1, ABCB1 | Measured mouse-vs-human ratios all < 1.5-fold. Recorded as **nulls**, useful as negative controls. |
| ABCC6 (MRP6), SLC28A1 (CNT1), SLCO2A1 (OATP2A1) | 1.7-10-fold higher in human but every value is below 0.5 pmol/mg. Too weak. |
| SLC15A2 (PEPT2), SLC16A10, SLC29A3, SLC10A1-like BSEP | One or both species at or below the quantification floor. |
| PEPT2 affinity (Song 2017, PMID 27836942) | Km in *Pichia pastoris* transformants, i.e. in vitro heterologous expression, not tissue level. |
| PPARalpha / peroxisomal beta-oxidation | Still **liver-only** evidence; no kidney comparison found. |
| Kim 2024 comparative human/mouse kidney single-cell (PMID 39121855) | **Developmental** (human 10.6-17.6 weeks vs mouse P0) against our healthy adult cortex. |
| Lake 2023 (PMID 37468583), Novella-Rausell 2023 (PMID 37275529) | Single-species each; deriving a cross-species direction would need our own cross-dataset normalisation, which is what the controls are meant to test. Circular. |
| Klotzer 2025 cross-species atlas (PMID 40775269) | Genuine >1M-cell human-plus-rodent atlas, but reports pathway *coordination* and conserved cell states, not per-gene species directions. |
| Hinchman 1990 GGT (PMID 1975172), Lash 1990 KYAT1 (PMID 2139845) | GGT1 dropped (section 3). KYAT1 retained as L14 but human-vs-**rat**, so not mouse-verified. |
| Ripp 1999 FMO3 (PMID 9884308) | Human kidney was never measured. |
| Ito 2006 CYP4 (PMID 16552476) | Localisation and inducibility, not baseline level; and no 1:1 ortholog. |
| Uricase / UOX (Lu 2019, PMID 31118497) | No functional human ortholog, and the mouse enzyme is hepatic. |

---

## 8. Conflicts on record

1. **SLC22A4 / OCTN1** — Basit 2019 (higher in human, pooled preclinical) vs Thakur 2024 (2.2-fold
   higher in mouse). Resolved in favour of the mouse column; flagged medium. **Least safe entry in
   the list.**
2. **SLC22A5 / OCTN2 and SLC5A2 / SGLT2** — Basit reports large cross-species ranges (7.4-fold,
   5.2-fold) but those are maxima across five species; human-vs-mouse is 1.2-fold for both. Dropped.
   Not a contradiction, a scope difference, and a warning against reading five-species ranges as
   pairwise effects.
3. **OAT1 vs OAT3 rank within human kidney** — Motohashi 2002 (PMID 11912245) puts OAT3 highest;
   Hilgendorf 2007 (PMID 17496207) puts OAT1 highest. Unresolved, but it concerns a within-human
   ranking and does not bear on the cross-species directions in L05 and L19, both of which now rest
   on Thakur's mouse column.
4. **ABCG2 / BCRP human ND** — inconsistent with other reports of detectable human renal BCRP;
   direction kept, magnitude flagged as uncertain.

---

## 9. Where to look next

1. **Thakur 2024 supplementary file** — holds strain, sex counts, donor metadata, the cortex/medulla
   marker estimation, and the per-sex abundance figures (S1-S3). Would let the 30 verified entries be
   re-expressed as **male-mouse vs human**, removing caveat 4.3, and would fix the strain question.
   Highest-value next step by a wide margin.
2. **Basit 2019 Table 1-2** and **Bleasby 2006 supplementary tables** via institutional access —
   both confirmed closed to anonymous access. Would let the Basit-versus-Thakur OCTN1 conflict be
   adjudicated and would restore mouse-specific directions for L02 (OCT3) and L07 (MRP1).
3. **Al-Majdoub 2021** *Clin Pharmacol Ther* 110:1389-1400, "Quantitative Proteomic Map of Enzymes
   and Transporters in the Human Kidney" (cited as ref 3 by Thakur) — an independent human-side
   quantification, useful for checking Thakur's human column.
4. **Klotzer 2025** supplementary differential-expression tables (PMID 40775269) — the best candidate
   for adult PT per-gene directions at the transcript level, which would address caveat 4.2.
5. Comparative mammalian organ RNA-seq including kidney (Brawand 2011 *Nature*; Merkin 2012
   *Science*; Breschi 2016 *Genome Biol*) — would give an unbiased, non-pharmacology control set and
   fix the remaining ADME bias of this list.
