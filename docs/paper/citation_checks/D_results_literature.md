# Full-text check of every paper cited in Fig. 5 and Fig. 6 (2026-11-03)

Scope: every published paper cited in support of a claim in
`docs/paper/fig5-human-vs-mouse-along-pt.md` ("Literature and external checks", "The most notable results"),
`docs/paper/fig6-aki-along-pt.md` (literature/artifact-check sections, verdict table, "The most notable results"),
and the supporting detail in `scratchpad/fig5_literature.md`. 66 papers. Each was pursued to full text
(Europe PMC JATS, PMC article HTML, publisher open access, institutional repository, Internet Archive copies of the
publishers' own free versions, J-STAGE); **51 were read in full**, and 15 stayed reachable only as abstracts and are
labelled "(abstract)". Those 15 are all pre-2016 subscription titles (physiology.org, Karger, Springer, SAGE,
Wiley) whose full text was never crawled by any open route.
Nothing in the repository was edited. Downloads are in `fulltext/dl_D/`, converted text in `fulltext/dl_D/txt/`.

**The verdict applies to our sentence, not to the paper.** "PARTLY" usually means the paper is sound and our
attribution overshoots it (wrong segment, wrong species, wrong sex, mRNA where the paper shows protein, a number
that is ours rather than theirs, or a statement that lives in the paper's supplement).

## Summary table

### Fig. 5 — human vs mouse along the PT

| Paper | PMID | Doc / claim (short) | Access | Full text | Verdict | Correction |
|---|---|---|---|---|---|---|
| Vrhovac 2015, Pflugers Arch | 25304002 | fig5 notable 2: "human SGLT2 protein covers S1 and S2, SGLT1 is the S3 transporter"; fig5lit adds "mouse SGLT1 is already high in S2" | Springer, Semantic Scholar, OpenAlex, CroRIS/FULIR, CORE all closed; EPMC abstract | N | PARTLY (abstract) | Human half exact, and "expression of both hSGLTs was sex-independent". The appended mouse-S2 clause is not in this paper — its abstract says rodent locations "largely resembled" the human ones. Drop it or source it to our own mouse data |
| Nagamori 2016, PNAS | 26739563 | fig5 "AGT1 is male-only in mouse"; fig6 "AGT1 is apical S3, male-only; b0,+AT is S1-high" | PMC article HTML | Y | SUPPORTED (AGT1) / PARTLY (b0,+AT) | AGT1 male-only, apical, outer stripe + medullary ray is their own data (8-wk C57BL/6J; S3 inferred from rBAT staining). b0,+AT-S1-high is cited prior work inside the paper, not measured there; their own sex data show b0,+AT equal in male and female. No human data |
| Breljak 2016, Am J Physiol Renal | 27053689 | fig5 "OAT2 protein is S1/S2-dominant in human" | physiology.org 403, repositories unreachable; EPMC abstract | N | SUPPORTED (abstract) | Holds for OAT1–OAT3 jointly, basolateral membrane only, no sex difference detected |
| Shen 2015, Drug Metab Dispos | 25904762 | fig5lit "OAT2 basolateral+apical in human, apical only in rat" | Wayback copy of the archived publisher PDF | Y | SUPPORTED | Monkey matches human. The IHC is **cortex only**, a commercial polyclonal, with no segment markers, no n and no sex or age for the human kidney — so it carries zero axial information |
| Ljubojević 2007, Am J Physiol Renal | 16885152 | fig5lit "Rodent OAT2 is S3 brush border, female > male and strongly androgen-inhibited ... the male mouse under-reports OAT2 at S3" | Wayback of the HighWire full text | Y | PARTLY | Rat arm fully supported (Wistar, mRNA + protein + gonadectomy series). The mouse arm is **protein only**, strain and age unstated, and the paper notes that earlier work found **no sex difference in mouse renal Oat2 mRNA** (they invoke post-transcriptional control) — so on our mRNA readout the "male mouse under-reports" step does not follow. Male mouse S3 is weakly positive, not negative |
| Edwards 1999, Pharmacology | 10352425 | fig5 "NEP polarity flips between rat and rabbit microdissection"; backs notable 1 (neprilysin) | Karger 403; EPMC abstract | N | PARTLY (abstract) | Not membrane polarity: microdissected **enzyme activity**, and it is the axial direction that flips. Rat: only the proximal straight tubule (86 ± 11.3 pmol/min/mm) and glomerulus measurable; rabbit: PCT 70.8 > PST 29.6. Species are rat and rabbit — neither is ours |
| McEvoy 2022, Nat Commun | 36496458 | fig5 "DCXR is S3-high in human single-cell data" | PMC XML | Y | PARTLY | "S3-abundant" is imported from their refs 9/10, so citing McEvoy for DCXR = S3 is partly circular; the paper shows cluster PT2 enrichment on a marker heat map, not an S1→S3 axis. 19 living donors, both sexes, 20,772 PT cells, and **no S2 cluster** |
| Sakaguchi 2024, J Physiol Sci | 38166558 | fig5 "human SLC17A3 falls S1 → S3"; backs notable 1 (NPT4) | PMC XML | Y | SUPPORTED | Table 1: 6.02 / 3.82 / 3.56. Caveats worth carrying: reanalysis of two public datasets, 3 adult male donors, 10,080 cells, mRNA only, and **S2 assigned by exclusion** |
| Cheval 2012, PLoS One | 23056504 | DATA: "human and mouse microdissected S1 vs S3 SAGE (one pooled library per segment)" | PMC XML | Y | SUPPORTED (resource) | Confirmed, with two corrections: **sex and mouse strain are never stated**, so these libraries cannot be called male; each library pools 8 mice or 9 humans; human tissue is tumour-nephrectomy cortex, 59 ± 10 y, against 8–10-wk mice; no S2 |
| Hosoyamada 2004, JASN | 14747372 | fig5lit "Mouse URAT1 is brush-border and male > female" | Wayback of the HighWire full text | Y | SUPPORTED, plus a usable axial detail | 3 male vs 3 female at 10 wk, mRNA male/female 2.3, IHC on males only. The authors say they could not assign PT segments, but the staining is "restricted to the proximal tubule and is not present in the medulla", starting at the exit of Bowman's capsule — a cortical/early-PT statement that argues against an S3 position rather than being silent |
| Sato 2011, Biochim Biophys Acta | 21074513 | fig5lit "rat URAT1 apical" | ScienceDirect 403; EPMC abstract | N | SUPPORTED (abstract) | Rat, no segment resolution |
| Karbach 2000, Am J Physiol Renal | 10997918 | fig5lit "rOCT2 mainly S2/S3" | Wayback of the HighWire full text | Y | SUPPORTED | The axial split is **protein** (basolateral, S2 medullary-ray straight portions + S3 outer stripe, S1 faint); the mRNA sits in S1, S2 and S3. Sprague-Dawley, **sex not stated and no sex comparison**, so the androgen half rests entirely on Urakami 2000. rOCT1 is the mirror case and is itself mRNA/protein-discordant |
| Urakami 2000, FEBS Lett | 10812069 | fig5lit OCT2 "is testosterone-driven ... so the mouse-early crossing is androgen-amplified" | Wiley/Elsevier closed; EPMC abstract | N | SUPPORTED (abstract) | Rat; whole-kidney and renal-slice readouts with **no segment data**, so it supports a level effect, not a crossing position |
| Motohashi 2013, J Pharm Sci | 23630107 | fig5lit "Human OCT2 protein is PT-wide with no segment resolution" | Wiley 403; EPMC abstract | N | SUPPORTED (abstract) | Serial sections of human **cortex only** — write "cortical PT without segment resolution"; outer-stripe S3 was not sampled |
| Zhang 2024, Int J Mol Sci | 38474316 | fig5lit "Rat SGLT5 rises progressively S1→S3" | PMC XML | Y | SUPPORTED | Method is reanalysis of public rat scRNA-seq (3 male + 3 female pools), microdissected bulk RNA-seq and published proteomics — not new ISH/IHC |
| Fukuzawa 2013, PLoS One | 23451068 | fig5lit "no human SGLT5 localization exists ... is BBM function in mouse" | PMC XML | Y | SUPPORTED | The paper states outright that no specific SGLT5 antibody exists; male C57BL/6 knockout plus brush-border vesicle uptake; no localization in any species |
| Levillain & Hus-Citharel 1998, Am J Physiol | 9841492 | fig5lit (Odc1) "activity PCT > PST" | Wayback of the HighWire full text | Y | SUPPORTED (both species), with an mRNA caveat | Microdissected **enzyme activity**, male OF-1 Swiss mice and male Sprague-Dawley rats; the segments are PCT / CPST (S2) / OSPST (S3), so S1 and S2 are not separated, and rat OSPST rests on one rat. Our readout is mRNA, and the paper states that **rat ODC mRNA is highest in the outer stripe**, opposite to its activity; only mouse mRNA follows the PCT pattern |
| Koibuchi 1993, Histochemistry | 8307774 | fig5lit (Odc1) "testosterone-induced mRNA stronger in pars convoluta" | Springer closed; abstract | N | PARTLY (abstract) | This is the **mRNA (in situ)** paper, and it says the "outer portion" of the cortical PT exceeds the "inner portion" — it never names the pars convoluta. Intact male BALB/c, one testosterone dose, read at 14 h |
| Koibuchi 1993, J Endocrinol | 8429279 | same Odc1 claim | Bioscientifica closed; abstract | N | PARTLY (abstract) | This one is **protein (IHC)**: pars convoluta strong, pars recta weak, and only after testosterone. Our sentence merges the two papers and calls both mRNA |
| Koibuchi 1990, Endocrinol Jpn | 2083533 | fig5lit (Odc1) "an earlier IHC put most signal in pars recta" | J-STAGE PDF, pages read as images | Y | SUPPORTED (with caveat) | The authors did not establish whether the stained cells were proximal or distal tubule, so "pars recta" is their own tentative reading; they counted stained cells, not intensity; testosterone did not move the location |
| Kim 2018, Exp Mol Med | 29674665 | fig5lit (Chrna4) "renal nicotinic work is all α7" | PMC XML | Y | SUPPORTED | No α4 anywhere; α7 is localized to proximal vs distal tubule only, with no S1/S2/S3 |
| Brie 2024, J Endocr Soc | 38370444 | fig5lit (Hao2) "canonical sex-dimorphic gene ... liver" | PMC XML | Y | SUPPORTED | Hao2 is **female**-predominant in mouse liver and growth-hormone-controlled, not androgen-controlled; no kidney data |
| Thomasson 2003, Br J Cancer | 14520461 | fig5lit (Lrig1) "No kidney axial literature (only RCC)" | PMC XML | Y | SUPPORTED | Human; 8 normal cortex samples show LRIG1 "most prominent in the tubuli", but no PT segment call |
| Chan 1979, J Histochem Cytochem | 39098 | fig5lit (Gpt) "ALT is the mirror image of AST along the proximal tubule and rises into the straight segment" | SAGE 403, reprint not archived; abstract | N | PARTLY (abstract) | Rat, young and adult, sex not stated, "4 parts of the proximal tubule" rather than S1/S2/S3. "Mirror image" is verbatim; **"rises into the straight segment" is our inference** — the 6-fold rise the abstract gives for proximal straight tubules is a postnatal-development result. ALT is about 100x lower than AST in absolute terms |
| Davies 1997, FEBS Lett | 9237666 | fig5lit (Pter) "Rat PT localization by in situ exists, no segment resolution" | Wiley closed; abstract | N | SUPPORTED (abstract) | Reword the second half to "none reported in the abstract", which is all we can check |
| Lin 2011, Am J Physiol Renal | 21228108 | fig5lit (Fut8) "renal FUT8 literature is injury/fibrosis only, no axial data" | physiology.org 403; abstract | N | SUPPORTED (abstract) | The abstract is cultured HK-2 cells only |
| Wang 2017, Am J Transl Res | 29118904 | same Fut8 claim | PMC | Y | SUPPORTED | Entirely cultured HK-2 cells; both Fut8 papers are cell-culture injury models |
| Ransick 2019, Dev Cell | 31689386 | fig5lit (labelled ABSTRACT) "existence of the male-biased mouse PT-S2/S3 program" | PMC | Y | PARTLY | Upgrade from ABSTRACT to full text. The sex differences run **both ways** (male- and female-biased genes in S2 and S3), the largest is in S3, and all PT clusters split by sex; n = 2 males and 2 females (C57BL/6J) |
| Chen 2021, JASN | 33769951 | DATA: "male mouse microdissected PTS1/PTS2/PTS3 TPM" | PMC | Y | SUPPORTED (definitions flagged) | Male C57BL/6, 6–8 wk; **PTS2 is the straight tubule in cortical medullary rays, not late PCT**; PTS3 is the outer stripe; 5/8/6 samples from 3/5/3 mice. Also: the Cyp24a1 and Slc16a1 TPM triples fig6 credits to Fuchs 2024 are this paper's |
| Chen 2023, JASN | 36758122 | DATA: "mouse PTS1/2/3 TPM, male and female" | PMC | Y | SUPPORTED (independence flagged) | **The male data were downloaded from the 2021 study** and match it exactly, so "Chen 2021 male" and "Chen 2023 male" are the same samples, not a replication. Female n = 6/6/7. Its proteomics is whole-kidney |
| Lee 2015, JASN | 25817355 | DATA: "male rat microdissected S1/S2/S3 RPKM" | PMC | Y | SUPPORTED | Male Sprague–Dawley 200–250 g, 7 samples per segment, median RPKM |
| Limbutara 2020, JASN (KTEA) | 32358040 | DATA: "rat S1/S2/S3 protein copies/cell" | PMC | Y | SUPPORTED | Male Sprague–Dawley; S1 cortical labyrinth, S2 medullary ray, S3 outer stripe; 3 each. S1 detects far fewer proteins (2,734 vs 4,750 in S2), so an S1 zero may be non-detection |
| Lake 2023, Nature | 37468583 | DATA: "human PT-S1/S2/S3 snRNA per donor (17 normal, 12 AKI, 23 CKD)"; fig5lit says 13/10/17 | PMC | Y | SUPPORTED (resource) | Two things to fix: fig5 and fig5lit give **different donor counts** for the same object; and in the paper's **single-cell** data PT-S1 and PT-S2 "could not be accurately resolved" and were kept merged, so per-donor S1-vs-S2 must come from nuclei. Lake's segment labels were themselves checked against Chen 2021, so Lake-confirms-Chen is mildly circular |
| Xiong 2023, Dev Cell | 37673062 | DATA + fig5 "Ar knockout confirms androgen control (Hao2, Sh3gl2, Lrig1)" | PMC | Y | PARTLY | Design matches (Six2-Cre nephron Ar deletion, both sexes, bulk RNA-seq + multiome). But the **text names only Hao2**; Sh3gl2, Lrig1, Ggact and Odc1 rest on supplementary tables, which this pass did not open — say "supplementary tables S2/S3" |
### Two papers our documents do not cite but should (checked in this pass)

| Paper | PMID | Doc / claim (short) | Access | Full text | Verdict | Correction |
|---|---|---|---|---|---|---|
| Schlosser 2020, Nat Genet | 31959995 | fig5 "human SMIT2 (SLC5A11) matches an independent donor" — **no citation given** | PMC article HTML | Y | SUPPORTED (and citable) | The claim has a published source we are not naming: "SLC5A11 showed high, specific expression in the S3 segment of the proximal tubule", from 4,524 adult human kidney nuclei (their ref 34 = Wu 2018 Cell Stem Cell, independent of Lake). Cite it instead of "an independent donor" |
| Lahjouji 2007, Biochim Biophys Acta | 17306760 | relevant to fig5 notable 1 ("SMIT2 ... late in the human PT and early in the mouse") | Elsevier closed; EPMC abstract | N | counterpoint (abstract) | The only mammalian SMIT2 protein localization: in **rabbit**, SMIT2 is cortical and "only SMIT2 is expressed at the luminal side of proximal convoluted tubules" — i.e. early PT, like our mouse limb. So the mouse-early half has a precedent and the human-late half is the new one. (Schlosser cites this same paper as agreeing with human S3, which is a loose reading) |


### Fig. 6 — control vs AKI along the PT

| Paper | PMID | Doc / claim (short) | Access | Full text | Verdict | Correction |
|---|---|---|---|---|---|---|
| Liu 2017, JCI Insight | 28931758 | DATA (GSE98622) "48 h against sham"; the bulk directions behind notable 1, and Cyp24a1 / Agt / Cyp51 / MCT1 | EPMC JATS + GEO sample list | Y | PARTLY (data description) | Male C57BL/6NCrl, 10–12 wk, **bilateral** 21-min clamp, whole kidney, n = 3 at 48 h. **There is no 48-h sham** — shams exist only at 4 h, 24 h and 12 mo, so name the reference used. The paper itself names only Slc7a13 (down) and Havcr1; every other fold change is our reanalysis and should be credited to the dataset |
| Kirita 2020, PNAS | 32571916 | AKI cell states; "Slc16a1 and Cyp24a1 absent from Kirita 2020" | EPMC JATS | Y | PARTLY | Male C57BL/6J, bilateral 18-min IRI, snRNA-seq, 4 h–6 wk. Both genes are absent from main text, legends and the Fig. 2A panel, but the Dataset S1–S3 DEG tables were not opened: say "not named in the main text". Useful: injured S3 is mostly 4–12 h and "had almost disappeared by 2 d" |
| Li 2022, Cell Metab | 36265491 | "PT carrier loss in injured states" | PMC | Y | SUPPORTED | IRI is **unilateral** 22 min, male, n = 2 per time point. The carriers named are Slc34a1, Slc5a12 and Slc7a13 — not our S3 partners — and nothing is segment-resolved |
| Wu 2024, Nat Commun | 38360882 | same; "injury causes downregulation of PT segment-defining markers" | EPMC JATS | Y | SUPPORTED | Male, bilateral 23-min IRI, 4 h/12 h/2 d/6 wk; platform is a 200-gene in-situ panel, so "could not separate PTS1/2/3" is a panel-limited statement |
| Xuanyuan 2025, Nat Commun | 40813851 | "injury spreads from the corticomedullary junction"; "healthy-S3 neighbourhood absent at 12 h–2 d" | EPMC JATS | Y | SUPPORTED | Male, bilateral 18-min IRI, Xenium 300-gene + Visium, two replicates per time point. The junction statement rests on maps (Fig. 2A), not a zone quantification — fig6's "maps only" is right |
| Humphreys 2008, Cell Stem Cell | 18371453 | "Proliferation peaks in outer stripe and S3 at 2 d" | Wayback copy of the publisher's own free PDF | Y | PARTLY | Only S3/outer medulla was counted (Ki67 47.4% at 2 d), there is **no S1/S2 comparison**, and 2 d is the earliest day sampled, so "peak" is "earliest point measured". Both sexes, unilateral and bilateral |
| Humphreys 2011, PNAS | 21576461 | same | PMC HTML | Y | PARTLY | Male, mostly **unilateral** 27-min IRI, sampled 12–48 h only and still rising at 48 h, so no peak is shown; the counted region pools "outer medulla and inner cortex" |
| Witzgall 1994, J Clin Invest | 7910173 | same, marked "abstract only" | JCI site PDF (OCR text layer) | Y | SUPPORTED (rat) | **Drop "abstract only"**. Male Sprague-Dawley **rat**, unilateral 40-min clamp; PCNA peaks at 48 h with S3 > S1/S2 |
| De Chiara 2022, Nat Commun | 36195583 | "Polyploid cells are mostly cortical S1/S2" | EPMC JATS | Y | SUPPORTED | Male, unilateral 30-min IRI; 13.3% in cortex vs 2.7% in outer stripe, "mostly ... S1 and S2"; the segment images are day 30, not day 2 |
| Lazzeri 2018, Nat Commun | 29632300 | "at 2 d only about 13% of Ki-67+ cells divide" | EPMC JATS | Y | PARTLY (number misread) | **13.2% is the share of all tubular cells in S/G2/M**, against 47.1% Ki-67+ — not the share of Ki-67+ cells that divide (that ratio is ~0.28), and S/G2/M is not division. Flow cytometry gives 9.8% S and 2% G2/M |
| Tran 2016, Nature | 26982719 | "PGC-1α and OXPHOS collapse after IRI"; Flvcr2 row "choline depleted after ischaemia" | EPMC JATS | Y | PARTLY (OXPHOS) / SUPPORTED (choline) | Male, bilateral 20-min IRI, 24 h. PGC1α falls, but **no post-ischaemic OXPHOS fall is shown** — the OXPHOS panels compare iNephPGC1α mice with controls. "Deficiency of betaine and choline" in post-ischaemic kidneys is verbatim (Fig. 1g,h) |
| Funk & Schnellmann 2012, Am J Physiol Renal | 22160772 | "PGC-1α and OXPHOS collapse after IRI (Funk 2012, abstract only)" | PMC HTML | Y | PARTLY | **Drop "abstract only"**, and cite it only for respiratory-chain **proteins** (NDUFB8, ATP synthase β, COX I/IV fall in cortex through 144 h): in this paper **PGC-1α protein rose** and **cortical ATP was elevated** at 24–72 h |
| Fujii 2019, JCI Insight | 31723053 | "zone-level ATP loss (abstract only)" | PMC HTML | Y | SUPPORTED (model caveat) | **Drop "abstract only"**. ATP falls in cortex and outer stripe but not inner stripe, by MALDI imaging MS — but the model is **unilateral 10-min** clip, milder than ours, and these are not ATP-sensor mice |
| Arnoux 2025, Kidney Int | 40850621 | "Arnoux 2025 names TAL as most affected", listed beside zone-level ATP loss | Geneva repository, published CC BY version | Y | PARTLY | TAL as the vulnerability hub is right, but the IRI result is a **gene-set score, not an ATP measurement**, and its Fig. 3a shows the electron-transport score **not reduced (slightly positive) in PTS1–PTS3** at 4–48 h. So Arnoux is a counterpoint to "the respiratory chain falls through S3", not support for it |
| Matsuzaki 2008, Drug Metab Dispos | 18180268 | "OCT2 and OAT down after IRI" | Wayback copy of the publisher PDF | Y | SUPPORTED (rat) | **Rat**, male, bilateral 30-min clamp, 48 h, whole kidney: rOCT2 protein markedly down, mRNA to 4%; rOAT1/rOAT3 protein also down |
| Schneider 2007, Am J Physiol Renal | 17244891 | same (OAT1/3) | Wayback copy of the publisher full text | Y | SUPPORTED (rat, female) | **Female** Sprague-Dawley rats, bilateral 45-min arterial clamp, cortex; mRNA down 6–24 h, restored by 72 h; protein recovers by 7 d |
| Fuchs 2024, J Clin Invest | 39688907 | "Cyp24a1 normally S1-restricted (85.9 / 2.4 / 0.0 TPM; Fuchs 2024)" | EPMC JATS | Y | PARTLY | The TPM triple is **Chen 2021**, not Fuchs — recredit it. Fuchs supports only "proximal tubules of the outer cortex" by RNAscope and never names S1/S2/S3; its scRNA PT cluster pools S1–S3; C57BL/6J, 12 wk, both sexes; no AKI model (adenine CKD only) |
| Egli-Spichtig 2018, Front Physiol | 30405444 | "Contradicted in folic-acid AKI (down 17-fold), so model-specific" | EPMC JATS | Y | SUPPORTED | Folic acid 250 mg/kg, 24 h, male and female C57Bl/6J, whole-kidney qPCR, n = 5–8. No IRI arm and no segmental data, so "model-specific" is our inference from one model at one time |
| Christov 2013, Kidney Int | 23657144 | "FGF23 rises in AKI" (beside the Cyp24a1 row) | EPMC JATS | Y | SUPPORTED, but no Cyp24a1 | Models are folic acid, glycerol and human cardiac surgery — **no IRI**; renal Cyp24a1 is never measured and plasma 1,25D is unchanged at 24 h, so this paper cannot carry an FGF23 → Cyp24a1 step |
| Matsusaka 2012, JASN | 22518004 | "Agt mRNA is S2/S3 (Matsusaka 2012)", behind notable 3 ("Agt suppressed in S2 where it is normally highest") | PMC article HTML | Y | NOT SUPPORTED as S2 | Its ISH puts Agt mRNA in the **S3 segment / outer stripe only**; S1 and S2 are where reabsorbed **protein** sits. The "normally highest in S2" half is supported by Chen 2021 microdissection (11.7 / 322.1 / 217.1 TPM), so cite Chen for the position and Matsusaka only for mRNA-in-S3-and-protein-in-S1/S2 |
| Kontogiannis & Burns 1998, Am J Physiol | 9458826 | "cortical Agt mRNA falls after ischaemia (Kontogiannis 1998, abstract only)" | Wayback of the HighWire full text | Y | SUPPORTED | **Drop "abstract only"**. Male Sprague-Dawley **rats**, bilateral 60-min pedicle occlusion; angiotensinogen is a Northern blot of whole cortical RNA, so there is no segmental resolution. Ischaemic/sham 0.42 at 6 h, 0.01 at 24 h, 0.02 at 72 h, 0.82 at 120 h |
| Yan 2022, Front Med | 35646989 | "PPARα pathway collapse after IRI (a mutant-mouse study)" | EPMC JATS | Y | PARTLY | Cyp4a10/Cyp4a14 are down in **Ptpn6+/mev versus wild-type littermates, both after IRI** — not IRI versus sham, and the genotypes do not differ without IRI. 34-min clamp, 24 h, n = 3, sex not reported |
| Naito 2009, Am J Pathol | 19095962 | "Ischaemia induces cholesterol synthesis (Naito 2009, abstract only)", the Cyp51 conflict | PMC article HTML | Y | PARTLY | **Drop "abstract only"** (full text is in PMC), and note the gene is **Hmgcr** (+30–35% mRNA), never Cyp51 — so the "phase-dependent conflict" with Liu 2017 compares two different genes. Male CD-1, 30-min unilateral ischaemia, contralateral kidney as control |
| Balzer 2022, Nat Commun | 35821371 | MCT1 "bulk induction about 3–5× (Liu 2017; Balzer 2022)" | EPMC JATS + supplementary | Y | PARTLY (supplement only) | Slc16a1 appears nowhere in the main text or figures; the number is ours from Supplementary Data 2 (control 6.6, short-IRI d1 20.2 = 3.1×, long-IRI d1 15.5 = 2.3×, n = 2 per injured group), back to baseline by day 14 |
| Yu 2025, Nat Commun | 41350284 | "MCT1 lactate entry drives ferroptosis (Yu 2025)" | EPMC JATS | Y | PARTLY | The knockout is **MCT4 (Slc16a3)**; MCT1 is supported only pharmacologically (AZD3965, siRNA) in human HK-2 cells and by inhibitors in vivo. Mouse allograft transplant + 35-min in-situ IRI; no segmental localization of MCT1 |
| Cai 2025, PNAS | 40956895 | "ZIP8 and ZIP14 are S3 metal transporters, ZIP8 knockout protects in IRI" | EPMC JATS + SI | Y | PARTLY | **ZIP14 is never mentioned**, and the S3 point is SI Fig. S3, a reanalysis of published scRNA data plus Human Protein Atlas IHC, not their own localization. The protective knockouts are a global A391T knock-in and a Six2-Cre kidney Slc39a8 deletion; I/R is nephrectomy + 40-min clamp at 24 h |
| Fujishiro 2012, Metallomics | 22534978 | the real source of "ZIP8 and ZIP14 are S3" | RSC/OUP 403; abstract (+ the same group's open 2019 review) | N | SUPPORTED (abstract) | ISH puts ZIP8 and ZIP14 in PT of the outer stripe, i.e. S3; the transport work uses an immortalized mouse PTEC line that is not segment-specific. This, not Cai 2025, is the citation for ZIP8/ZIP14 = S3 |
| Melo Ferreira 2021, JCI Insight | 34003797 | "Fgb in cortical S1/S2 is known" | EPMC JATS | Y | PARTLY | The statement is about an **injured-PT cluster** rich in Fga/Fgb/Fgg that also carries S1/S2 markers, not about basal Fgb in healthy S1. 129/SvEv mice, bilateral IRI, 6 h; sex not stated |
| Becker 2010, Am J Physiol Renal | 20392797 | uncited, behind "MCT1 is basolateral" | publisher closed; abstract | N | SUPPORTED (abstract) | Basolateral in "proximal tubule cells"; no segment assignment |
| Yanase 2008, Histochem Cell Biol | 18751721 | uncited, behind "MCT1 is ... S1-high" | publisher closed; abstract | N | SUPPORTED (abstract) | "MCT1 and MCT2 immunoreactivities were respectively found baso-laterally in S1 and thick ascending limbs" — **this**, not Becker, is the S1 source |


---

## What this changes in the two documents

**Attributions that are wrong as written (fix before the paper quotes them).**
1. fig6, Agt row: "Agt mRNA is S2/S3 (Matsusaka 2012)". Matsusaka's ISH puts Agt **mRNA in S3 / the outer stripe
   only**, with S1 and S2 carrying reabsorbed protein. The "normally highest in S2" part of notable result 3 is
   Chen 2021 microdissection (11.7 / 322.1 / 217.1 TPM), so split the citation.
2. fig6, Cyp24a1 row: the TPM triple "85.9 / 2.4 / 0.0" is **Chen 2021**, not Fuchs 2024. Fuchs supports only
   "proximal tubules of the outer cortex" by RNAscope and never names S1/S2/S3.
3. fig6, OXPHOS row: Tran 2016 shows **no post-ischaemic OXPHOS fall** (its OXPHOS panels are the PGC1α-mouse
   contrast), Funk 2012's **PGC-1α protein rose** and its **cortical ATP rose**, and Arnoux 2025's Fig. 3a shows the
   electron-transport score **not reduced in PTS1-PTS3** after IRI. Cite Funk for respiratory-chain proteins, Tran
   for PGC1α, and move Arnoux from support to counterpoint. This makes "New" more defensible, not less.
4. fig6, Cyp51 row: Naito 2009's gene is **Hmgcr**, never Cyp51, so the "phase-dependent conflict" compares two
   different genes.
5. fig6, Lazzeri number: 13.2% is the share of **all tubular cells in S/G2/M** (against 47.1% Ki-67+), not the share
   of Ki-67+ cells that divide.
6. fig6, ZIP row: Cai 2025 never mentions **ZIP14**, and its S3 panel is an SI reanalysis of published data;
   Fujishiro 2012 is the actual ISH source for ZIP8/ZIP14 in the outer stripe.
7. fig6, MCT1 row: Balzer 2022 never mentions Slc16a1 (the 3x is ours from Supplementary Data 2, n = 2 per injured
   group), and Yu 2025's knockout is **MCT4/Slc16a3**, with MCT1 supported only pharmacologically. The
   "basolateral and S1-high" statement should cite Yanase 2008.
8. fig5/fig5lit, SGLT row: "mouse SGLT1 is already high in S2" is **not** in Vrhovac 2015.
9. fig5, Edwards row: "NEP polarity flips" should be "the within-PT axial direction of NEP **activity** differs
   between rat and rabbit"; it is microdissected enzyme activity, and neither species is ours.
10. fig5lit, Odc1 row: the two Koibuchi 1993 papers are being merged — the Histochemistry one is mRNA and says
    "outer portion" (never "pars convoluta"), the J Endocrinol one is the IHC that says pars convoluta. And
    Levillain's full text reports that **rat ODC mRNA is highest in the outer stripe**, opposite to its activity,
    which matters because our readout is mRNA.
11. fig5lit, Gpt row: "ALT ... rises into the straight segment" is our inference, not Chan 1979's abstract; its
    6-fold PST rise is a postnatal-development result.

**Labels to drop.** "Abstract only" is wrong for Witzgall 1994, Funk 2012, Fujii 2019, Naito 2009 and
Kontogiannis 1998 (all five full texts were read here), and Ransick 2019 can be upgraded from ABSTRACT in fig5lit.

**Data-description fixes.**
- Liu 2017 / GSE98622 has **no 48-h sham** (shams are at 4 h, 24 h and 12 mo), so "48 h against sham" must name the
  reference; and every fold change except Slc7a13 and Havcr1 is our reanalysis, not the paper's result.
- Chen 2023's **male samples are the Chen 2021 samples**, so "both sexes" is one male set plus one female set, not
  two independent mouse datasets.
- Cheval 2012 never states sex or mouse strain; it cannot be described as male. Its libraries pool 8 mice or 9
  humans, and the human tissue is tumour-nephrectomy cortex at 59 +/- 10 y against 8-10-week mice.
- Chen 2021's **PTS2 is the straight tubule in cortical medullary rays**, not late PCT.
- Lake 2023: fig5 and fig5lit give different donor counts for the same object (17/12/23 against 13/10/17), and in
  the paper's single-cell data PT-S1 and PT-S2 "could not be accurately resolved". Lake's own segment labels were
  checked against Chen 2021, so Lake-confirms-Chen is mildly circular.
- Xiong 2023 names only **Hao2** in its text; Sh3gl2, Lrig1, Ggact and Odc1 come from its supplementary tables.
- Limbutara 2020 (KTEA) detects far fewer proteins in S1 (2,734) than S2 (4,750), so an S1 zero may be
  non-detection.

**Two citations we are missing.** The unsourced "human SMIT2 (SLC5A11) matches an independent donor" has a
published source (Schlosser 2020, PMID 31959995, human nuclei). And the only mammalian SMIT2 protein localization
(Lahjouji 2007, rabbit) places it in the **proximal convoluted** tubule, which matches our mouse limb — so the
mouse-early half of that reversal has precedent and only the human-late half is new.

**Nomenclature.** HGNC lists MCT4 as an alias of **both** SLC16A3 and SLC16A4, and MCT5 only for SLC16A4. fig6's
correction is right by the standard (Halestrap) numbering, but "MCT4" in older papers is genuinely ambiguous, so
keep the gene symbol beside the protein name.

**What survives untouched.** The handoff-direction evidence (Li 2022, Wu 2024, Xuanyuan 2025, Matsuzaki 2008,
Schneider 2007), the S3-proliferation precedents (Witzgall 1994 in rat, Humphreys 2008/2011 with the caveat that
neither compares S1/S2 against S3), the endoreplication caveat (De Chiara 2022), Fujii 2019's zone-level ATP loss,
Egli-Spichtig 2018's folic-acid contradiction, Sakaguchi 2024's SLC17A3 values, Cheval 2012 and the Knepper
microdissection series as resources, and every "no axial literature exists" statement in fig5lit that was checkable.

## Per-paper records


### Group A — Fig. 5 — transporters, gene localization and the shared data resources

*Five of these (Shen 2015, Ljubojevic 2007, Karbach 2000, Hosoyamada 2004, plus Vrhovac, Breljak, Edwards and Motohashi) were re-attempted in Group E below. Where Group E obtained full text, its record supersedes the abstract-based one here.*

#### Vrhovac 2015 — SGLT1/SGLT2 localization in human kidney (PMID 25304002)
1. **Access:** INACCESSIBLE — abstract only. Tried: `link.springer.com/content/pdf/10.1007/s00424-014-1619-7.pdf` (200 but text/html stub), Springer article page (303 → IdP auth), Semantic Scholar `openAccessPdf` (CLOSED), OpenAlex (`oa_status: closed`, no repository full text), Europe PMC (no PMCID; "Subscription required"), bib.irb.hr/720071 → CroRIS waiting page, FULIR (fulir.irb.hr) title/keyword search (no record), CORE search (403). Full text **N**.
2. **Our claim:** fig5 — "human SGLT2 protein covers S1 and S2, SGLT1 is the S3 transporter, Vrhovac 2015". fig5lit line 85 — "Also protein: hSGLT2 occupies S1 *and* S2 brush border while hSGLT1 is the S3 transporter (Vrhovac 2015, PMID 25304002, ABSTRACT); mouse SGLT1 is already high in S2."
3. **Verdict: PARTLY (abstract).** The human statement is verbatim correct:
   - (Abstract) "hSGLT2 and SGLT1 were localized to the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively."
   - (Abstract) "Novel affinity-purified antibodies against human SGLT1 (hSGLT1) and SGLT2 (hSGLT2) were used" — antibody validation details are in the inaccessible full text.
   - (Abstract) "the expression of both hSGLTs was sex-independent" — useful for us: the human arm of this handoff is not an androgen story.
   **Correction:** the appended clause "mouse SGLT1 is already high in S2" is not supported by anything available from this paper. The abstract says the opposite direction of similarity: "The renal locations of both transporters largely resembled those in rats and mice", with the only stated rodent difference being that "the renal expression of hSGLT1 was absent in thick ascending limb of Henle (TALH) and macula densa" in human, unlike rodents. A mouse-S2 SGLT1 claim needs a rodent source (or should be labelled as our own mouse microdissection crossing at 0.50/0.53).
4. **Axial wording we can borrow:** "the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively"; "the renal locations of both transporters largely resembled those in rats and mice".

#### Nagamori 2016 — AGT1/SLC7A13 is the second rBAT partner (PMID 26739563)
1. **Access:** PMC article HTML, `pmc.ncbi.nlm.nih.gov/articles/PMC4725474/` (pre-downloaded `PMC4725474.html`, converted). Full text **Y**.
2. **Our claim:** fig5 — "part of the mouse-early position is androgen: AGT1 is male-only in mouse and the female crossing is later". fig5lit — "mouse AGT1 is S3-apical and male-only (Nagamori 2016, PMID 26739563, FULL TEXT — that paper has no human data)". fig6 — "AGT1 is apical S3, male-only; b0,+AT is S1-high (Nagamori 2016)".
3. **Verdict: SUPPORTED** for AGT1 (own data), **PARTLY** for b0,+AT (cited, not measured).
   - (Results, Localization in Mouse Kidney) "revealed strong AGT1 immunoreactivity on the renal tubules in the outer stripe of the outer medulla and medullary ray"; "The AGT1 protein was mainly localized on the apical membrane of the tubules".
   - (Results) "the antibodies specifically recognized a band for AGT1 in the male kidney but not in the female (Fig. 1A). mRNA for AGT1 was also detected only in the male kidney".
   - (Discussion) "It indicates that AGT1 is present in the S3 segment of the proximal tubules because the strong rBAT staining is a marker for the S3 segment (11)."
   **Corrections / nuances:** (i) the S3 call is *inferred* from the OSOM + medullary-ray position and the rBAT marker, not from direct S2/S3 discrimination; (ii) two independent antibodies plus in-situ hybridization agree, and a previous antiserum had put AGT1 basolaterally in proximal straight and distal convoluted tubules — the paper argues that was off-target; (iii) **b0,+AT S1-high and the rBAT S3→S1 decline are prior literature inside this paper** ("rBAT is the most abundant in the S3 segment of proximal tubules, and its expression declines toward the S1 segment (11, 12). In contrast, the expression of b0,+AT is highest in the S1 segment"), while their own sex result is "rBAT protein, as well as AGT1, is more abundant in male than female whereas the amount of b0,+AT has no biological sex difference"; (iv) IHC used 8-week-old C57BL/6J mice; (v) confirmed: **no human data** anywhere (human appears only in citations).
4. **Axial wording we can borrow:** "rBAT is the most abundant in the S3 segment of proximal tubules, and its expression declines toward the S1 segment"; "the expression of b0,+AT is highest in the S1 segment and decreases toward the S3 segment".

#### Breljak 2016 — NaDC3 and OAT1-3 along the human nephron (PMID 27053689)
1. **Access:** INACCESSIBLE — abstract only (Europe PMC `resultType=core`). Tried: `physiology.org/doi/pdf/10.1152/ajprenal.00113.2016` (403 to curl), Semantic Scholar and OpenAlex (closed, no repository copy), bib.irb.hr/825186 → CroRIS waiting page, Göttingen DSpace-CRIS records (503 / 401), CORE search (403). Full text **N**.
2. **Our claim:** fig5 — "OAT2 protein is S1/S2-dominant in human, so OCTN2 → OAT2 is unresolved". fig5lit — "human OAT2 protein sits in S1/S2-dominant cortical PT with the outer-stripe S3 unstained (Breljak 2016, PMID 27053689, ABSTRACT...)".
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "OAT1-OAT3-related immunostaining in human kidneys was detected only in the BLM of cortical proximal tubules; all three OATs were stained more intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment in the outer stripe remained unstained."
   - (Abstract) "The specificity of antibodies was validated with human embryonic kidney-293 cells stably transfected with functional OATs."
   - (Abstract) "sex differences in their expression could not be detected."
   **Nuance to add when quoting:** the result is for OAT1, OAT2 and OAT3 jointly and is **basolateral-membrane** staining; the S3 statement is split (weak in medullary rays, absent in the outer stripe), which is a stronger "early" statement than "S1/S2-dominant" alone.
4. **Axial wording we can borrow:** "stained more intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment in the outer stripe remained unstained"; "the BLM of cortical proximal tubules".

#### Shen 2015 — OAT2/SLC22A7 creatinine transporter (PMID 25904762)
1. **Access:** INACCESSIBLE — abstract only. Tried `dmd.aspetjournals.org/content/dmd/43/7/984.full.pdf` (403), Semantic Scholar (rate-limited then closed), OpenAlex (closed). Full text **N**.
2. **Our claim:** fig5lit — "basolateral+apical in human, apical only in rat — PMID 25904762, ABSTRACT".
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "OAT2 protein was localized to both basolateral and apical membranes of human and cynomolgus monkey renal proximal tubules, but appeared only on the apical membrane of rat proximal tubules."
   **Nuance:** monkey behaves like human (worth one clause, since it makes the rat the outlier); the IHC carries **no segment resolution**, so it cannot say where along the PT.
4. **Axial wording we can borrow:** "species difference in the transporter handling of creatinine" (the paper's own framing is species, not axis).

#### Ljubojević 2007 — OAT2 in rat and mouse kidney, sex hormones (PMID 16885152)
1. **Access:** INACCESSIBLE — abstract only (physiology.org 403; no OA copy). Full text **N**.
2. **Our claim:** fig5lit — "Rodent OAT2 is S3 brush border, female > male and strongly androgen-inhibited (PMID 16885152, ABSTRACT), so the male mouse under-reports OAT2 at S3."
3. **Verdict: PARTLY (abstract)** — right for rat, thinner for mouse.
   - (Abstract) "By IC, the antibody exclusively stained brush border (BB) of the proximal tubule S3 segment (S3) in the OS and medullary rays (F > M)."
   - (Abstract) "strong downregulation by testosterone and weak upregulation by estradiol and progesterone treatment"; "gender differences (F > M) that appear in puberty".
   - (Abstract) "the protein in adult M and F mice ... In mice, the expression pattern largely resembled that in rats."
   **Correction:** mouse evidence is protein in adult males and females only, reported as "largely resembled" the rat; there is no mouse segment IHC series or hormone manipulation. Write "rat (mouse protein similar)" rather than "rodent". Also note the baseline is low and sex-independent before puberty.
4. **Axial wording we can borrow:** "brush border (BB) of the proximal tubule S3 segment (S3) in the OS and medullary rays"; "expression of OAT2 mRNA was predominant in the outer stripe (OS) tissue".

#### Edwards 1999 — neutral endopeptidase along the rat and rabbit nephron (PMID 10352425)
1. **Access:** INACCESSIBLE — abstract only (Karger blocks curl and WebFetch; Semantic Scholar/OpenAlex closed). Full text **N**.
2. **Our claim:** fig5 — "Partial precedents: NEP polarity flips between rat and rabbit microdissection (Edwards 1999, abstract)", supporting notable result 1 (neprilysin early in human, late in mouse).
3. **Verdict: PARTLY (abstract)** — the species contrast is real; the word "polarity" is wrong.
   - (Abstract) "In the rat, only the proximal straight tubule and glomerulus had measurable NEP activity of 86 +/- 11.3 pmol/min/mm tubule length and 5.8 +/- 1.5 pmol/min/glomerulus, respectively."
   - (Abstract) "In the rabbit, significant activity was observed in both the proximal convoluted tubule (70.8 +/- 7.2 pmol/min/mm) and proximal straight tubule (29.6 +/- 2.3 pmol/min/mm)".
   - (Abstract) "The results demonstrate that in both the rat and rabbit NEP is restricted to the glomerulus and proximal tubule."
   **Correction:** this is microdissected **enzyme activity** (not protein polarity, not mRNA), and the flip is in the **within-PT axial direction**: rat = straight segment only (late), rabbit = convoluted > straight (early). Species are rat and rabbit — the paper contains **no human and no mouse** data, so it is a precedent for a species-dependent NEP axis, not direct support for "early in human, late in mouse". Suggested wording: "the direction of the NEP gradient within the proximal tubule already differs between rodent species: activity is confined to the straight tubule in rat but is highest in the convoluted tubule in rabbit (Edwards 1999, abstract)".
4. **Axial wording we can borrow:** "only the proximal straight tubule and glomerulus had measurable NEP activity"; "significant activity ... in both the proximal convoluted tubule ... and proximal straight tubule".

#### McEvoy 2022 — single-cell map of healthy human kidney (PMID 36496458)
1. **Access:** Europe PMC JATS XML (`PMC9741629`), converted. Full text **Y**.
2. **Our claim:** fig5 — "DCXR is S3-high in human single-cell data (McEvoy 2022) and S1-high in rat". Ledger — Fig. 1e "PT2 shows increased expression of S3-abundant genes (DCXR, AGXT, SLC22A7, SLC7A13)".
3. **Verdict: PARTLY** — the quote is exact, but what it establishes is cluster enrichment of literature-defined S3 markers, not an independent S3 measurement.
   - (Results, Single-cell map) "PT segment-specific separation is evident; PT1, 4, and 6 are enriched for PT segment 1 (S1) marker SLC5A2 and S1/2-abundant genes (SLC7A7, ANK2, SLC4A4, SLC6A19, SLC22A8), while PT2 shows increased expression of S3-abundant genes (DCXR, AGXT, SLC22A7, SLC7A13) (Fig. 1e)9,10."
   - (Fig. 1 legend) "d UMAP plot of compartment-specific analysis of 20772 proximal tubular cells, comprising 6 clusters. e Heat map showing the expression levels of cluster marker genes."
   - (Methods, Human specimens) "20 living donor kidney samples (10 male donors and 10 female donors) were processed for sequencing. One sample was subsequently removed due to quality control reasons".
   **Corrections:** (i) "S3-abundant" is imported from their refs 9/10 — DCXR's S3 identity is assumed by this paper, so citing it for "DCXR is S3-high in human" is circular unless we say "consistent with the S3 marker set used by McEvoy 2022"; (ii) scRNAseq of 19 living donors, **both sexes** (that paper's main theme is sex-biased PT expression), so it is not a male-only comparator; (iii) **no S2 cluster** — PT1/4/6 are S1 and "S1/2", PT2 is S3, PT3 is dissociation-stress, PT5 is VCAM1+ injury-like; (iv) the rat S1-high claim is ours, not theirs — stated correctly in the doc.
4. **Axial wording we can borrow:** "PT segment-specific separation is evident"; "PT1, 4, and 6 are enriched for PT segment 1 (S1) marker SLC5A2 and S1/2-abundant genes ... while PT2 shows increased expression of S3-abundant genes".

#### Sakaguchi 2024 — three cell populations for urate excretion in human kidney (PMID 38166558)
1. **Access:** Europe PMC JATS XML (`PMC10763458`), converted. Full text **Y**.
2. **Our claim:** fig5 — "human SLC17A3 falls S1 → S3 (Sakaguchi 2024)" (supporting "the urate exporter NPT4 sits early in the human PT"). Ledger — "Table 1: human SLC17A3 PT_S1 6.02, PT_S2 3.82, PT_S3 3.56 ... 3 adult male donors".
3. **Verdict: SUPPORTED** (values and donors both check out).
   - (Table 1) "SLC17A3 NPT4 Secretion AE Facilitate (Δψ-driven) Organic anions (e.g., PAH, urate) **6.02 3.82 3.56 2.18** Unknown" — columns are Av. (PT_S1), Av. (PT_S2), Av. (PT_S3), Av. (LOH (DL)).
   - (Methods, snRNA-seq data sets) "the data sets must be from kidneys of adult human males"; "the data sets used in this study consisted of three males and a total of 10,080 cells".
   - (Methods, Annotation of segments) "As there was no report of specific marker genes for the S2 segment, the cluster which did not fall into S1 or S3 segments in three PT clusters was assigned as the S2 segment."
   **Resource caveats for our wording:** reanalysis of **public** datasets (GSE118184, GSE131882), not new tissue; mRNA only — **no protein** is shown (the apical-efflux assignment in Table 1 is cited from functional papers); three discrete PT clusters with S2 defined by exclusion; values are log-normalised average expression, so "6.02 → 3.56" is not a fold change.
4. **Axial wording we can borrow:** "across the cell populations ... along the three PT segments (S1–S3) and DL"; "PT clusters are anatomically divided into three segments".

#### Cheval 2012 — Of Mice and Men: divergence of gene expression patterns in kidney (PMID 23056504)
1. **Access:** Europe PMC JATS XML (`PMC3463552`), converted. Full text **Y**.
2. **Our claim (resource):** fig5lit table — "Cheval 2012 PLoS One, PMID 23056504 | DATA | human **and** mouse microdissected S1 vs S3 SAGE (one pooled library per segment)".
3. **Verdict: SUPPORTED (resource)** as described, with two things to add to our description.
   - (Results, Transcriptome Database) "We analyzed SAGE libraries previously generated in our laboratory from glomeruli (Glom), initial and terminal portions of the proximal tubule (S1 and S3), medullary and cortical thick ascending limbs ... dissected from human [8] and mouse kidneys [9]."
   - (Results) "biological variability was minimized by pooling samples from several individuals (8 mice and 9 humans) in each library"; "The number of sequenced tags varied from 43,000 to 100,000 tags per library ... tag abundances were normalized to 10,000 tags in each library."
   - (Discussion) "human samples emanated from healthy kidney fragments of donors who had undergone surgery for removal of kidney tumors [8]"; "the patients analyzed in this study were much older (59±10 years), in relative terms, than the mice (8–10 weeks)".
   **Flags:** (i) no S2 — only S1 ("initial") and S3 ("terminal") portions of the PT, so Cheval cannot arbitrate S2 positions; (ii) **sex is never stated for either species, and the mouse strain is not given** — our documents should not describe these libraries as male; (iii) human tissue is tumour-nephrectomy cortex from 59 ± 10-year-olds versus 8–10-week-old mice, i.e. an age/comorbidity mismatch the authors themselves raise as a caveat; (iv) one library per structure means no within-group variance and no replicate statistics (the authors note SAGE "precludes multiple testing and assessment of assay variability"); the only reproducibility estimate is two independent mouse OMCD libraries (slope 1.019, r = 0.742).
4. **Axial wording we can borrow:** "initial and terminal portions of the proximal tubule (S1 and S3)"; "renal functions result from the spatial integration of elementary processes originating in the glomerulus and the successive segments constituting the nephron".

#### Hosoyamada 2004 — URAT1 function and localization in mouse kidney (PMID 14747372)
1. **Access:** INACCESSIBLE — abstract only (journals.lww.com not reachable by curl/WebFetch; Semantic Scholar and OpenAlex closed). Full text **N**.
2. **Our claim:** fig5lit — "Mouse URAT1 is brush-border and male > female (PMID 14747372, ABSTRACT)".
3. **Verdict: PARTLY (abstract)** — correct as stated, but it carries no segment information.
   - (Abstract) "the brush border or intracellular membrane of proximal tubules was stained by an affinity-purified antibody that recognized mouse URAT1 (mURAT1) expressed on Xenopus oocyte."
   - (Abstract) "The 70-kD protein was N-glycosylated and was identified as a Triton X-100 insoluble brush border membrane protein."
   - (Abstract) "RST mRNA and protein levels were higher in male kidneys than female."
   **Correction:** "proximal tubules" only — **no S1/S2/S3 assignment**, so this paper cannot support (or contradict) a position for the SLC22A23 → URAT1 crossing. Keep the doc's existing honest note that no human URAT1 segment-level protein data were found.
4. **Axial wording we can borrow:** "brush border or intracellular membrane of proximal tubules" (no axial statement available).

#### Sato 2011 — rat Urat1 (Slc22a12) (PMID 21074513)
1. **Access:** INACCESSIBLE — abstract only. Semantic Scholar/OpenAlex report a BRONZE PDF at `sciencedirect.com/science/article/pii/S0005273610003755/pdf`, but ScienceDirect answers 403 (tdm-reservation). Full text **N**.
2. **Our claim:** fig5lit — "rat URAT1 apical (PMID 21074513, ABSTRACT)".
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "rUrat1 was localized at the apical membrane in proximal tubular epithelial cells in rat kidney."
   - (Abstract) "Uric acid transport by rUrat1 was chloride (Cl-) susceptible with a Km of 1773μM" (brush-border membrane vesicles plus IHC).
   **Nuance:** no segment resolution; rat, not mouse.
4. **Axial wording we can borrow:** "localized at the apical membrane in proximal tubular epithelial cells" (no axial statement available).

#### Karbach 2000 — OCT1 and OCT2 localization in rat kidney (PMID 10997918)
1. **Access:** INACCESSIBLE — abstract only (physiology.org 403; the Europe PMC "free after 12 months" link points to the dead `intl-ajprenal.physiology.org` host). Full text **N**.
2. **Our claim:** fig5lit — "Rodent protein agrees that OCT2 is late (rOCT2 mainly S2/S3, PMID 10997918, ABSTRACT)".
3. **Verdict: SUPPORTED (abstract),** with a useful extra detail and one caution.
   - (Abstract) "By in situ hybridization, significant amounts of both rOCT1 and rOCT2 mRNA were detected in S1, S2, and S3 segments of proximal tubules."
   - (Abstract) "expression of the rOCT1 protein was mainly observed in S1 and S2 segments of proximal tubules, with lower expression levels in the S3 segments. At variance, rOCT2 protein was mainly expressed in the S2 and S3 segments."
   - (Abstract) "Both transporters were localized to the basolateral cell membrane."
   **Notes:** the axial information is **protein only** — the mRNA is present in all three segments, which is worth saying because our own measurement is mRNA. **No sex comparison is reported in this paper**; the androgen half of our sentence rests entirely on Urakami 2000. Species is rat.
4. **Axial wording we can borrow:** "rOCT2 protein was mainly expressed in the S2 and S3 segments"; "mainly observed in S1 and S2 segments of proximal tubules, with lower expression levels in the S3 segments".

#### Urakami 2000 — hormonal regulation of OCT2 in rat kidney (PMID 10812069)
1. **Access:** INACCESSIBLE — abstract only (Wiley/Elsevier closed; Semantic Scholar and OpenAlex closed). Full text **N**.
2. **Our claim:** fig5lit — "and is testosterone-driven (PMID 10812069, ABSTRACT), so the mouse-early crossing is androgen-amplified".
3. **Verdict: SUPPORTED (abstract)** for the hormone effect; it cannot speak to position.
   - (Abstract) "Treatment of male and female rats with testosterone significantly increased the expression levels of rOCT2 mRNA and protein in the kidney, whereas estradiol treatment moderately decreased the expression levels of rOCT2."
   - (Abstract) "There was no regulation of renal rOCT1 mRNA expression by testosterone or estradiol."
   - (Abstract) "renal rOCT2 expression in males is much higher than that in females".
   **Corrections:** species is **rat**, not mouse; the readouts are whole-kidney mRNA/protein and renal-slice TEA accumulation, so there is **no segment-level information** — the paper supports "OCT2 level is androgen-driven in rodents", not "the crossing position is androgen-shifted". Our clause should separate the two.
4. **Axial wording we can borrow:** none available (no axial measurement); only "renal rOCT2 expression in males is much higher than that in females".

#### Motohashi 2013 — OCT, OAT and MATE protein localization in human kidney (PMID 23630107)
1. **Access:** INACCESSIBLE — abstract only (`jpharmsci.org/article/S0022354915309333/pdf` 403; Semantic Scholar/OpenAlex closed). Full text **N**.
2. **Our claim:** fig5lit — "Human OCT2 protein is PT-wide with no segment resolution (PMID 23630107, ABSTRACT)".
3. **Verdict: SUPPORTED (abstract),** with one tightening.
   - (Abstract) "the expression and distribution of apical MATE1 and MATE2-K, and basolateral OAT1, OAT3, and OCT2 were compared using serial sections of human kidney cortex."
   - (Abstract) "Levels of OAT, OCT2, and MATE mRNA in the proximal tubules were greatly higher compared with glomerulus. The results quantitatively indicated that these transporters were localized to proximal tubules in the renal cortex."
   - (Abstract) "MATE1 and MATE2-K protein were detected in proximal epithelial cells in which OCT2 protein was expressed at the basolateral membranes."
   **Correction:** only **cortex** was examined (laser microdissection + serial-section IHC), so the honest statement is "cortical proximal tubule, without S1/S2/S3 resolution" rather than "PT-wide" — the outer-stripe S3 was not sampled, which matters because that is exactly where a late human OCT2 would sit.
4. **Axial wording we can borrow:** "serial sections of human kidney cortex"; "localized to proximal tubules in the renal cortex".

#### Zhang 2024 — fructose transporters in the rat nephron (PMID 38474316)
1. **Access:** Europe PMC JATS XML (`PMC10931557`), converted. Full text **Y**.
2. **Our claim:** fig5lit — "Rat SGLT5 rises progressively S1→S3 (PMID 38474316, FULL TEXT)".
3. **Verdict: SUPPORTED.**
   - (Abstract) "key fructose transporters SGLT5 and GLUT5 progressively increased from S1 to S3, and both were significantly upregulated in S3 vs. S1/S2 (Slc5a10: 1.9 log2FC, p < 1 × 10−299 ...)".
   - (Results, Transport and Metabolism of Fructose) "Both SGLT2 (Slc5a2) and GLUT2 (Slc2a2) were upregulated in PT.S1 while SGLT5 (Slc5a10) was downregulated. ... Finally, SGLT5 (Slc5a10) and GLUT5 (Slc2a5) were both upregulated in PT.S3".
   - (Results) "all four fructose transporters were expressed in both S2 and S3, and fructose transporter mRNA and protein expression tightly paralleled one another" (Fig. 5: microdissected-segment proteomics in "million copies per cell").
   **Method / sex, for accuracy:** this is a **reanalysis of public data**, not ISH or new IHC — rat scRNAseq GSE137869 ("three pools of male ... and three pools of female ... Sprague Dawley rats", merged), rat microdissected bulk RNAseq PRJNA24440, published rat PT-subsegment proteomics, and human KPMP snRNAseq for the human comparison. So the rat statement is mixed-sex transcriptomics with supporting proteomics, and the same direction is reported in human (Appendix Fig. A8).
   **Bonus caution the paper hands us:** (Abstract) "PT-S2 express a mixture of 18% ± 9% S1, 58% ± 8% S2, and 19% ± 5% S3 transcripts" and "anatomical labels applied to scRNAseq clusters may be misleading" — a citable argument for a continuous axis over discrete segment labels.
4. **Axial wording we can borrow:** "progressively increased from S1 to S3"; "a transcriptome cluster may obscure the heterogeneity of the cell population within a nephron segment".

#### Fukuzawa 2013 — SGLT5 reabsorbs fructose in the kidney (PMID 23451068)
1. **Access:** Europe PMC JATS XML (`PMC3581502`), converted. Full text **Y**.
2. **Our claim:** fig5lit — "no human SGLT5 localization exists (PMID 23451068, FULL TEXT, is BBM function in mouse)".
3. **Verdict: SUPPORTED.**
   - (Discussion) "Although we have no direct evidence showing the presence of this transporter in proximal tubules because a specific SGLT5 antibody is not available, we looked at fructose uptake in a BBMV assay."
   - (Abstract) "Sodium-dependent fructose uptake disappeared in renal brush border membrane vesicles from SGLT5-deficient mice, and the increased urinary fructose in SGLT5-deficient mice indicated that SGLT5 was the major fructose reabsorption transporter in the kidney."
   - (Discussion) "In situ hybridization shows that SGLT2, which exists in segments S1 and S2 [39], and SGLT5 are co-expressed in epithelial tissue of the glomerulus and distal and proximal tubules [40]" — a **cited** ISH result, not theirs, and with no S1/S2/S3 resolution.
   **Resource details:** Slc5a10-knockout on C57BL/6 background, male mice (metabolic studies at 17–21 weeks); human appears only as cited in-vitro work on human SGLT5 ("human SGLT5 is exclusively expressed in the kidney ... it mediates the transport of fructose and mannose [25]"). So: no protein localization in any species, no segment assignment, no human tissue.
4. **Axial wording we can borrow:** "renal brush border membrane vesicles"; "we have no direct evidence showing the presence of this transporter in proximal tubules because a specific SGLT5 antibody is not available".

### Group B — Fig. 5 — gene-level precedents and the microdissection / atlas resources

#### Levillain 1998 — Ornithine decarboxylase along the mouse and rat nephron (PMID 9841492)
1. **Access:** Europe PMC metadata showed subscription only. journals.physiology.org returned 403 to curl (doi/pdf and the legacy ajprenal URL), and WebFetch also failed. Semantic Scholar has `openAccessPdf` CLOSED. A WebSearch on the exact title found only index and mirror pages (unibo mirror 403). Judged from the PubMed abstract (NCBI efetch). **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim (fig5lit, Odc1 row):** "Rodent literature agrees with the mouse direction: activity PCT > PST (PMID 9841492, ABSTRACT)".
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "In Triton X-100-treated tubules from control mice, ODC was exclusively found in proximal tubules (PCT > PST)."
   - (Abstract) "the highest ODC activity was found in rat and mouse PCT, a segment devoid of arginase."
   - Notes: the readout is enzyme activity (L-[1-14C]Orn decarboxylation) in microdissected viable tubules, not mRNA. Both species were used, and the PCT is highest in both. The abstract does not state sex or strain. A later Pflügers Arch 2000 paper by the same group used female mice; it is not the cited one.
4. **Axial wording:**
   - "In control mice, Orn decarboxylation occurred mainly in the proximal convoluted tubule (PCT)."
   - "In rats, Orn decarboxylation also occurred predominantly in the proximal tubule."

#### Koibuchi 1993 — Heterogeneous expression of ODC gene in the proximal tubule of the mouse kidney following testosterone (PMID 8307774)
1. **Access:** Springer link.springer.com returned a JS client challenge to curl and redirected WebFetch to the login IdP. Semantic Scholar CLOSED; the Crossref TDM link also hit the challenge. Judged from the PubMed abstract. **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim (fig5lit, Odc1 row):** "testosterone-induced mRNA stronger in pars convoluta (PMIDs 8307774, 8429279, ABSTRACT)".
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "Testosterone (n = 5) or vehicle (n = 5) was subcutaneously injected (1 mg/animal) into male BALB/c mice (8 weeks in age) 14 h before sacrifice."
   - (Abstract) "The hybridization signal was greater in the outer portion of the proximal tubule than in the inner portion."
   - **Correction:** This is the mRNA (in situ, 35S riboprobe) paper. It says "outer portion" vs "inner portion" of the cortical PT, not "pars convoluta". Write "stronger in the outer cortical PT" unless the full text equates this with the pars convoluta. The design uses intact (not castrated) males and a single supraphysiological testosterone dose, read out 14 h later. It describes an induced pattern, not basal androgen maintenance.
4. **Axial wording:**
   - "A marked increase in ODC mRNA was exclusively detected in the proximal tubule of the renal cortex in the testosterone-treated animals."

#### Koibuchi 1993 — Induction of ODC immunoreactivity in the male mouse kidney following testosterone: an axial heterogeneity in the PT (PMID 8429279)
1. **Access:** Crossref lists the joe.bioscientifica.com HTML and the journals.bioscientifica.com PDF. Both returned 403 to curl (Cloudflare) and to WebFetch. Semantic Scholar CLOSED. Judged from the PubMed abstract. **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim:** as above ("testosterone-induced mRNA stronger in pars convoluta").
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "In the pars recta, cells were stained weakly and homogeneously, while in the pars convoluta, the luminal surface of the cells showed stronger immunoreactivity."
   - (Abstract) "few cells that were obviously immunoreactive to ODC were observed in the control animals"
   - (Abstract) "Male BALB C mice at 8 weeks of age were used throughout. Fourteen hours before death, they received a subcutaneous injection of testosterone (1 mg/animal) or solvent"
   - **Correction:** This is **protein (immunocytochemistry)**, not mRNA. Our sentence merges the two papers. Suggested text: "testosterone-induced ODC mRNA was higher in the outer cortical PT (ISH, 8307774), and induced protein stained the pars convoluta more strongly than the pars recta (IHC, 8429279)". Both studies used intact adult males with a single 1 mg testosterone dose, read out at 14 h.
4. **Axial wording:**
   - Title phrase: "an axial heterogeneity in the proximal tubule".
   - "many granule-like particles that were strongly ODC immunoreactive were observed inside the lumen of the pars convoluta."

#### Koibuchi 1990 — Localization by IHC of renal ODC in the mouse with and without testosterone (PMID 2083533)
1. **Access:** J-STAGE PDF (https://www.jstage.jst.go.jp/article/endocrj1954/37/4/37_4_555/_pdf). The text layer was garbled, so all 7 pages were rendered with pdftoppm and read. **Full text: Y.**
2. **Our claim (fig5lit):** "though an earlier IHC put most signal in pars recta (PMID 2083533, ABSTRACT)".
3. **Verdict: SUPPORTED, with caveat.** The status can be upgraded from ABSTRACT to full text.
   - (Methods, Animals) "Male BALB C mice (8 weeks of age)". Group 1 received 1 mg testosterone sc 15 h before sacrifice; group 2 received the solvent.
   - (Results) "it was not established whether ODC-immunoreactive cells mainly localized in the proximal tubule or in the distal tubule."
   - (Results) "More immunoreactive cell bodies were observed in the pars recta."
   - (Abstract) "no significant immunohistochemical difference was observed either in immunoreactivity or its localization."
   - Notes: the model is simply intact males with and without testosterone; there is no castrated or knockout model. The pars recta assignment is by distribution pattern and the authors call it tentative ("seemed to be confined mainly to the pars contorta and the pars recta"). It counts immunoreactive *cells*, which are scattered among non-reactive cells, rather than staining intensity. Testosterone raised activity about 6-fold (Fig. 1) without changing the IHC pattern.
4. **Axial wording:**
   - (Results) "More ODC-immunoreactive cell bodies were observed in the inner portion of the cortex."
   - (Discussion) "ODC-immunoreactive cells are localized mainly in the pars contorta and the pars recta of the proximal tubule."

#### Kim 2018 — Proximal tubular α7 nAChR attenuates ischemic AKI (PMID 29674665)
1. **Access:** pre-converted PMC text (txt/PMC5938048.txt; Europe PMC JATS). **Full text: Y.**
2. **Our claim (fig5lit, Chrna4 row):** "No axial literature in either species; renal nicotinic work is all α7 (e.g. PMID 29674665, ABSTRACT)".
3. **Verdict: SUPPORTED.** The full text has no mention of α4/Chrna4, β2 or any other subunit; all work is α7 (agonist DMBX-A/GTS-21, antagonist MLA).
   - (Results) "In sham and IR mice, α7nAChR was expressed in proximal tubular cells, but not in distal tubular cells (Fig. 4a)."
   - (Methods) "Male C57BL/6 mice (7 weeks old)". In vitro work used HK-2 cells.
   - The PT is identified only as AQP1+ (vs calbindin+ distal). There is no S1/S2/S3 or convoluted/straight resolution. This paper alone cannot back the general statement "all renal nicotinic work is α7"; it is one example.
4. **Axial / injury wording:**
   - "the proximal tubule contains highly metabolically active nephron segments" (Discussion; generic, no axial data).

#### Brie 2024 — GH axis and sexually dimorphic liver signatures (PMID 38370444)
1. **Access:** pre-converted PMC text (txt/PMC10872697.txt). **Full text: Y.**
2. **Our claim (fig5lit, Hao2 row):** "Hao2 is a canonical sex-dimorphic gene elsewhere (PMID 38370444, ABSTRACT, liver)".
3. **Verdict: SUPPORTED.** Direction: **female-predominant**. Species: **mouse liver**. Mechanism studied: the GH pattern via STAT5b, not androgen.
   - (Abstract) "focusing on the male-predominant genes Lcn13, Asns, and Cyp7b1, and the female-predominant genes Hao2, Pgc1a, Hamp2, Cyp2a4, and Cyp2b9."
   - (Introduction) "Of note, HAO2 is a female-predominant gene, and HCC is more prevalent in men".
   - (Results) "In HiGH females, Hao2, Cyp2a4, and Cyp2b9 increased (Fig. 2B and 3B)". Hepatocyte GHR loss raised Hao2 in males.
   - Notes: there are no kidney data. The female-high liver direction matches our female-mouse kidney direction, but the regulator differs (GH pattern in liver; AR in kidney per Xiong 2023).
4. **Wording:** none axial (liver only).

#### Thomasson 2003 — LRIG1 and EGFR in renal cell carcinoma (PMID 14520461)
1. **Access:** pre-converted PMC text (txt/PMC2394322.txt). **Full text: Y.**
2. **Our claim (fig5lit, Lrig1 row):** "No kidney axial literature (only RCC, PMID 14520461, ABSTRACT)".
3. **Verdict: SUPPORTED.** The paper does include normal human cortex, but with no axial information.
   - (Abstract) "Eight matched samples of uninvolved kidney cortex were also evaluated."
   - (Results) "The staining pattern of LRIG1 in the kidney cortex was very similar to that of EGFR with a most prominent staining in the tubuli (Figure 2)"
   - (Results) "EGFR in the kidney cortex displayed a pattern with staining of tubular structures, distal tubules staining slightly more intensely than proximal tubules."
   - Notes: species human; qRT-PCR and IHC. The paper makes no proximal-vs-distal or PT-subsegment statement for LRIG1 itself.
4. **Wording:** none axial.

#### Chan 1979 — Two aminotransferases and D-amino acid oxidase within the nephron of young and adult rats (PMID 39098)
1. **Access:** Europe PMC and Semantic Scholar list a free (BRONZE) SAGE PDF, but journals.sagepub.com (pdf, epdf, landing) returned Cloudflare 403 to curl and WebFetch. The old jhc.org URLs are dead. WebSearch found only index pages (and a zendy viewer, not used). Judged from the PubMed abstract. **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim (fig5lit, Gpt row):** "Rodent precedent: ALT is the mirror image of AST along the proximal tubule and rises into the straight segment (PMID 39098, ABSTRACT)".
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "Alanine aminotransferase (lower by a factor of 100 in absolute terms) was distributed as the mirror image of aspartate aminotransferase within proximal and distal tubules."
   - (Abstract) "Adult aspartate aminotransferase activities (percentage of the highest) ranged from 100 in the distal straight segment to 25 in the late part of the proximal straight segment"
   - (Abstract) "During development alanine aminotransferase increased 6-fold and D-amino acid oxidase, 4.5-fold in proximal straight tubules"
   - **Correction:** The species is rat (adult plus developing); sex is not given in the abstract; the readout is enzyme activity in microdissected nephron parts, with 4 PT parts measured. The abstract does not state that ALT rises from convoluted to straight PT. "Rises into the straight segment" is inferred from the "mirror image" wording and from AST being low in late PST. Write instead: "ALT activity is distributed as the mirror image of AST within the proximal and distal tubule (rat; abstract)". Alternatively, flag the rise as our inference until the PT-part numbers are read.
4. **Axial wording:**
   - "measured in glomeruli, 4 parts of the proximal tubule, 2 parts of the distal tubule and in patches from the thin limb area and the papilla."

#### Davies 1997 — rpr-1, a resiniferatoxin-binding phosphotriesterase-related protein in rat kidney tubules (PMID 9237666)
1. **Access:** The Europe PMC "free PDF" link (onlinelibrary.wiley.com pdfdirect) and the febs.onlinelibrary.wiley.com pdf/full pages returned 403 (Cloudflare) to curl and WebFetch. The ScienceDirect PII page returned 403, the Wiley TDM API returned 400, and Semantic Scholar is CLOSED. Judged from the PubMed abstract. **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim (fig5lit, Pter row):** "Rat PT localization by in situ exists (rpr-1, PMID 9237666, ABSTRACT), no segment resolution".
3. **Verdict: SUPPORTED (abstract)** for rat PT localization by ISH.
   - (Abstract) "Northern blotting and in situ hybridisation show expression primarily in proximal tubules of the kidney, in which rpr-1 distribution correlates with resiniferatoxin-binding activity."
   - The claim "no segment resolution" cannot be verified: the abstract does not mention segments, but the ISH figure might. Phrase it as "no segment-level statement in the abstract". Sex and age are not stated.
4. **Wording:** none axial in the abstract.

#### Lin 2011 — Blocking core fucosylation of TGF-β1 receptors attenuates EMT of renal tubular cells (PMID 21228108)
1. **Access:** physiology.org returned 403; the legacy intl-ajprenal HighWire URL is dead; Semantic Scholar is CLOSED. Judged from the PubMed abstract. **Full text: N. INACCESSIBLE — abstract only.**
2. **Our claim (fig5lit, Fut8 row):** "Renal FUT8 literature is injury/fibrosis only (PMIDs 21228108, 29118904, ABSTRACT), no axial data."
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "in the process of EMT of cultured human renal proximal tubular epithelial (HK-2) cells."
   - (Abstract) "blocking it with Fut8 small interfering RNA greatly reduced the phosphorylation of Smad2/3 protein"
   - Caveat: Wang 2017 (same group) cites this paper as ref. 21 for "the expression of core fucosylation was significantly upregulated in rats with unilateral ureter obstruction (UUO) [21,22]". The full text may therefore contain rat UUO tissue staining that the abstract does not mention. Unverified, but it would be injury tissue, not normal axial data.
4. **Wording:** none axial.

#### Wang 2017 — Inhibiting core fucosylation protects against albumin-induced PTEC injury (PMID 29118904)
1. **Access:** pre-converted PMC HTML text (txt/PMC5666051.txt). **Full text: Y.**
2. **Our claim:** as above.
3. **Verdict: SUPPORTED.** The study is entirely in vitro, with no tissue.
   - (Methods) "Immortalized human renal proximal tubular epithelial cells (HK-2) were purchased from the American Type Culture Collection"
   - (Fig. 1 legend) "BSA activates injury-induced pathways and increases the level of FUT8 expression."
   - Correction (minor): both cited FUT8 papers are HK-2 cell-culture injury models (TGF-β1 EMT; BSA overload). Neither is a tissue fibrosis study. Write "cell-culture injury models only".
4. **Wording:** none axial.

#### Ransick 2019 — Single cell profiling reveals sex, lineage and regional diversity in the mouse kidney (PMID 31689386)
1. **Access:** pre-converted PMC author-manuscript text (txt/PMC6948019.ncbi.txt). **Full text: Y.**
2. **Our claim (fig5lit, marked ABSTRACT):** "existence of the male-biased mouse PT-S2/S3 program".
3. **Verdict: PARTLY.** Sex differences exist and concentrate in S2/S3, but they are **two-way** (male- and female-restricted genes), not a male-biased program.
   - (Results) "Strikingly, while most clusters were evenly represented by all kidney samples (Figure S2B), all proximal tubule clusters separated by sex (Figure 2D)."
   - (Results) "The highest disparity in gene activity between the sexes was in the S3 region."
   - (Results) "Interestingly, S2 and S3 sex-specific expression domains form a tight, non-overlapping boundary at the S2-S3 junction (Figure S2D)."
   - Segment definitions (Results): S1/S2/S3 were defined by markers "that exhibited comparable expression between the sexes: Spp2/Slc5a2 (S1), Slc34a1 (S1–2) and Slc22a7 (S3)". Examples by segment:
     - S3: Cyp7b1 male-only, Slc7a12 female-only.
     - S2: Cyp2e1 male-only, Cyp4a14 female-only.
     - Extending into S2 at lower levels: Acsm3 male, Prlr female.
   - Design (STAR Methods): "Adult kidneys were harvested from 9–10-week C57BL/6J mice (2 males and 2 females)." The Results text says 8–9 weeks. Assays: 10x scRNA-seq on cortex, OM and IM zones, plus RNAscope.
   - Caveat: "male S3 cells were markedly underrepresented in the scRNA-seq profile (~10% of expected...)".
   - Suggested wording: "sexually dimorphic PT programs concentrated in S2 and S3 (male- and female-specific genes; Ransick 2019), with S1 the least dimorphic (Xiong 2023)". Xiong 2023 is the paper that states S1 has little sex difference.
4. **Axial wording (S3 cortex vs outer stripe):** The paper does not split cortical S3 from outer-stripe S3. It places S3 in the outer stripe, while noting the PT was recovered mainly from the cortical zone:
   - "Juxtamedullary nephrons transitioned in the outer medulla from Aqp1+/Slc7a13+ cells in the S3 region in the outer stripe"
   - "proximal tubules, distal convoluted tubules and the cortical collecting duct were recovered from zone 1"
   - On injury: "the S3 region is particularly vulnerable to injury (Park et al., 2004; Neugarten et al., 2018)."

#### Part 2: data resources (group B)

#### Chen 2021 — mRNAs and isoforms across all 14 renal tubule segments of mouse (PMID 33769951)
1. **Access:** pre-converted PMC HTML text (txt/PMC8017530.txt). Values were cross-checked against GSE150338 `tpm.txt.gz` and `knepper_pt.tsv` in the scratchpad. **Full text: Y.**
2. **Our description:** "male mouse microdissected PTS1/PTS2/PTS3 TPM". fig6 quotes "Cyp24a1 … (85.9 / 2.4 / 0.0 TPM; Fuchs 2024)" and "S3 is normally near-silent (0.2 TPM)" for Slc16a1.
3. **Verdict: SUPPORTED; segment definitions and one attribution flagged.**
   - (Methods) "Male, 6- to 8-week-old, C57BL/6 mice (Taconic) were used for renal tubule microdissection." Each sample pooled "Four to eight tubules". Quantified as TPM via RSEM.
   - (Fig. 1 legend) "PTS1, the initial segment of the proximal tubule; PTS2, proximal straight tubule in cortical medullary rays; PTS3, last segment of the proximal straight tubule in the outer stripe"
   - (Results) "Some transcripts are exclusively localized to the PTS1 (turquoise cluster, Figure 2C)". This cluster includes Cyp24a1, Gc, Cyp27b1, Slc5a2 and Slc5a12, which supports "Cyp24a1 normally S1-restricted".
   - Replicates from the GEO column names: PTS1 has 5 samples from 3 mice, PTS2 8 from 5 mice, PTS3 6 from 3 mice. Overall the study had 93 samples from 53 mice.
   - **Flags:**
     - (a) S2 here is the cortical straight tubule in medullary rays, not the late convoluted tubule. "S1" is the initial PT segment. The S1→S2 contrast in this dataset is therefore convoluted vs straight-in-medullary-ray.
     - (b) The quoted means match `knepper_pt.tsv` exactly. A plain mean over the GEO sample TPMs gives Cyp24a1 89.7/2.4/0.0 and Slc16a1 46.5/10.8/0.2, so the source of the averaging should be stated. Cyp24a1 PTS1 varies 25–183 TPM across samples.
     - (c) fig6 line 521 attributes "85.9 / 2.4 / 0.0 TPM" to "Fuchs 2024"; the numbers are Chen 2021 (GSE150338).
4. **Axial wording:**
   - "The three PT segments (PTS1, PTS2, and PTS3) exhibited comparable expression patterns, which differed markedly from those of other tubule segments"

#### Chen 2023 — Multiomics reveal sex differences in mouse renal proximal subsegments (PMID 36758122)
1. **Access:** pre-converted PMC HTML text (txt/PMC10125651.txt). **Full text: Y.**
2. **Our description:** "mouse PTS1/2/3 TPM, male and female"; fig5: "Chen 2021 male and Chen 2023 male + female mouse microdissection".
3. **Verdict: SUPPORTED, but the male arm is not independent of Chen 2021.**
   - (Methods) "Six to 8-week-old male and female C57BL/6 mice from Taconic were used throughout the study. The estrus cycles for female mice were not synchronized."
   - (Results) "11 samples [n=5, male; n=6, female] were from PTS1; 14 [n=8, male; n=6, female] were from PTS2; and 13 (n=6, male; n=7, female] were from PTS3"
   - (Methods) "RNA-seq data from male mouse kidneys from our previous study31 were downloaded from https://esbl.nhlbi.nih.gov/MRECA/Nephron/and GSE150338."
   - The male n (5/8/6) equal the PTS1/2/3 sample counts in GSE150338. Treat "Chen 2021 male" and "Chen 2023 male" as the same samples, and only the female arm as new.
   - Also: per-segment ATAC-seq. Proteomics is whole-kidney TMT: "We detected more than 5300 proteins across four male and five female kidneys". It is not segment-resolved.
   - (Abstract) "the major differences in the proximal tubules between the sexes were in the S2/S3 segments"
   - Volcano counts: 79 sex-specific transcripts in PTS1, 1177 in PTS2, 766 in PTS3.
4. **Axial wording:**
   - "A number of TFs showed sex-biased expression in proximal tubules, particularly in PTS2 and PTS3"
   - "Slc22a6 (Oat1) ... was primarily expressed in male PTS2, while Slc22a8 (Oat3) ... was in female PTS2."

#### Lee 2015 — Deep sequencing in microdissected renal tubules (PMID 25817355)
1. **Access:** pre-converted PMC HTML text (txt/PMC4625681.txt). **Full text: Y.**
2. **Our description:** "male rat microdissected S1/S2/S3 RPKM".
3. **Verdict: SUPPORTED.**
   - (Concise Methods) "Male Sprague–Dawley rats weighing 200–250 g". Assay: oligo-dT single-cell-style RNA-seq of microdissected tubules.
   - (Table 1) S1, S2 and S3 have 7 replicates each, with dissected lengths of 0.6–11.8, 2.2–13.0 and 1.8–23.6 mm.
   - (Results) "Gene expression levels of all genes determined by median reads per kilobase of exon model per million mapped reads (RPKM)". The reported value is the **median** RPKM.
   - The S1/S2/S3 anatomical criteria are only in the Fig. 1A image and the Supplement, not in the text. The same lab's later definitions (Chen 2021, Limbutara 2020) are labyrinth / medullary ray / outer stripe.
4. **Axial wording:** none beyond "proximal region [S1, S2, and S3]".

#### Limbutara 2020 — Quantitative proteomics of all 14 renal tubule segments in rat (KTEA) (PMID 32358040)
1. **Access:** pre-converted PMC HTML text (txt/PMC7269347.txt). **Full text: Y.**
2. **Our description:** "rat S1/S2/S3 protein copies/cell".
3. **Verdict: SUPPORTED.**
   - (Methods) "male Sprague–Dawley rats age 4–8 weeks". Separate rats were used for cortical and medullary dissections.
   - (Table 1) S1 is the "First segment of proximal tubule" from the cortical labyrinth; S2 is from the cortical medullary ray; S3 is from the outer stripe of the outer medulla. Each has 3 replicates.
   - (Results) "Expression levels were quantified as copies per cell (cpc) using the Proteomic Ruler technique."
   - (Results) "Over the 44 samples, pooled from 90 rats, the average number of quantified proteins was 4234."
   - Caution: proteins quantified per segment are S1 2734 ± 103, S2 4750 ± 415, S3 3508 ± 287. A protein missing from S1 may be under-detection, which biases S1→S3 slopes for low-abundance proteins.
4. **Axial wording:** Table 1 origins: "Cortical labyrinth", "Cortical medullary ray", "Outer stripe of outer medulla".

#### Lake 2023 — Atlas of healthy and injured cell states and niches in the human kidney (PMID 37468583)
1. **Access:** pre-converted PMC text (txt/PMC10356613.txt; Europe PMC JATS). **Full text: Y.** Supplementary tables were not read.
2. **Our description:**
   - fig5: "Lake 2023 human PT-S1/S2/S3 snRNA per donor (PMID 37468583; 17 normal, 12 AKI, 23 CKD donors with ≥ 20 cells per segment)".
   - fig5lit: "13 normal cortex donors (7 male), plus 10 AKI and 17 CKD donors".
3. **Verdict: SUPPORTED as a resource.** The donor counts come from CELLxGENE v1.5 with our own filters, so they need not match the paper.
   - Cohort sizes in the paper:
     - (Abstract) "healthy reference kidneys (45 donors) and diseased kidneys (48 patients)"; this counts all technologies.
     - (Main) "from 58 reference tissues (35 donors) and 52 diseased tissues (36 patients) that covered the spectrum of conditions from healthy to AKI and CKD"; this counts the omics integration.
     - (Methods) scCv3 used "15 CKD,12 AKI and 18 living donor biopsy cores".
     - Per-assay snCv3 donor counts are only in Supplementary Table 1.
   - Assays: snCv3 (10x snRNA-seq), scCv3 (10x scRNA-seq) and SNARE2 (joint RNA + ATAC), plus Slide-seq2, Visium and 3D imaging. The reference tissues are nephrectomies, deceased donors and living-donor biopsies; the snCv3 dissections covered cortex, outer medulla and inner medulla.
   - PT-S1/S2/S3 were defined as snCv3 clusters annotated with a standardized marker nomenclature (Supplementary Tables 4–5, ED Fig. 2a). They were cross-checked by correlation with "mouse bulk segmental RNA-seq data from Chen et al., 2021" (ED Fig. 2e) and with Ransick 2019 (ED Fig. 2h).
   - **S1 vs S2 separability:** (Methods) "Cell types that could not be accurately resolved (PT-S1/PT-S2) were kept merged." This applies to scCv3 cells integrated onto snCv3. S1 and S2 are separate only in the snRNA (snCv3) data. Per-donor S1-vs-S2 values must come from snCv3 nuclei, not scCv3 cells.
   - **Flags:**
     - (a) Our two documents give different donor counts (17/12/23 vs 13/10/17). Reconcile them, for example by stating the filters: all regions vs cortex-only, ≥ 20 cells.
     - (b) Lake's PT segment labels were partly validated against Chen 2021 mouse segments. Using Lake to "independently" confirm a Chen-based mouse/human contrast is mildly circular for segment identity, though not for expression direction.
4. **Axial / injury wording:**
   - "We specified normal cell populations as the end points for trajectory inference (S1–S3 in PT and M-TAL in TAL)"
   - "UMAP of PT subclasses (PT-S1-S3, aPT) with projected RNA velocities" (ED Fig. 9 legend)

#### Xiong 2023 — Direct androgen receptor control of sexually dimorphic gene expression in the mammalian kidney (PMID 37673062)
1. **Access:** pre-converted PMC author-manuscript text (txt/PMC10873092.ncbi.txt). **Full text: Y.** Supplementary Tables S2 and S3 (the AR-responsive gene lists) were not read here.
2. **Our description:**
   - fig5lit: "per-PT-segment female-vs-male and nephron Ar-knockout-vs-WT tables".
   - fig5: "Ar knockout confirms androgen control (Hao2, Sh3gl2, Lrig1)".
3. **Verdict: PARTLY.** The design matches our description. Of our genes, only Hao2 is named in the paper text.
   - (Results) "we generated nephron-specific removal of Ar (Six2-Ar-KO) and Esr1 (Six2-Esr1-KO), then assayed kidneys at 8 weeks."
   - (Results) "Hao2 was highly expressed in female PT-S2 and PT-S3, but only low levels of expression were detected in homologous male PT segments (Fig. 6C)"
   - (Results) "a marked increase in PT-S3 Hao2 expression was observed on AR removal from the male nephrons (Fig. 6C)"
   - Design details:
     - Six2TGC/+ Cre × Ar flox (B6.129S1-Artm2.1Reb/J, JAX 018450). "Six2-CRE activity in nephron progenitor cells removes any potential sex hormone input from nephron progenitors and their nephron derivatives from the onset of embryonic kidney development".
     - Both sexes were knocked out: "nephron removal of either Ar or Esr1 had little effect on gene expression in the female kidney".
     - Assays: whole-kidney bulk RNA-seq ("performed between 8-12 weeks", Fig. 3 legend) and 10x single-nucleus RNA+ATAC multiome of WT and Six2-Ar-KO "at 8-10 weeks of age". The Methods give 9–12 weeks for the C57BL/6J multiome mice.
     - Also included: whole-body Sox2-Ar-KO, castration, ovariectomy, testosterone injection, and GTEx plus 5-donor human snRNA.
   - Per-segment tables: differential expression was computed "for each PT segment separately" (sSeq on meta-cells; Tables S3.1–3.2).
   - **Not named anywhere in the main text or legends:** Sh3gl2, Lrig1, Ggact, Odc1, Slc7a13 (AGT1) and Slc22a7 (OAT2). Their AR-responsiveness in fig5 must be sourced to Table S2.10 / S3.2 specifically, and the support for Sh3gl2 and Lrig1 still needs checking there.
4. **Axial wording:**
   - "PT segment 1 (PT-S1) cluster comprised S1 cells from all genotypes and sexes, indicative of a low level of AR-dependent variability between male and female sexes"
   - "each PT segment adopts a distinct segmental identity, with a segment-specific female ground-state that is masculinized by the direct action of Ar"
   - "PT-S3 showed the highest number of sex-biased genes (Fig. S4C)."

### Group C — Fig. 6 — injury, proliferation, metabolism and carrier loss

#### Liu 2017 — Molecular characterization of the transition from acute to chronic kidney injury following ischemia/reperfusion (PMID 28931758)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC5612583.txt). Table 1 is an image, downloaded from the PMC CDN and read visually. I also read the GEO GSE98622 series record and its sample list (ncbi.nlm.nih.gov/geo, text view). Full text: Y.
2. **Our claim:** "the Liu 2017 IRI whole-kidney bulk time course (GSE98622, 48 h against sham)"; "direction confirmed in Liu 2017 bulk (Slc7a13 −4.4, Slc3a1 −1.4, Slc22a2 −1.2, Slc6a13 −1.0, Slc46a1 −1.0)"; Cyp24a1 "IRI bulk up (Liu 2017)"; Agt "Liu 2017 −1.0"; Cyp51 "Liu 2017 gives −2.2 at 48 h"; MCT1 "bulk induction about 3–5× (Liu 2017; Balzer 2022)".
3. **Verdict: PARTLY (data description).**
   - (Methods, "IRI.") "Warm renal ischemia-reperfusion injury was performed on 10- to 12-week-old (25–28 g) C57BL/6CN male mice from Charles River."
   - (Methods, "IRI.") "...both renal pedicles, which were clamped for 21 minutes using nontraumatic microaneurysm clips"; body temperature "maintained at 36°C".
   - (Results, "RNA-seq analysis identified...") "whole-kidney total mRNA sequencing of 3 to 4 biological replicates at each time point"; "complemented by sham IRI surgery controls (4 hours, 24 hours, and 12 months)".
   - GEO sample titles: IRI48h-R1..R3; SHAM4h, SHAM24h and SHAM12m with 3 each; NORM3m, NORM9m and NORM15m. The assay is Illumina TruSeq polyA RNA-seq (HiSeq 2000), analysed with FPKM and DESeq2.
   - **Corrections.**
     - (a) No sham exists at 48 h. "48 h against sham" must name the reference (SHAM24h, or the pooled 4 h/24 h shams).
     - (b) Of our genes, the paper names only Slc7a13 (Table 1: module VII, downregulated PT genes) and Havcr1 (module III, "dramatically increased at day 2").
     - (c) Slc3a1, Slc22a2, Slc6a13, Slc46a1, Slc16a1, Cyp24a1, Agt and Cyp51 do not appear in the main text or Table 1. I did not check the supplementary module lists.
     - (d) All of our fold changes come from our own reanalysis of GSE98622. Write "in GSE98622 (Liu 2017 data)", not "Liu 2017 reports".
     - Minor inconsistencies in the paper itself: the pilot used 9–11 wk mice and the Methods say 10–12 wk; the strain is written C57BL/6NCrl and C57BL/6CN. The Methods sentence "mice with clamping of the pedicles were used as sham surgery controls" is a typo; GEO says "except for clamping".
4. **Axial / injury wording we can borrow:**
   - (Results) "fibrosis was restricted to the outer medullary region adjacent to damaged S3 segments"
   - (Results) "significant tubular damage in the outer medullary region, as evidenced by cell flattening, sloughing of tubular epithelium"

#### Kirita 2020 — Cell profiling of mouse acute kidney injury reveals conserved cellular responses to injury (PMID 32571916)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC7355049.txt), plus the PMC Fig. 2 image read visually. Full text: Y (main text and figures). Supplementary Datasets S1–S3 were not read.
2. **Our claim:** cited in the manuscript for single-cell AKI states. An earlier review said "Slc16a1 and Cyp24a1 are absent from Kirita 2020", which supports fig6's "Not reported" positional claims.
3. **Verdict: PARTLY.**
   - (Methods, "Surgery.") "For bilateral IRI, 8- to 10-wk-old male mice were anesthetized with isoflurane"; "clamping the renal pedicle with a nontraumatic microaneurysm clamp (Roboz, Rockville, MD) for 18 min."
   - (Results) "Mice were euthanized at 4 and 12 h and 2, 14, and 42 d after bilateral ischemia–reperfusion injury (IRI)". The Fig. 1 legend gives "n = 3 mice per group." The assay is snRNA-seq.
   - (Results, "Proximal Tubule Responses...") "three healthy subclusters (the S1, S2, and S3 segments of the proximal tubule), one repairing subcluster, and three injured subclusters". These are defined by unsupervised subclustering, not spatially.
   - (Results) "“Injured S1/2” and “injured S3” were primarily composed of cells from the 4- and 12-h time points and expressed Myc ... and Havcr1"; "the proximal tubule injury clusters had almost disappeared by 2 d".
   - **Corrections.**
     - Havcr1 is present (main text and Fig. 2A). The Fig. 2A marker dot plot shows Slc5a12, Slc22a30, Slc7a13, Bcat1, Top2a, Myc, Havcr1, Krt20, Hspa1a, Vcam1, Dcdc2a, Sema5a and Ccl2.
     - Slc16a1 and Cyp24a1 are not named in the main text, the legends or the Fig. 2A panel. The DEG datasets (S1–S3) were not checked, so state "not named in Kirita 2020's main text or figures".
     - No sentence singles out S3 as most injured. The paper says PT "suffers the most injury due to high metabolic activity".
     - The text says the injured clusters "shared a part of differentially expressed genes (DEGs) of 'healthy S1' and 'healthy S2' ... respectively". That pairs injured S3 with healthy S2, so do not cite it for S3 identity.
4. **Axial / injury wording we can borrow:**
   - (Results) FR-PTC "down-regulated expression of terminal differentiation markers such as Slc5a12, Slc22a30, and Slc7a13 even at late time points"
   - (Results) "most cells progressed to healthy S1/S2 or S2/S3 proximal tubule segments"

#### Li 2022 — Comprehensive single-cell transcriptional profiling defines shared and unique epithelial injury responses during kidney fibrosis (PMID 36265491)
1. **Access:** PMC article (pre-converted txt/PMC9742301.ncbi.txt). Full text: Y.
2. **Our claim:** "PT carrier loss in injured states (Li 2022; Wu 2024; ...)". Earlier quote: injured PT states show "downregulation of healthy PT marker genes such as solute-linked carriers (e.g., Slc34a1, Slc5a12, Slc7a13)".
3. **Verdict: SUPPORTED.**
   - (Results, "Generation of two mouse models...") "We performed uni-IRI and UUO surgeries on 8–9 week old adult male C57BL6/J mice"; "(0, 6 hours and 2, 7, 14, 28 days post uni-IRI ... n = 2 for each timepoint)".
   - (Methods) "ischemia was induced on the left kidney by clamping the renal pedicle with a nontraumatic micro aneurysm clamp (RS-5420, Roboz) for 22 minutes".
   - (Results, "Diverse Proximal Tubule Injury States") "All injured PT cell states were also characterized by downregulation of healthy PT marker genes such as solute-linked carriers (e.g., Slc34a1, Slc5a12, Slc7a13), suggesting cell dedifferentiation."
   - **Notes.**
     - IRI is included, but it is **unilateral** (22 min) and paired with UUO.
     - The relevant time points are 6 h and 2 d: "Type1 injured cells comprised 80% of all PT cells at 6 hrs after injury, with this proportion falling rapidly to 5% of the total by D2".
     - The assay is sci-RNA-seq3 on nuclei. The named carriers are not our S3 handoff partners, and there is no segment or position resolution.
4. **Axial / injury wording we can borrow:**
   - (Results) "This revealed the expected healthy S1/S2/S3 PT clusters, as well as multiple injury states (high expression of Havcr1 and Nrg1)"

#### Wu 2024 — High resolution spatial profiling of kidney injury and repair using RNA hybridization-based in situ sequencing (PMID 38360882)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC10869771.txt). Full text: Y.
2. **Our claim:** same as Li 2022. Earlier quote: PTS1/2/3 could not be separated in integrated data "because injury causes downregulation of PT segment-defining markers".
3. **Verdict: SUPPORTED.**
   - (Results, "dRNA HybISS detects spatial dynamics...") "this difference was not captured when clustering the combined dataset because injury causes downregulation of PT segment-defining markers (Fig. 5A and Supplementary Fig. 15)".
   - (Methods, "Animals") "All in vivo experiments were performed on 8- to 10-week-old C57BL6/J male mice". (Methods, "Bilateral IRI surgery") "Both right and left kidneys were clamped with nontraumatic microaneurysm clamps (RS-5420, Roboz) for 23 min."
   - (Methods) "mice were euthanized after Bi-IRI for sham control or 4 h, 12 h, 2 days, and 6 weeks post-surgery."
   - The platform is dRNA HybISS (Cartana), a 200-gene panel. The healthy comparison includes female mice, but the IRI time course is male.
4. **Axial / injury wording we can borrow:**
   - (Results) "the spatial distribution of the injPT cell state and the injury marker Havcr1 expression extends throughout the cortex at early timepoints". This places early KIM-1 in the cortex, not only in S3.

#### Xuanyuan 2025 — Multimodal spatial transcriptomic characterization of mouse kidney injury and repair (PMID 40813851)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC12354736.txt). Full text: Y.
2. **Our claim:** "injury spreads from the corticomedullary junction (Xuanyuan 2025, maps only)"; "Xuanyuan 2025: healthy-S3 neighbourhood absent at 12 h–2 d".
3. **Verdict: SUPPORTED.**
   - (Results, "Spatial transcriptomics accurately recovers...") "These injured PT cells were initially localized within the corticomedullary junction and progressively expanded into the cortical region by 12 hours (Fig. 2A)."
   - (Results, "Spatial organization and cell neighborhood dynamics...") "CN3, primarily composed of healthy PTS3 cells, was absent during the mid-injury phase (hour12, day2), consistent with the vulnerability of S3 proximal tubule segments to ischemic injury (Fig. 4B, C)."
   - (Methods, "Animal surgery") "clamping the renal pedicles ... for 18 minutes"; "Only male mice were used in this study". It is bilateral, in C57BL/6J mice aged 8–10 wk.
   - (Results) The six time points are "sham, 4 hours, 12 hours, 2 days, 14 days, and 6 weeks" with "two technical replicates for each time point". The panel is Xenium with 300 genes, plus Visium CytAssist.
   - **Notes.**
     - The corticomedullary statement is in the text, but it rests on maps of selected regions (Fig. 2A) with no per-zone quantification, so "maps only" is fair.
     - Inj-PT rises from 15.3% (4 h) to 33.2% (12 h) of total cells.
     - The Methods time list omits 14 d, which is internally inconsistent.
4. **Axial / injury wording we can borrow:**
   - (Results) "Proximal tubule cells were prominently located in the cortex and outer stripe of the medulla, with S3 segments displaying intense brush border PAS staining and high expression of Slc7a13."
   - (Results) PT subclusters: "three healthy segments (Healthy S1, S2, S3), three injured segments (Inj S1, S2, S3), and one failed-repair PT population"

#### Humphreys 2008 — Intrinsic epithelial cells repair the kidney after injury (PMID 18371453)
1. **Access:**
   - Not in PMC: elink gives no PMC record, and a PMC title search returns 0.
   - cell.com and ScienceDirect return 403. Semantic Scholar lists it as Bronze OA at cell.com.
   - Obtained the publisher's own free PDF via a Wayback Machine snapshot: web.archive.org/web/20220804070333/https://www.cell.com/cell-stem-cell/pdf/S1934-5909(08)00055-6.pdf (13 pages, text layer).
   - The PII given in the task (...647) is wrong; the correct PII is S1934590908000556.
   - Full text: Y.
2. **Our claim:** "Proliferation peaks in outer stripe and S3 at 2 d (Humphreys 2008, 2011; Witzgall 1994, abstract only)".
3. **Verdict: PARTLY.**
   - (Results) "In uninjured kidneys from wild-type mice, 0.4% ± 0.1% of tubular cells in the S3 segment of the proximal tubule express the proliferation marker Ki67."
   - (Results) "Two days after injury, 47.4% ± 6.0% express Ki67, and this fraction falls to 5.4% ± 2.3% at 7 days after injury and 2.6% ± 0.9% 15 days after injury"
   - (Fig. 3A legend) "Ki67-positive tubular nuclei (red) were identified in outer medulla of kidney sections at baseline and at 2, 7, or 15 days after kidney injury."
   - **Correction.**
     - Counts were made only in S3 / outer medulla, with no cortex S1/S2 comparison, so this paper does not show that S3 exceeds other segments.
     - 2 d is the earliest time sampled; there is no 1 d or 3 d point.
     - Model: (Methods) clamps "removed after 25 min (males) or 35 min (females)", so both sexes are allowed. Unilateral IRI was used for BrdU and histology, and bilateral IRI for creatinine. The model for the wild-type Ki67 time course is not stated. The transgenics are mixed 129/B6.
4. **Axial / injury wording we can borrow:**
   - (Results) "50.5% ± 5.3% of tubular epithelial cells in the S3 segment of the proximal tubule coexpressed RFP and the proliferation marker Ki67"
   - (Fig. 3D legend) "In more severely injured outer medullary tubules, flattened RFP-positive cells lining the tubule are Ki67 positive"

#### Humphreys 2011 — Repair of injured proximal tubule does not involve specialized progenitors (PMID 21576461)
1. **Access:** PMC article HTML (pre-converted txt/PMC3107336.txt). Full text: Y.
2. **Our claim:** as for Humphreys 2008.
3. **Verdict: PARTLY.**
   - (Methods) "Male C57Bl6 mice, aged 8 to 10 wk"; clamps "removed after 27 min".
   - (Fig. 2 legend) "unilateral IRI was performed, and mice were killed at either 12, 24, 36, or 48 h after injury." Bilateral IRI appears only in Fig. S1, at 48 h.
   - (Results) "substantial IdU incorporation could be detected beginning 24 h after injury and increasing thereafter, from 0.3 ± 0.1% in uninjured cortex to 26.9 ± 4.1% of tubular cells by 48 h"
   - (Methods) "The combined outer medulla and inner cortex were assessed because this is the region with most proliferative activity after IRI."
   - **Correction.** Proliferation is still rising at the last time point (48 h), so no peak is identified. The region is pooled (outer medulla + inner cortex), with no S1/S2 vs S3 split. The paper supports "high and rising at 2 d in outer medulla/inner cortex", not "peaks at 2 d".
4. **Axial / injury wording we can borrow:**
   - (Results) "at least 55% (CldU% + Ki67% − copositive%) of outer medullary epithelia have either reentered the cell cycle (Ki67 expression) or completed the S-phase"
   - (Results) "75.0 ± 18.1% of CldU-positive cells and 72.3 ± 8.9% of IdU-positive cells coexpressed Kim-1"

#### Witzgall 1994 — Localization of PCNA, vimentin, c-Fos, and clusterin in the postischemic kidney (PMID 7910173)
1. **Access:** the PMC citation_pdf_url (pmc.ncbi.nlm.nih.gov/articles/PMC294357/pdf/jcinvest00034-0311.pdf) returned a "Preparing to download" HTML page, and the Europe PMC render returned 403. The JCI site PDF (www.jci.org/articles/view/117214/files/pdf, 14 pages) worked and has an OCR text layer, which I read. Full text: Y.
2. **Our claim:** as for Humphreys 2008, with "Witzgall 1994, abstract only".
3. **Verdict: SUPPORTED (in rat).**
   - (Methods) "190-275-g male Sprague-Dawley rats"; "A left flank incision was made and the renal artery and vein were clamped"; clamp removed "40 min later". This is unilateral with the contralateral kidney intact.
   - (Results, "Expression of PCNA") "PCNA-positive nuclei could already be detected after 24 h of reperfusion, but the maximum level of expression was seen in the 48-h reperfused kidneys (Fig. 2)."
   - (Fig. 3 legend) "The number of stained nuclei is much greater in the S3 segment of the proximal tubule in the outer stripe, than in S1 and S2 segments located in the cortex."
   - **Correction.**
     - The species is rat, the model is unilateral 40-min ischemia, and the marker is PCNA (late G1/early S).
     - Fig. 2 legend: "At 120 h all cells were negative for PCNA."
     - Remove "abstract only".
4. **Axial / injury wording we can borrow:**
   - (Fig. 6 legend) "Some proximal tubules show an abrupt transition from negative to positive vimentin staining ... interpreted to represent a transition from the S2 (negative) to S3 (positive)"
   - (Discussion) "a specific cellular response (S3 > S1, S2) to ischemia"; also (Fig. 6 legend) "Occasional vimentin-positive cells ..., presumably S3-like cells, were found in otherwise negative regions of S2"

#### De Chiara 2022 — Tubular cell polyploidy protects from lethal AKI but promotes consequent CKD (PMID 36195583)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC9532438.txt). Full text: Y.
2. **Our claim:** "S-phase transcripts can mean endoreplication. Polyploid cells are mostly cortical S1/S2 (De Chiara 2022)".
3. **Verdict: SUPPORTED.**
   - (Results) "Bi-coloured polyploid TC localized mostly in the cortex (13.3 ± 2.5% over total coloured TC in the cortex vs 2.7 ± 0.4% in the outer stripe of outer medulla"
   - (Results) "Bi-coloured polyploid TC mostly localized to S1 and S2 segments of the cortex sparing the injured outer medulla"
   - (Results) "at day 3, 17 ± 2.4% of TC (i.e., 38.8 ± 9.2% of the cells that were cycling at day 2) had acquired ≥4C DNA content"
   - **Model.**
     - (Methods) "male mice underwent a unilateral ischemia reperfusion injury (IRI, ischemic injury) of 30 min, and were then sacrificed at day 2, 3, 5 and 30 days post-IRI." This is left renal artery clamping with the contralateral kidney untouched.
     - The Confetti segment images are from day-30 kidneys (Fig. 1i–m).
     - Segments were assigned by AQP1 plus "morphology and localization in distinct kidney areas".
     - The glycerol model was used only in females.
4. **Axial / injury wording we can borrow:**
   - (Results) "ischemic necrosis affects specifically S3 segments of the proximal tubules in the outer medulla raising the question why AKI involves widespread cell cycle marker positivity all over the cortex"

#### Lazzeri 2018 — Endocycle-related tubular cell hypertrophy and progenitor proliferation recover renal function after AKI (PMID 29632300)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC5890293.txt). Full text: Y.
2. **Our claim:** "at 2 d only about 13% of Ki-67+ cells divide (Lazzeri 2018)".
3. **Verdict: PARTLY. The direction is right, but the number and its denominator are misread.**
   - (Results, "Cell cycle markers do not predict cell division after AKI") "few TECs appeared in S/G2/M (mVenus+) at day 2 after IRI despite widespread tubular Ki-67 positivity (13.2 ± 3.4% vs. 47.1 ± 9.2%"
   - (Results, "Pax8+ TECs endocycle...") "Pax8/FUCCI2 mice 2 days after IRI showed 9.8 ± 0.7% of Pax8+ TECs in S phase, but only 2 ± 0.5% in G2/M phase"
   - (Methods) "Male mice ..."; "The duration of ischemia was 30 min"; left renal artery; "The right contralateral kidney was maintained untouched as a control." These data are IRI. Glycerol AKI was used in females elsewhere in the paper.
   - **Correction.** 13.2% is the fraction of all FUCCI2-labelled tubular cells in S/G2/M, compared with 47.1% of tubular cells that are Ki-67+. It is not the fraction of Ki-67+ cells that divide, and S/G2/M is not division. Suggested wording: "2 d after unilateral IRI, 47% of tubular cells are Ki-67+, but only 13% are in S/G2/M by a FUCCI reporter and 2% in G2/M by flow cytometry (Lazzeri 2018)."
4. **Axial / injury wording we can borrow:**
   - (Results) "Pax2+ cells in the OSOM are the cells that regenerate injured tubules after AKI by mitotic cell division, while other TECs of the S1 and S2 segments rather enter endocycles". These are day-30 data: endocycling 22.3% in cortex vs 7.8% in OSOM.

#### Tran 2016 — PGC1α-dependent NAD biosynthesis links oxidative metabolism to renal protection (PMID 26982719)
1. **Access:** Europe PMC JATS XML (pre-converted txt/PMC4909121.txt). Full text: Y, including the Extended Data legends.
2. **Our claim:** "PGC-1α and OXPHOS collapse after IRI (Tran 2016; ...)"; for Flvcr2, "Choline transporter; choline depleted after ischaemia (Tran 2016)".
3. **Verdict: PARTLY for PGC-1α/OXPHOS; SUPPORTED for choline.**
   - (Results) "Following transient local ischemia, renal function worsened, PGC1α expression declined, tubular mitochondria swelled, and a pronounced accumulation of acylglycerols developed in tubules". Extended Data 1c: "Renal PGC1α expression 24 h after sham or IRI (n = 5 animals per group)".
   - (Results) "Four were shared between settings, with all four lower in PGC1α−/− and post-ischemic kidneys (Fig 1g,h, Extended Data 3a,b)." Also: "Deficiency of betaine and choline, two osmolytes essential for cell volume maintenance..."
   - (Methods) "Ischemia-reperfusion injury (IRI) was performed on 8–12 week-old males ... applying a microvascular clamp to each renal pedicle for 20 minutes". This is bilateral. Kidneys were "collected 24 h after the AKI model".
   - **Correction.**
     - The paper shows PGC1α falling after IRI, but it does **not** show OXPHOS genes or proteins falling after IRI. Its COX activity, COX IV and Ndufs1/Cycs/Atp5o data compare iNephPGC1α mice with controls.
     - The choline result is from whole-kidney homogenate metabolomics in sham vs IRI (n = 6/group, 24 h). It is not segment-resolved, and Flvcr2 is not mentioned.
4. **Axial / injury wording we can borrow:**
   - (Results) mitochondrial abundance has an "anatomical distribution favoring cortex and outer stripe of the outer medulla"

#### Funk & Schnellmann 2012 — Persistent disruption of mitochondrial homeostasis after acute kidney injury (PMID 22160772)
1. **Access:** PMC article HTML (pre-converted txt/PMC3340936.txt). Full text: Y, so fig6's "abstract only" is wrong.
2. **Our claim:** "PGC-1α and OXPHOS collapse after IRI (Tran 2016; Funk 2012, abstract only)".
3. **Verdict: PARTLY. The respiratory-chain loss is supported; PGC-1α loss is contradicted.**
   - (Methods) "Eight-week-old male C57BL6 mice weighing 25–30 g were subjected to bilateral renal pedicle ligation"; "occluded with a vascular clamp for 20 min, and mice were euthanized at 24, 72, or 144 h".
   - (Results) "Immunoblot analysis on renal cortical lysates revealed that NDUFB8, ATP synthase β, and COX I protein levels decreased after I/R injury and did not recover over 144 h". ND6, COX I and ATP synthase β mRNA also fell from 24 h.
   - (Results) "PGC-1α, PGC-1β, Tfam, and NRF-1 mRNA tended to decrease slightly immediately after injury. In contrast to mRNA expression, protein levels of PGC-1α and Tfam increased after injury"
   - (Results) "tissue ATP levels were elevated 24 and 72 h after I/R injury and returned to control levels at 144 h (Fig. 7B)."
   - **Correction.** Cite Funk 2012 only for loss of respiratory-chain proteins in cortex and PT (24–144 h). Do not cite it for PGC-1α collapse (the protein rises) or for ATP loss (cortical ATP rises).
4. **Axial / injury wording we can borrow:**
   - (Results) "Immunoreactivity of COX IV within the kidney cortex, particularly with proximal tubules, became progressively less intense 72 and 144 h after injury"

#### Fujii 2019 — Xanthine oxidase inhibitor ameliorates postischemic renal injury in mice by promoting resynthesis of adenine nucleotides (PMID 31723053)
1. **Access:** PMC article HTML (pre-converted txt/PMC6948864.txt). Full text: Y, so "abstract only" is wrong.
2. **Our claim:** "zone-level ATP loss (Fujii 2019, abstract only; ...)". Earlier review: "ATP falls in cortex and outer stripe but not inner stripe".
3. **Verdict: SUPPORTED, with a model caveat.**
   - (Abstract) "After 10-minute ischemia and 24-hour reperfusion (10mIR), in the cortex and outer stripes of the outer medulla, ATP substantially decreased, but not in the inner stripes of the outer medulla and inner medulla."
   - (Methods) "Male C57BL/6 mice were purchased at 8 weeks of age"; "Renal ischemia was produced by left renal artery occlusion using artery clips for 1 or 10 minutes."
   - (Results) "Total adenylates in the cortex and OSOM that were depleted even after the 24-hour reperfusion were notably recovered by the febuxostat treatment."
   - **Notes.**
     - The method is semiquantitative MALDI imaging mass spectrometry calibrated by CE-MS, not ATP-sensor mice.
     - The model is a mild, unilateral 10-min ischemia: "The brush border of the proximal tubular cells was well maintained".
     - Readouts are by zone, not by segment.
4. **Axial / injury wording we can borrow:**
   - (Results) "The physiological accumulation of adenosine in the OSOM disappeared after the clipping procedure"

#### Arnoux 2025 — Integrated spatial and functional metabolic profiling identified the TAL as a mitochondrial vulnerability hub in AKI (PMID 40850621)
1. **Access:** UNIGE Archive ouverte record archive-ouverte.unige.ch/unige:187670, downloaded as "Article (Published version) CC BY-4.0" via access.archive-ouverte.unige.ch/access/metadata/e8567904-.../download (18-page PDF). Fig. 3 was rendered at 400 dpi and read visually. Full text: Y.
2. **Our claim:** "zone-level ATP loss (...; Arnoux 2025 names TAL as most affected)". Earlier review: "all regions decline 12–48 h with the TAL most affected".
3. **Verdict: PARTLY.**
   - (Results, "Kidney metabolism is spatially altered...") "In the early phase (12–48 hours postinjury, depending on the nephron segment), all regions exhibited a marked reduction in metabolic activity, with the TAL showing the most pronounced decline (Figure 3a)."
   - (Methods) "we used gene sets from the Reactome database as input for the singscore algorithm." (Fig. 3a legend) "changes in metabolic processes relative to control (CTRL) at 4 hours, 12 hours, 48 hours, and 6 weeks after ischemia-reperfusion injury (IRI) across distinct renal cell types."
   - (Methods) "Eight-week-old male C57Bl/6JRj mice"; right nephrectomy plus a left thread clamp; "After 25 minutes of ischemia, the clamp was released". This model was used for the propofol experiment.
   - **Correction.**
     - The IRI result is a transcriptomic pathway score per cell type, not a measurement of ATP (or any metabolite) by zone. The text does not name the source of the time series. Its time points (4 h, 12 h, 48 h, 6 wk) match an snRNA-seq IRI course, and GSE139107 is cited only as the annotation reference.
     - In Fig. 3a (read visually), the ETC scores for PTS1, PTS2 and PTS3 are **not reduced** at 4–48 h (light positive), while TAL, DCT, CNT, PC and IC are reduced. TCA in PTS3 rises at 12–48 h.
     - So Arnoux supports "TAL most affected" but **argues against** PT/S3 respiratory-chain loss at the cluster level. It is a contrast to state, not support.
4. **Axial / injury wording we can borrow:**
   - (Results) "the S3 segment of the proximal tubule and TAL were abundant in the OSOM"
   - (Results) "the TAL undergoes progressive metabolic decline, highlighting its intrinsic vulnerability, whereas proximal tubules exhibit adaptive metabolic enhancement, especially within the cortical region."

#### Matsuzaki 2008 — Altered pharmacokinetics of cationic drugs caused by down-regulation of rOCT2 (Slc22a2) and rMATE1 in I/R-induced AKI (PMID 18180268)
1. **Access:**
   - dmd.aspetjournals.org, the Elsevier/ScienceDirect page (PII S0090955624016180) and Semantic Scholar (CLOSED) failed.
   - Obtained the publisher PDF, which was free after 12 months, via a Wayback snapshot: web.archive.org/web/20170815180115/http://dmd.aspetjournals.org/content/dmd/36/4/649.full.pdf (6 pages, text layer).
   - Full text: Y.
2. **Our claim:** "OCT2 and OAT down after IRI" (no author named). Earlier review: "OCT2 protein down in rat I/R (Matsuzaki 2008)".
3. **Verdict: SUPPORTED (rat).**
   - (Methods) "Male Sprague-Dawley rats, initially weighing 200 to 210 g"; "Renal ischemia was induced using vascular clamps ... over both pedicles for 30 min."; "All experiments were performed under surgical anesthesia at 48 h after I/R."
   - (Results) "rOCT2 protein expression was markedly suppressed in I/R rat kidney compared with that in sham-operated rat kidney, whereas there was no significant difference in the expression of rOCT1."
   - (Results) "the levels of rOCT1 and rOCT2 mRNA in I/R rat kidneys were significantly depressed to 48 and 4%, respectively." Also: "The expressions of rOAT1 and rOAT3 protein were significantly depressed in I/R rat kidney".
   - **Notes.** Western blots are on whole-kidney homogenates (n = 3) and qPCR is whole kidney. Both mRNA and protein fall, and there is no regional or segment resolution.
4. **Axial / injury wording we can borrow:**
   - (Discussion) rOCT1 and rOCT2 "are localized to the basolateral membranes of renal proximal tubules"

#### Schneider 2007 — Downregulation of OAT1 and OAT3 correlates with impaired PAH secretion after ischemic ARF in rats (PMID 17244891)
1. **Access:**
   - journals.physiology.org returns 403, and Semantic Scholar lists it as CLOSED.
   - Obtained the publisher's full-text HTML, which was free after 12 months, via a Wayback snapshot: web.archive.org/web/20080907095533/http://ajprenal.physiology.org/cgi/content/full/292/5/F1599.
   - Full text: Y. Figures were not viewed; the legends and Results text were read.
2. **Our claim:** "OCT2 and OAT down after IRI". Earlier review: "OAT1/3 down (Schneider 2007)".
3. **Verdict: SUPPORTED (rat, female).**
   - (Methods) "I/R injury was induced by bilateral clamping of renal arteries for 45 min. Female Sprague-Dawley rats (200 to 250 g body wt)"
   - (Results) "the relative amount of mRNA for both OAT1 and OAT3 is strongly decreased 6 h and is totally restored after a 72-h postischemic reperfusion period."
   - (Results) "the relative expression of OAT1 is strongly reduced 6 h after reperfusion and then gradually recovers with time." Protein matches sham "7 days (168 h) after intervention".
   - **Correction.** The animals are female rats, and the model is a bilateral 45-min arterial clamp. The tissue is cortex extracts, with semiquantitative RT-PCR plus Western blot. mRNA is down at 6–24 h and recovered by 72 h, which matches our 12 h–2 d window. Protein lags.
4. **Axial / injury wording we can borrow:**
   - (Discussion) OAT1/3 are "the basolaterally located, rate-limiting ... exchange proteins". There is no segment wording beyond proximal tubule.

### Group D — Fig. 6 — Cyp24a1, Agt, the PPAR-alpha genes, MCT1 and the ZIPs

#### Fuchs 2024 - Intestinal Cyp24a1 regulates vitamin D locally, independent of renal Cyp24a1 (PMID 39688907)

1. **Access:** Europe PMC JATS full text (`/PMC11827884/fullTextXML`), converted to `txt/PMC11827884.txt`. Full text Y.
2. **Our claim:** fig6: "Cyp24a1 extends across S2 (control S1-restricted, AKI flat and high) | Normally S1-restricted (85.9 / 2.4 / 0.0 TPM; Fuchs 2024)". Earlier review quoted "Cyp24a1 was basally expressed in proximal tubules of the outer cortex" and expansion "across the cortex along all segments of the proximal tubule".
3. **Verdict: PARTLY.** The two quoted sentences are verbatim correct, but the paper assigns no PT segment and contains no TPM values.
   - (Results, basal and inducible expression) "In the kidney, Cyp24a1 was basally expressed in proximal tubules of the outer cortex (Figure 1A)."
   - (Results, same paragraph) "Cyp24a1 mRNA expression markedly increased and expanded across the cortex along all segments of the proximal tubule"
   - (Fig. 1 legend) "PT, proximal tubule (S1, S2, S3 segments)" - the scRNA-seq cluster pools the three segments, so no segmental resolution there either.
   - **Corrections.** Method is RNAscope Multiplex Fluorescent v2 (plus a BaseScope exon-5 probe for recombination), not a reporter. "Outer cortex" is an anatomical, not a segmental, statement - S1 is a reasonable reading but the paper does not say it. The 85.9 / 2.4 / 0.0 TPM triple is Chen 2021 microdissection and must be cited there. Sex: "WT mice of both sexes used for this study were procured from Jackson Laboratories and analyzed at about 12 weeks of age"; results from both sexes are reported combined where they did not differ, and the adenine-CKD experiments are male. Strain C57BL/6J. The induction shown is calcitriol (10 ug/kg, 6 h), not injury; the paper has no AKI or IRI model.
4. **Axial / injury wording we can borrow:**
   - "basally expressed in proximal tubules of the outer cortex"
   - "expanded across the cortex along all segments of the proximal tubule" - close to our own "spreads from S1 across S2", in a calcitriol rather than an injury setting.

#### Egli-Spichtig 2018 - FGF23 expression increased in multiple organs in folic-acid AKI (PMID 30405444)

1. **Access:** Europe PMC JATS full text (`/PMC6206018/fullTextXML`), `txt/PMC6206018.txt`. Full text Y.
2. **Our claim:** fig6: "Contradicted in folic-acid AKI (Egli-Spichtig 2018, down 17-fold), so model-specific".
3. **Verdict: SUPPORTED.**
   - (Results, iFGF23, vitamin D and inflammation in FA-AKI) "plasma 1,25(OH)2D increased by 3-fold in FA-treated mice which was accompanied by 6-fold increase in renal Cyp27b1 and 17-fold reduction in Cyp24a1 mRNA expression"
   - (Methods, Animals) "Eight-week-old male and female C57Bl/6J mice from in house breeding were used for all experiments. AKI was induced by intraperitoneal injection of folic acid (250 mg/kg ...)"; "After 24 h, mice were anesthetized ... before blood and organs were collected."
   - (Fig. 1 legend) "relative renal (F) Cyp27b1, (G) Cyp24a1 ... mRNA expression in vehicle (open triangle) and FA (black squares) treated mice after 24 h. Gus was used as housekeeping gene ... Five to eight mice per group."
   - **Corrections / limits.** Both sexes, 8 wk, C57Bl/6J; folic acid 250 mg/kg, single time point 24 h; whole-kidney qPCR, so no segmental or spatial information at all. The paper contains no ischemia-reperfusion data (no "ischem" or "reperfus" string anywhere), so the contrast with our 12 h-2 d IRI is model AND time-point AND assay. The paper's own framing is that FGF23 fails to regulate renal vitamin D metabolism in FA-AKI, which it attributes to loss of renal Klotho.
4. **Axial / injury wording we can borrow:** none - the paper has no axial language. Useful injury phrasing: "FA-AKI reduced Klotho mRNA expression and triggered renal inflammation".

#### Christov 2013 - Plasma FGF23 levels increase rapidly after acute kidney injury (PMID 23657144)

1. **Access:** Europe PMC JATS full text (`/PMC3766419/fullTextXML`), `txt/PMC3766419.txt`. Full text Y.
2. **Our claim:** fig6, unattributed in the row: "FGF23 rises in AKI" (used as mechanistic support for Cyp24a1 induction).
3. **Verdict: SUPPORTED for the FGF23 rise; NOT SUPPORTED as a Cyp24a1 link.**
   - (Abstract / Results) "plasma FGF23 levels increased significantly from baseline already after 1 hour of AKI, with an 18-fold increase at 24 hours"
   - (Results, human) "plasma FGF23 levels rose 15.9-fold by 24 hours following cardiac surgery in patients who developed AKI"
   - (Methods, Experimental Animals) "Male animals were studied. CD1 and C57BL/6 mice at 11-13 weeks of age"
   - **Models and timing.** Folic acid 240 mg/kg i.p. in CD1 mice (24 h, with a 0-8 h time course); glycerol rhabdomyolysis ("pigment nephropathy") in C57BL/6; genetic arms PPR-cKO, PTH-null, VDR-null; human nested case-control, 14 cardiac-surgery patients, 4 with AKI. **There is no IRI model.**
   - **Cyp24a1 is never measured** (no occurrence in the text). Moreover "There was no difference in mean plasma 1,25D levels between AKI and control animals 24 hours after FA-injection", and the rise is explicitly PTH-, 1,25D-, VDR- and phosphate-independent. So the chain FGF23 up therefore Cyp24a1 up in AKI is not supported by this paper, and Egli-Spichtig measures it going the other way in the same model. Recommend dropping "FGF23 rises in AKI" as mechanistic support for the Cyp24a1 expansion, or stating explicitly that the vitamin D axis in AKI does not follow the FGF23 expectation.
4. **Axial / injury wording we can borrow:** none axial. Useful: "FGF23 rises earlier in the course of AKI than phosphate or NGAL".

#### Matsusaka 2012 - Liver angiotensinogen is the primary source of renal angiotensin II (PMID 22518004)

1. **Access:** PMC article HTML with browser user agent, re-converted to `Dd/PMC3380650.txt`. Full text Y.
2. **Our claim:** fig6: "Agt mRNA is S2/S3 (Matsusaka 2012)"; notable result 3: "Agt is suppressed in S2 where it is normally highest".
3. **Verdict: NOT SUPPORTED for S2.** The paper assigns Agt mRNA to S3 only; S1/S2 is the protein.
   - (Results, Kidney Agt KO) "In control mice, in situ hybridization showed that Agt mRNA was expressed in proximal tubule cells of the S3 segment as reported previously."
   - (Fig. 1C legend) "In control kidneys, Agt mRNA is expressed in the outer stripe of the outer medulla."
   - (Fig. 3 legend) "Unlike Agt mRNA, which is predominantly expressed in the S3 segment (Figure 1C), Agt protein is mainly localized in the S1 and S2 segments."
   - (Introduction) "the intensity of renal Agt mRNA, confined to the proximal straight tubule"
   - **Corrections.** Method: ISH on control mouse kidney (and immunostaining/Western for protein). Animals: adult male mice on a C57BL/6N background, testosterone-treated ("adult male mice ... were treated with testosterone and used for the present studies"). The S1/S2 protein signal is megalin-dependent reabsorption of liver-derived Agt, not local synthesis: "Agt staining is strictly limited to proximal tubule cells expressing megalin". Visium HD measures mRNA, so the comparator for our S2 measurement is the S3/proximal-straight-tubule mRNA pattern. Our sentence "suppressed in S2 where it is normally highest" therefore rests on our own control curve, not on this paper; recommended rewording is that Agt mRNA is a late-PT (S3 / proximal straight tubule) transcript whose within-S2 rise our controls show and injury flattens.
4. **Axial / injury wording we can borrow:**
   - "expressed in proximal tubule cells of the S3 segment"; "confined to the proximal straight tubule"
   - "The cortex that mainly contains S1 and S2 segments of the proximal tubule was photographed" - a clean way to say cortex equals S1/S2.

#### Kontogiannis and Burns 1998 - Role of AT1 angiotensin II receptors in renal ischemic injury (PMID 9458826)

1. **Access:** INACCESSIBLE - abstract only. journals.physiology.org returns HTTP 403 to curl and to WebFetch (both the DOI PDF and the legacy ajprenal path); not in PMC (ID converter: "Identifier not found in PMC"); Semantic Scholar openAccessPdf status CLOSED; a title search turned up no repository or preprint copy. Judged from the Europe PMC core abstract.
2. **Our claim:** fig6: "cortical Agt mRNA falls after ischaemia (Kontogiannis 1998, abstract only)".
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "After bilateral renal pedicle occlusion for 60 min, serum creatinine increased, peaking at 72 h"
   - (Abstract) "Both renal cortical angiotensinogen mRNA and proximal tubular AT1 receptor mRNA were significantly reduced early after reperfusion, returning to sham levels by 120 and 72 h, respectively."
   - **Corrections.** Species is **rat**, not mouse ("the effect of renal ischemia/reperfusion on components of the intrarenal renin-angiotensin system in rats"). Ischemia is 60 min bilateral, severe by our standards; time points span 24-120 h, so "early after reperfusion" is probably 24 h, overlapping our 12 h-2 d window. The measurement is whole renal cortex, with no segmental resolution, and the mRNA method is not stated in the abstract - mark the method unverified. Note also that intrarenal angiotensin II rises while Agt mRNA falls.
4. **Axial / injury wording we can borrow:** "binding of 125I-labeled ANG II was preserved in glomeruli but was decreased in whole cortex and outer medulla early after reperfusion" - cortex versus outer medulla as injury zones.

#### Yan 2022 - PPARalpha signaling in SHP-1 insufficiency and renal IRI (PMID 35646989)

1. **Access:** Europe PMC JATS full text (`/PMC9134314/fullTextXML`), `txt/PMC9134314.txt`. Full text Y.
2. **Our claim:** fig6: "Cyp4a10 and Cyp4a14 flattened within S2 | PPARalpha pathway collapse after IRI (Yan 2022, a mutant-mouse study)"; notable result 3 lists "Cyp4a14 and Cyp4a10 (the PPARalpha programme) are flattened".
3. **Verdict: PARTLY.** Both genes are named, but the contrast is genotype within IRI, not IRI versus sham.
   - (Results, KEGG) "the most significantly enriched (FDR < 0.005) KEGG pathways for the 62 genes that were significantly down-regulated more than twofold were the peroxisome-proliferator-activated receptor alpha (PPARalpha) signaling pathway and retinol metabolism"
   - (Results, same paragraph) "we identified the nine PPARalpha signaling pathway-related genes CD36, Fabp1, Acaa1b, Ehhadh, Cyp4a10, Cyp4a14, Cyp4a31, Pck1, and Hmgcs2"
   - (Results) "We also compared the expression of these genes in the control kidneys without IRI in Ptpn6+/mev and wild-type mice and found that there were no significant changes"
   - **Corrections.** The downregulation is in Ptpn6+/mev (SHP-1-insufficient heterozygotes) **relative to wild-type littermates, both 24 h after IRI**. The paper never reports a wild-type IRI versus wild-type sham comparison for these genes, so it does not establish that IRI alone flattens the Cyp4a programme; it establishes that losing SHP-1 deepens the loss, and that the PPARalpha agonist fenofibrate restores Acaa1b, Ehhadh, Cyp4a10 and Cyp4a14 and lessens injury. Model: "three 7-week-old Ptpn6+/mev mice and three wild-type littermates" on C57BL/6, right nephrectomy then "The pedicle of the left kidney was clamped for 34 min", "At 24 h post-IRI, the mice were sacrificed"; n = 3 per group, whole-kidney bulk RNA-seq plus RT-PCR. **Sex is not reported anywhere** - flag it. No segmental or spatial data.
4. **Axial / injury wording we can borrow:**
   - "proximal TECs use fatty acids as their main energy source due to the high energy demand"
   - "cell cycle function-associated genes are generally upregulated at the peak of injury (24 h after IRI)" - matches our replication-programme timing.

#### Naito 2009 - Renal ischemia-induced cholesterol loading along the HMG CoA reductase gene (PMID 19095962)

1. **Access:** PMC article HTML with browser user agent, re-converted to `Dd/PMC2631318.txt`. Full text Y - **fig6's "abstract only" tag is wrong and should be removed.**
2. **Our claim:** fig6: "Cyp51 flattened within S2 | Ischaemia induces cholesterol synthesis (Naito 2009, abstract only), but Liu 2017 gives -2.2 at 48 h | Phase-dependent conflict; state it".
3. **Verdict: PARTLY.** Ischemia does activate the cholesterol-synthesis pathway, but the paper's gene is Hmgcr, not Cyp51, and the mRNA effect is small.
   - (Abstract) "Ischemia/reperfusion (I/R) induced marked cholesterol loading, which corresponded with elevated Pol II recruitment to HMGCR and increased expression levels of both HMGCR protein and mRNA."
   - (Results) "at both 4 and 72 hours after ischemia, HMGCR mRNA levels were elevated by ~30 to 35%"
   - (Results) "at 3 days after unilateral renal ischemia ... cholesterol levels were elevated by ~35%" with "an approximate sixfold increase in cholesterol ester content"
   - (Methods, Ischemia-Reperfusion Injury) "Male CD 1 mice (30 to 35 g; Charles River Laboratories ...)"; "subjected to 30 minutes of left renal pedicle occlusion"; "Approximately 72 hours after surgery, they were re-anesthetized, and bilateral nephrectomy was performed."
   - **Corrections.** Species/strain/sex: male CD-1 mice (not C57BL/6). Model: 30-min **unilateral** left ischemia, contralateral kidney as internal control (validated against sham animals); readouts at 4 h and 72 h; tissue is renal cortex. **Cyp51 is never mentioned**, and no sterol-pathway gene other than Hmgcr is measured, so this reference supports only "ischemia activates cholesterol synthesis", at 4-72 h, by about 30-35% at the mRNA level. The "phase-dependent conflict" with Liu 2017 (Cyp51 -2.2 at 48 h) therefore compares two different genes as well as two phases; say so.
4. **Axial / injury wording we can borrow:**
   - (Fig. 1 legend) "Extensive proximal tubular necrosis and cast formation was observed in the renal cortex (A) and in the outer medullary stripe (B)."
   - "a failure of feedback inhibition of mevalonate/cholesterol signaling pathway(s) clearly existed" - the mechanism for a sterol gene staying up under cholesterol excess.

#### Balzer 2022 - Single-cell analysis of adaptive versus fibrotic kidney regeneration (PMID 35821371)

1. **Access:** Europe PMC JATS full text (`/PMC9276703/fullTextXML`), `txt/PMC9276703.txt`, plus the Europe PMC supplementary-file bundle (`Dd/balzer_supp/`, Supplementary Data 2 = `41467_2022_31772_MOESM5_ESM.xlsx`). Full text Y.
2. **Our claim:** fig6: "MCT1 is basolateral and S1-high; S3 is normally near-silent (0.2 TPM); bulk induction about 3-5x (Liu 2017; Balzer 2022)"; earlier review "Balzer 2022 shows ~3x at day 1".
3. **Verdict: PARTLY - the number exists only in supplementary data.**
   - Slc16a1/MCT1/monocarboxylate appears **nowhere** in the main text, figures or figure legends; the only lactate-related mention is an LDH cytotoxicity assay.
   - (Supplementary Data 2, bulk RNA-seq TPM, row Slc16a1) controls 8.03, 10.88, 9.20, 4.31, 5.57, 5.73, 3.64, 5.56 (mean 6.6); IRI short 1 d 29.41 and 11.07 (mean 20.2, 3.1x); IRI long 1 d 26.37 and 4.62 (mean 15.5, 2.3x); day 3 7.4 and 10.8 (1.1-1.6x); day 14 3.7 and 7.1 (0.6-1.1x).
   - (Methods) "Gene and isoform expression levels (TPM) were estimated using RSEM."
   - (Fig. 1a legend) "Male C57BL/6 mice were subjected to short and long ischemia and followed for 1, 3 and 14 d (n = 12). Bulk and scRNA-seq was performed in n = 2 mice per condition and time point and compared to n = 6 controls"
   - (Results) "We analyzed changes following short (23 min) or long (30 min) bilateral ischemia and collected samples 1, 3, and 14 d post-ischemia"
   - **Corrections.** Cite it as "our calculation from Balzer 2022 Supplementary Data 2", with n = 2 per injured group, a 2.7-fold spread between the two day-1 replicates, and 8 controls. The induction is transient (day 1 only) and whole-kidney, so it cannot speak to which segment carries it. The paper's own PT labels are segmental only in the coarse sense: "PT S1 proximal convoluted tubule; PT S3, proximal straight tubule". Note that the usable time point (day 1, 23-30 min bilateral) is a good match for our 12 h-2 d IRI.
4. **Axial / injury wording we can borrow:**
   - "PT S1 proximal convoluted tubule; PT S3, proximal straight tubule" (cluster nomenclature)
   - "a model of adaptive and maladaptive repair by carefully titrating renal ischemia time" - useful framing for injury dose.

#### Yu 2025 - NPM1 lactylation and ferroptotic trigger waves in graft kidney (PMID 41350284)

1. **Access:** Europe PMC JATS full text (`/PMC12783763/fullTextXML`), `txt/PMC12783763.txt`. Full text Y.
2. **Our claim:** fig6: "MCT1 lactate entry drives ferroptosis (Yu 2025)".
3. **Verdict: PARTLY - correct in direction, but MCT1 is the weaker half of the pair.**
   - (Results) "NC-CM supplemented with sodium lactate recaptured the function of IR-CM to induce ferroptosis in normal HK-2 cells; this effect was abolished by AZD3965, a selective monocarboxylate transporter 1 (MCT1) inhibitor blocking lactate import"
   - (Results) "In vivo validation further revealed that inhibitors of MCT4 and MCT1 suppressed ferroptosis ... and protected renal function"
   - (Results) "lactate transport rate ... was decreased, following MCT4 inhibition in lactate-secreting cells or MCT1 inhibition in lactate-importing cells"
   - (Methods, Animals) "Six-week-old male C57BL/6 mice"; (Methods, In situ kidney IRI model) "the right renal pedicle was occluded for 35 min using a microvascular clamp and released to induce reperfusion injury"
   - **Corrections.** The transporter that is knocked out is **MCT4 (Slc16a3)**, tubule-specific via AAV9-KSP1.3 and Slc16a3CKO; MCT1 evidence is pharmacological (AZD3965, syrosingopine for MCT4) and siRNA in the human HK-2 line, plus inhibitors in vivo. Primary in vivo system is **allograft kidney transplantation** in male C57BL/6 mice analysed on postoperative day 2, with an in-situ warm IRI model (35 min, after left nephrectomy) alongside; human data are urinary and serum lactate in 150 + 64 transplant recipients, not tissue. Mechanism: lactate taken up by neighbouring cells drives AARS1-mediated NPM1 K257 lactylation, which suppresses SLC7A11 and triggers ferroptosis. **No segmental (S1/S2/S3) localization of MCT1 anywhere**, so this paper makes MCT1 induction in injured S3 plausible as an uptake route but says nothing about where along the PT it happens.
4. **Axial / injury wording we can borrow:**
   - "highly glycolytic IR-RC secreted lactate, triggering ferroptotic waves in adjacent normal RC populations"
   - "ferroptotic trigger waves" / "ferroptosis propagation" - a neighbour-to-neighbour spread language that fits a positional readout.

#### Cai 2025 - SLC39A8-mediated zinc dyshomeostasis potentiates kidney disease (PMID 40956895)

1. **Access:** Europe PMC JATS full text (`/PMC12478044/fullTextXML`), `txt/PMC12478044.txt`, plus the SI appendix from the Europe PMC supplementary bundle (`Dd/cai_supp/pnas.2426352122.sapp.pdf`, text in `sapp.txt`, Fig. S3 page rendered to PNG). Full text Y.
2. **Our claim:** fig6: "ZIP5 stays on further into S2 | ZIP8 and ZIP14 are S3 metal transporters, ZIP8 knockout protects in IRI (Cai 2025)".
3. **Verdict: PARTLY.**
   - (SI Fig. S3 legend) "SLC39A8 is predominantly localized in the renal proximal tubular cells and is conserved in human and mice." Panels A and B are reanalyses of two published mouse scRNA datasets (the panel B dot plot has rows PT(S1-S2) and PT(S3), with the large red dot on PT(S3)); panel C is Human Protein Atlas IHC.
   - (Results) "MUT mice showed significantly attenuated AKI following ischemia-reperfusion (I/R) injury (Fig. 1 C-K)"
   - (Results) "we generated kidney-specific Slc39a8 knockout (KKO) mice by crossing Slc39a8flox/flox mice with Six2-Cre mice"
   - (SI Methods, renal I/R model) "the renal pedicles were clamped for 40 minutes to induce transient ischemia. The mice were sacrificed 24 h after sham or I/R surgery"
   - **Corrections.** (a) **ZIP14 / Slc39a14 is never mentioned** - drop it from this citation and attribute it to Fujishiro 2012. (b) The S3 point is a reanalysis of other people's single-cell data plus HPA IHC, in the SI, not an experiment in this paper; the main text says only "predominantly localized in the renal proximal tubular cells". (c) The protective genotypes are the global A391T loss-of-function knock-in and the kidney-specific (Six2-Cre) Slc39a8 knockout; the opposite direction uses a Kap-Cre Slc30a1 (zinc exporter) knockout. (d) Injury models: I/R (left nephrectomy then 40-min right pedicle clamp, 24 h), glycerol rhabdomyolysis, UUO (14 d) and 5/6 nephrectomy - no cisplatin. (e) Strain C57BL/6J; "Male C57BL/6J mice were purchased from GemPharmatech" covers the diet experiments and the WT arm, but **sex is not stated for the knock-in or knockout cohorts**.
4. **Axial / injury wording we can borrow:** the dataset labels "PT(S1-S2)" and "PT(S3)" (SI Fig. S3B) are the only axial wording; nothing usable in the main text.

#### Fujishiro 2012 - ZIP8, ZIP14 and DMT1 in cadmium and manganese transport in mouse kidney PT cells (PMID 22534978)

1. **Access:** INACCESSIBLE - abstract only. RSC HTML and the OUP PDF both return HTTP 403 to curl and to WebFetch; not in PMC; Semantic Scholar openAccessPdf CLOSED. Judged from the Europe PMC core abstract plus the same group's open review (Biol Pharm Bull 42:1076, 2019, J-STAGE HTML, saved as `Dd/bpb2019.txt`), which restates the 2012 experiments.
2. **Our claim:** fig6 attributes "ZIP8 and ZIP14 are S3 metal transporters" to Cai 2025; Fujishiro 2012 is the actual source.
3. **Verdict: SUPPORTED (abstract).**
   - (Abstract) "In situ hybridization of these transporters revealed that ZIP8 and ZIP14 are highly expressed in the proximal tubules of the outer stripe of the outer medulla."
   - (Abstract) "These results suggest that ZIP8 and ZIP14 expressed in the S3 segment of proximal tubules play significant roles in the absorption of Cd(2+) and Mn(2+) in the kidney."
   - (2019 review, Roles of ZIP8 in kidney) "mRNA of ZIP8 is highly detected in the boundary region between the cortex and outer medulla in the proximal tubules, where the S3 segment (straight parts before the Henle loops) of the proximal tubules are concentrated"
   - **Corrections.** The localization is mouse kidney mRNA by ISH, so S3 is inferred from the outer-stripe position rather than from segment markers. The transport experiments use an **immortalized mouse PTEC line in a trans-well system, not a segment-specific line**: the 2019 review states "a study with immortalized PTECs derived from each segment of S1, S2, and S3 is now in progress at our laboratory", i.e. in 2012 no S1/S2/S3-derived lines existed. Mouse strain, sex and age are not retrievable from the abstract or the review. The claim is about cadmium and manganese handling, not zinc gradients.
4. **Axial / injury wording we can borrow:** "the S3 segment (straight parts before the Henle loops)"; "endocytic reabsorption of low-molecular-weight proteins ... is carried out mainly in the S1 and S2 segments, while the S3 segment is rich in transporters".

#### Melo Ferreira 2021 - Spatial and single-cell integration in kidney injury (PMID 34003797)

1. **Access:** Europe PMC JATS full text (`/PMC8262485/fullTextXML`), `txt/PMC8262485.txt`. Full text Y.
2. **Our claim:** fig6: "Fgb in cortical S1/S2 is known (Melo Ferreira 2021)", used to argue that our Fgb S1 signal is a known pattern (fig6 then drops Fgb as control noise anyway).
3. **Verdict: PARTLY.** The paper's Fgb statement is about an injured PT cluster, not about healthy cortical expression.
   - (Results, mapping of single-cell clusters) "Beyond the traditional epithelial cell types, we identified 2 clusters of injured PT cells highly expressing either Havcr1 or fibrinogen (Fga, Fgb, Fgg)."
   - (same paragraph) "The PT fibrinogen cluster also had some expression of Havcr1, but was rich in Fga, Fgb, and Fgg and markers of the more cortical S1 and S2 PT markers."
   - (Methods, Murine models) "From age-matched 8- to 10-week-old 129/SvEv mice (Taconic Biosciences) ... both renal pedicles were exposed and clamped for 22 minutes or 30 minutes through a midline incision and then released"; "Kidneys were excised upon mouse euthanization 6 hours after each procedure"
   - **Corrections.** Strain is **129/SvEv**, not C57BL/6; bilateral IRI, 22 min for the spatial data and 30 min for the scRNA-seq, both read at **6 h** - earlier than our 12 h-2 d. Mouse sex is not stated (the human reference nephrectomy is a 59-year-old female). The second model is cecal ligation and puncture, not a PT injury model. So the paper supports "fibrinogen genes mark an injured PT state that carries S1/S2 markers at 6 h", not "Fgb is normally expressed in cortical S1/S2". Our own conclusion that control Fgb is at the detection floor is consistent with that.
4. **Axial / injury wording we can borrow:**
   - "markers of the more cortical S1 and S2 PT markers" / "the S3 PT, suggesting injury to the outer stripe region of the medulla"
   - "Both a cortical PT S3 cluster (PT-S3-C) and a medullary PT S3 cluster (PT-S3-OS) were identified" - precedent for splitting S3 by position, useful when we say injury differs along S3.

#### Optional extras behind "MCT1 is basolateral and S1-high"

**Becker 2010, Localization of members of MCT monocarboxylate transporter family Slc16 in the kidney (PMID 20392797).** INACCESSIBLE - abstract only (journals.physiology.org 403 to curl and WebFetch, not in PMC, no repository copy found). Abstract: "MCT1, MCT2, MCT7, and MCT8 localized to basolateral membranes of the epithelial cells lining the nephron. MCT1 and MCT8 were detected in proximal tubule cells". Mouse kidney; supports **basolateral** and proximal-tubule, but gives **no S1/S2/S3 assignment**. It also reports that metabolic acidosis downregulates MCT1 mRNA with unchanged protein - a reminder that MCT1 mRNA moves with acid-base state, which is worth a sentence given that IRI is an acidotic insult.

**Yanase 2008, Cellular expression of SMCT (Slc5a8) and the MCT family in the mouse kidney (PMID 18751721).** INACCESSIBLE - abstract only (Springer, not in PMC, Semantic Scholar CLOSED). Abstract: "MCT1 and MCT2 immunoreactivities were respectively found baso-laterally in S1 and thick ascending limbs of Henle's loop", and "SMCT was present at the brush border in S2 and S3 of proximal tubules". **This, not Becker, is the source for "MCT1 is basolateral and S1-high"** (mouse, immunohistochemistry, protein). It also gives the SMCT1 S2/S3 apical placement that pairs with our lactate-handling story.

---

#### Group D closing note: items to fix in fig6

1. "Agt mRNA is S2/S3 (Matsusaka 2012)" -> Agt mRNA is S3 / proximal straight tubule by ISH; S1/S2 is protein reabsorbed via megalin. Our "normally highest in S2" then rests on our own controls, not on Matsusaka.
2. "Normally S1-restricted (85.9 / 2.4 / 0.0 TPM; Fuchs 2024)" -> the TPM values are Chen 2021; Fuchs supports only "proximal tubules of the outer cortex" (RNAscope) and the calcitriol-driven spread "along all segments".
3. "FGF23 rises in AKI" (Christov) cannot carry a Cyp24a1 step: Christov never measures Cyp24a1, 1,25D is unchanged, and in the same folic-acid model Cyp24a1 falls 17-fold.
4. "PPARalpha pathway collapse after IRI (Yan 2022)" -> the comparison is SHP-1-insufficient versus wild-type, both after IRI; no wild-type IRI-versus-sham result for Cyp4a10/Cyp4a14 in that paper.
5. "Naito 2009, abstract only" -> full text is in PMC; the gene is Hmgcr (mRNA +30-35% at 4 h and 72 h), not Cyp51, in male CD-1 mice after 30-min unilateral ischemia.
6. "ZIP8 and ZIP14 are S3 metal transporters ... (Cai 2025)" -> ZIP14 is not in Cai; the S3 assignment belongs to Fujishiro 2012 (ISH, outer stripe), and Cai's own S3 panel is a reanalysis of published scRNA data in its SI.
7. "bulk induction about 3-5x (... Balzer 2022)" -> Balzer never mentions Slc16a1; the ratio (about 3.1x at day 1, short IRI; 2.3x long IRI; n = 2) is our own calculation from Supplementary Data 2 and should be labelled as such.
8. "MCT1 lactate entry drives ferroptosis (Yu 2025)" -> the knockout is MCT4 (Slc16a3); MCT1 evidence is inhibitor and siRNA, in human HK-2 cells and in a transplant model read at day 2.
9. "Fgb in cortical S1/S2 is known (Melo Ferreira 2021)" -> that cluster is an injured PT state at 6 h in 129/SvEv mice; the paper does not report basal Fgb in healthy S1/S2.
10. Add Yanase 2008 as the citation for "MCT1 is basolateral and S1-high" (mouse IHC); Becker 2010 gives basolateral PT but no segment.

### Group E — Second-pass retrieval of papers the first pass could not open

#### Vrhovac 2015 — SGLT1/SGLT2 localization in human kidney (PMID 25304002)
1. **Access:** INACCESSIBLE — abstract only. Routes tried this round: Wayback `link.springer.com/article/10.1007/s00424-014-1619-7` (captures 2015, 2020, 2022, **2023-07-27 fetched**) — the archived page carries the full abstract, the reference list and 11 figure *thumbnails* with no captions and no Results/Methods; `link.springer.com/article/.../figures/1` and `/figures/5` (Springer "Client Challenge" to curl, 303 → `idp.springer.com/authorize` to WebFetch); Europe PMC core record (no PMCID, "Subscription required" only); FULIR EPrints search (`fulir.irb.hr/cgi/search?q=...`, no record — that repository belongs to a different Zagreb institute); `repozitorij.imi.hr` (Dabar bot check, HTTP 418); Internet Archive Scholar (JS session challenge); title+pdf web search (paywalled Springer plus one unverified mirror, not used). Full text **N**.
2. **Our claim:** fig5 — "human SGLT2 protein covers S1 and S2, SGLT1 is the S3 transporter, Vrhovac 2015". fig5lit — "Also protein: hSGLT2 occupies S1 *and* S2 brush border while hSGLT1 is the S3 transporter (ABSTRACT); mouse SGLT1 is already high in S2."
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "hSGLT2 and SGLT1 were localized to the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively."
   - (Abstract) "The renal locations of both transporters largely resembled those in rats and mice"; "the expression of both hSGLTs was sex-independent."
   - (Abstract) "Different to rodents, the renal expression of hSGLT1 was absent in thick ascending limb of Henle (TALH) and macula densa."
   **Correction / limits:** the four things the task asked to check — tissue origin (nephrectomy vs donor), number of kidneys, sex counts, and how S1/S2/S3 were identified in human sections — are all in the unreachable Results/Methods; our documents must not state them. Likewise the abstract never says whether hSGLT1 is weakly present in S1/S2, so "SGLT1 is the S3 transporter" should be written as "localized to S3" rather than "absent earlier". The clause "mouse SGLT1 is already high in S2" needs a rodent source (the paper's only rodent sentence is the "largely resembled" one), or should be labelled as our own mouse measurement.
4. **Axial wording we can borrow:** "the brush border membrane (BBM) of proximal tubule S1/S2 and S3 segments, respectively"; "the renal locations of both transporters largely resembled those in rats and mice".

#### Breljak 2016 — NaDC3 and OAT1-3 along the human nephron (PMID 27053689)
1. **Access:** INACCESSIBLE — abstract only. Routes tried this round: `archive.org/wayback/available` → one capture of `ajprenal.physiology.org/content/311/1/F227` (2016-07-14), fetched and converted: it is the publisher landing page with the complete abstract plus "View Full Text / Please sign in"; the Wayback copies of `F227.full` and `F227.full.pdf` return 404; CDX over the prefix `ajprenal.physiology.org/content/311/1/` lists ~60 captures for that issue but only the bare `F227` URL (other articles in the issue do have `.full`/`.full.txt` captures, so this is luck of the crawl); CDX for `journals.physiology.org/doi/full|pdf/10.1152/ajprenal.00113.2016` → no captures (only a 302 of the bare DOI path); direct `physiology.org/doi/pdf/...` 403 to curl with browser headers and to WebFetch, despite Europe PMC labelling it "Free pdf (Unpaywall)"; Göttingen DSpace-CRIS record is metadata only. Full text **N**.
2. **Our claim:** fig5 — "OAT2 protein is S1/S2-dominant in human, so OCTN2 → OAT2 is unresolved". fig5lit — "human OAT2 protein sits in S1/S2-dominant cortical PT with the outer-stripe S3 unstained".
3. **Verdict: SUPPORTED (abstract).** The task's caution — whether OAT2 is actually stained — resolves in our favour: OAT2 is one of the three OATs and appears by name in both the localization and the variability sentences.
   - (Abstract) "OAT1–OAT3-related immunostaining in human kidneys was detected only in the BLM of cortical proximal tubules; all three OATs were stained more intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment in the outer stripe remained unstained."
   - (Abstract) "Expression of NaDC3, OAT1, OAT2, and OAT3 proteins exhibited considerable interindividual variability in both male and female kidneys, and sex differences in their expression could not be detected."
   - (Abstract) "The specificity of antibodies was validated with human embryonic kidney-293 cells stably transfected with functional OATs."
   **To add when quoting:** the signal is **basolateral**, so it cannot be set against an apical transporter without saying so; both sexes were used and no sex difference was found; the number of kidneys and how S1/S2/S3 were told apart are not recoverable, so write "cortical proximal tubules, S1/S2 > S3 in medullary rays, outer-stripe S3 unstained" and stop there.
   **Independent corroboration we do hold in full text:** Ljubojević 2007 (full text, this batch) states "the transporter in the human kidney (hOAT2) was immunolocalized to the basolateral membrane of cortical proximal convoluted tubules", and Shen 2015 (full text, this batch) finds human OAT2 basolateral **and** apical — so human OAT2 membrane polarity is itself contested, which strengthens the doc's "unresolved".
4. **Axial wording we can borrow:** "stained more intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment in the outer stripe remained unstained"; "only in the BLM of cortical proximal tubules".

#### Edwards 1999 — neutral endopeptidase along the rat and rabbit nephron (PMID 10352425)
1. **Access:** INACCESSIBLE — abstract only (Europe PMC core). Routes tried this round: CDX for `www.karger.com/Article/FullText/28304` and `/Article/Pdf/28304` (captures exist but all are 301 redirects), CDX for the legacy `content.karger.com/ProdukteDB/produkte.asp?Aktion=ShowFulltext&ArtikelNr=28304` (collapses to a 301 of the bare script), karger.com 403 to curl/WebFetch, exact-title search (aggregator abstract pages only). Full text **N**.
2. **Our claim:** fig5 — "Partial precedents: NEP polarity flips between rat and rabbit microdissection (Edwards 1999, abstract)", supporting notable result 1 (neprilysin early in human, late in mouse).
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "In the rat, only the proximal straight tubule and glomerulus had measurable NEP activity of 86 +/- 11.3 pmol/min/mm tubule length and 5.8 +/- 1.5 pmol/min/glomerulus, respectively."
   - (Abstract) "In the rabbit, significant activity was observed in both the proximal convoluted tubule (70.8 +/- 7.2 pmol/min/mm) and proximal straight tubule (29.6 +/- 2.3 pmol/min/mm)".
   - (Abstract) "The results demonstrate that in both the rat and rabbit NEP is restricted to the glomerulus and proximal tubule."
   **Correction:** "polarity" is wrong — the readout is **activity in microdissected segments** and the flip is axial (rat: straight only; rabbit: convoluted > straight, ~2.4×). Glomerular activity is measurable but small and in different units, so it cannot be compared with the tubular numbers. No human and no mouse data, so this is a precedent for a species-dependent NEP axis, not direct support for "early in human, late in mouse". Suggested wording: "the direction of the NEP gradient within the proximal tubule already differs between rodent species — activity is confined to the straight tubule in rat but highest in the convoluted tubule in rabbit (Edwards 1999, abstract)". The segment list, the number of animals and whether a fluorogenic substrate was used are not in the abstract; the inhibitor data (phosphoramidon, thiorphan IC50) and the phorbol-ester result are.
4. **Axial wording we can borrow:** "only the proximal straight tubule and glomerulus had measurable NEP activity"; "NEP is restricted to the glomerulus and proximal tubule".

#### Ljubojević 2007 — renal OAT2 in rats and mice is regulated by sex hormones (PMID 16885152)
1. **Access:** **Wayback of the old HighWire full text**, `http://web.archive.org/web/20100701122043/http://ajprenal.physiology.org:80/cgi/content/full/292/1/F361` (found via the availability API after the modern `/content/292/1/F361` path returned nothing). Full text **Y**. Saved `dl_D/E/ljubojevic2007.{html,txt}`.
2. **Our claim:** fig5lit — "Rodent OAT2 is S3 brush border, female > male and strongly androgen-inhibited (ABSTRACT), so the male mouse under-reports OAT2 at S3."
3. **Verdict: PARTLY** — right for rat; the mouse half needs restating, and the mRNA-vs-protein distinction matters for us.
   - (Results, mouse) "In the M mouse kidney, cortical tubules and other structures exhibited no significant staining ... whereas the brush-border of proximal tubule S3 segments in the outer stripe was weakly stained."
   - (Fig. 8 legend) "Densitometric evaluation of the 66-kDa band (F) exhibited 10-fold higher abundance of this protein in F animals (n = 4 in each experimental group, P < 0.001)."
   - (Discussion) "The reason for a discrepancy in the expression of renal OAT2 mRNA and protein in mice is unknown and indicates possible sex hormone-dependent posttranscriptional regulation."
   **Corrections:** (i) the hormone experiments (castration, ovariectomy, testosterone/estradiol/progesterone) are **rat only**, Wistar, prepubertal 20–25 d / young 6 wk / adult 11–12 wk; (ii) the mouse data are protein only — IHC plus whole-kidney BBM blots, with **no mouse mRNA measured here**, and the mouse strain and age are never given; (iii) earlier mouse work found **no** sex difference in renal Oat2 mRNA, so on an mRNA readout like ours the "male mouse under-reports OAT2" step is not supported — say "mouse OAT2 *protein* is ~10× lower in males"; (iv) male mouse S3 is weakly positive whereas male rat S3 is essentially negative; (v) in both rodents the rat/mouse localization is apical (BBM) while human OAT2 is basolateral — the authors end by doubting rodents as a model for human OAT2.
4. **Axial / sex wording we can borrow:** "brush border of proximal tubule S3 segments ... causing zonal (outer stripe > cortex) and gender differences (F > M)"; "only the brush border of S3 segments in medullary rays was fragmentary and weakly stained".

#### Karbach 2000 — OCT1 and OCT2 localization in rat kidney (PMID 10997918)
1. **Access:** **Wayback of `ajprenal.physiology.org/cgi/content/full/279/4/F679`** (snapshot 2008-09-07). Full text **Y**. Saved `dl_D/E/karbach2000.{html,txt}`.
2. **Our claim:** fig5lit — "Rodent protein agrees that OCT2 is late (rOCT2 mainly S2/S3)".
3. **Verdict: SUPPORTED**, with the mRNA caveat now confirmed from the Results rather than inferred from the abstract.
   - (Results, rOCT2 protein) "In the medullary rays, IR of ab2 localized rOCT2 to the straight portions of the S2 segments ... In S1 segments of proximal tubules it was difficult to detect IR of ab2 above background."
   - (Results) "our data suggest that rOCT2 is transcribed in all segments of the proximal tubule but that rOCT2 protein is mainly expressed in the S2 and S3 segments."
   - (Results, rOCT2 mRNA) "The hybridization signal of rOCT2 in the proximal tubules was less intense in the cortex than in the outer stripe."
   **Notes for accuracy:** protein is **basolateral** in both transporters; the axial contrast is protein-level — mRNA is present in S1, S2 and S3 (only cortex-vs-outer-stripe graded), which matters because our own measurement is mRNA. rOCT1 is the opposite case and is itself discordant: protein S1 and S2 strong, S3 weaker, while its **mRNA** is strongest in the pars recta ("consistent with a localization in late S2 and S3 segments"). Animals: Sprague-Dawley rats, 300–350 g, **sex never stated**, and no sex comparison — the androgen half of our sentence rests entirely on Urakami 2000. The paper also rejects an earlier report that assigned rOCT2 mRNA to S3 exclusively, on probe-specificity grounds.
4. **Axial wording we can borrow:** "rOCT2 protein is mainly expressed in the S2 and S3 segments"; "localized rOCT2 to the straight portions of the S2 segments".

#### Levillain & Hus-Citharel 1998 — ornithine decarboxylase along the mouse and rat nephron (PMID 9841492)
1. **Access:** **Wayback of `ajprenal.physiology.org/cgi/content/full/274/6/F1020`** (snapshot 2010-04-25; eight captures exist 2003–2010). Full text **Y**. Saved `dl_D/E/levillain1998.{html,txt}`.
2. **Our claim:** fig5lit, Odc1 row — "activity PCT > PST".
3. **Verdict: SUPPORTED** in both species, with one caveat that cuts against an mRNA readout.
   - (Results, mouse permeabilized tubules) "Under these conditions, Orn was exclusively decarboxylated in the proximal tubule: PCT > CPST OSPST (Table 1)." (the source renders "≈" as a space)
   - (Results, rat) "The highest rate of Orn decarboxylation occurred in the PCT, followed by the CPST and the OSPST."
   - (Discussion) "But the highest level of ODC mRNA was found in the outer stripe of the outer medulla of the rat kidney (43), more precisely, in the OSPST (S3)."
   **Method / sex:** male OF-1 Swiss mice (35–40 g) and male Sprague-Dawley rats (150 g); collagenase-digested kidney, hand microdissection, L-[1-14C]Orn with DFMO as the ODC-specific control, plus Triton-permeabilized tubules for ODC activity proper. Segments are named **PCT, CPST (= S2), OSPST (= S3)** — S1 and S2 are not separated, so the paper cannot place a within-convoluted gradient (the authors say a progressive decline along the PCT "seems highly likely" but "could not be demonstrated"). Rat OSPST rests on one rat (n = 10 samples) and the authors call it provisional. ODC activity is ~2× higher in every rat segment than the mouse equivalent; rat but not mouse MTAL and medullary collecting ducts have activity.
   **Caveat to carry into the Odc1 row:** activity and mRNA disagree in rat (activity PCT-high, mRNA outer-stripe-high) while in mouse ODC mRNA is PCT ("the ODC gene was highly expressed in the PCT ... the nephron segment which contained the highest ODC activity"). Our measurement is mRNA, so cite this paper for mouse and flag the rat mismatch.
4. **Axial wording we can borrow:** "ODC activity is located in the proximal tubule, with a decreasing intensity from the PCT to the CPST"; "the intensity of the reaction was high in the proximal tubule ... and seemed to decrease from the cortex to the medulla".

#### Kontogiannis & Burns 1998 — AT1 receptors in renal ischemic injury (PMID 9458826)
1. **Access:** **Wayback of `ajprenal.physiology.org/cgi/content/full/274/1/F79`** (snapshot 2007-05-14; ten captures 2003–2007). Full text **Y**. Saved `dl_D/E/kontogiannis1998.{html,txt}`. **fig6's "abstract only" label should be removed.**
2. **Our claim:** fig6 — "cortical Agt mRNA falls after ischaemia (Kontogiannis 1998, abstract only)", supporting the notable result that Agt is suppressed.
3. **Verdict: SUPPORTED**, now with the numbers and the method.
   - (Methods, Animal model and surgery) "Male Sprague-Dawley rats (250-350 g) were used for all studies ... To induce renal ischemia, both renal pedicles were occluded for 60 min with vascular clamps."
   - (Results) "in ischemic renal cortex at 6, 24, and 72 h after reperfusion, angiotensinogen mRNA was markedly decreased [angiotensinogen mRNA intensity (ischemic/sham) at 6 h: 0.42 ± 0.14 ... 24 h: 0.01 ± 0.01 ... 72 h: 0.02 ± 0.01]."
   - (Fig. 4 legend) "Renal ischemia/reperfusion decreases renal cortical angiotensinogen mRNA."
   **Precision for our wording:** species **rat**, bilateral 60-min pedicle occlusion, time points 0/1/3/6/24/72/120 h, n = 4–5 per point; angiotensinogen was measured by **Northern blot of whole renal cortical RNA** — cortex only and no segmental resolution (the Percoll-isolated proximal-tubule fraction was used for AT1 and AT2 receptor mRNA, not angiotensinogen). Recovery by 120 h (0.82 ± 0.15, NS). The drop is essentially complete at 24–72 h, which brackets our 12 h–2 d window.
4. **Injury wording we can borrow:** "renal ischemic injury causes an early increase in intrarenal ANG II levels, associated with reduction of mRNA for angiotensinogen and proximal tubular AT1 receptors"; "Renal angiotensinogen gene expression recovered at 120 h after reperfusion".

#### Chan 1979 — two aminotransferases and D-amino acid oxidase within the nephron (PMID 39098)
1. **Access:** INACCESSIBLE — abstract only, but now from the publisher's own archived page: **Wayback of `www.jhc.org/cgi/content/abstract/27/3/751`** (2008-12-05). The archived `/cgi/reprint/27/3/751` is a frameset whose PDF target (`/cgi/reprintframed/...`, `/cgi/reprint/27/3/751.pdf`) is not archived; `journals.sagepub.com/doi/pdf/10.1177/27.3.39098` has only 302/403 captures and 403s live. Full text **N**.
2. **Our claim:** fig5lit, Gpt row — "ALT is the mirror image of AST along the proximal tubule and rises into the straight segment".
3. **Verdict: PARTLY (abstract).**
   - (Abstract) "Alanine aminotransferase (lower by a factor of 100 in absolute terms) was distributed as the mirror image of aspartate aminotransferase within proximal and distal tubules."
   - (Abstract) "Adult aspartate aminotransferase activities (percentage of the highest) ranged from 100 in the distal straight segment to 25 in the late part of the proximal straight segment".
   - (Abstract) "During development alanine aminotransferase increased 6-fold and D-amino acid oxidase, 4.5-fold in proximal straight tubules".
   **Correction:** the "mirror image" half is verbatim; "rises into the straight segment" is **our inference** — the abstract gives no ALT value for any single PT part, and the 6-fold rise in proximal straight tubules is a developmental, not an axial, statement. Keep the magnitude caveat (ALT is 100× smaller than AST). Species rat, "young and adult", sex not stated; the axial unit is "4 parts of the proximal tubule", not S1/S2/S3; the method is quantitative enzyme assay on dissected freeze-dried segments (Lowry microanalysis), not histochemistry.
4. **Axial wording we can borrow:** "4 parts of the proximal tubule, 2 parts of the distal tubule"; "distributed as the mirror image of aspartate aminotransferase within proximal and distal tubules".

#### Hosoyamada 2004 — function and localization of urate transporter 1 in mouse kidney (PMID 14747372)
1. **Access:** **Wayback of `jasn.asnjournals.org/cgi/content/full/15/2/261`** (snapshot 2008-12-03; six captures 2007–2008 — the modern `/content/15/2/261` captures are the paywalled redesign). Full text **Y**. Saved `dl_D/E/hosoyamada2004.{html,txt}`.
2. **Our claim:** fig5lit — "Mouse URAT1 is brush-border and male > female".
3. **Verdict: SUPPORTED**, and the full text adds an axial statement the abstract hides.
   - (Results, IHC) "Staining of mURAT1 was restricted to the proximal tubule and is not present in the medulla (Figure 4C)."
   - (Results, IHC) "Localization of mURAT1 was observed from the exit of Bowman's capsule ... Therefore, it was difficult to indicate the relationship between the localization pattern of mURAT1 and the segments of proximal tubule."
   - (Fig. 6 legend) "The ratio of band intensity of mURAT1 mRNA/-actin mRNA was 1.17 ± 0.05 in male and 0.50 ± 0.10 in female (mean ± SD)."
   **Corrections / additions:** (i) no S1/S2/S3 assignment is claimed, and the authors say so explicitly — but "restricted to the proximal tubule and not present in the medulla", starting at the exit of Bowman's capsule, is a **cortical/early-PT** statement that weighs against an outer-stripe S3 position; the Discussion reads it as "coinciding with physiologic studies of the reabsorption site of urate in early proximal tubule"; (ii) the IHC signal is mixed — "dotlike patterns in the cytoplasm and linear patterns in the brush border membrane" — and only the N-glycosylated 70-kD form is in the BBM fraction; (iii) sex: three 10-wk males vs three 10-wk females, whole-kidney Western and Northern, mRNA male/female ratio 2.3; the IHC used 10-wk **males only**; (iv) mouse strain is not stated (animals bought from a Saitama supplier); (v) functional work is Xenopus oocytes (Km 1213 ± 222 µM).
4. **Axial wording we can borrow:** "restricted to the proximal tubule and is not present in the medulla"; "observed from the exit of Bowman's capsule ... the reabsorption site of urate in early proximal tubule".

#### Motohashi 2013 — OCT, OAT and MATE protein localization in human kidney (PMID 23630107)
1. **Access:** INACCESSIBLE — abstract only. Routes tried this round: `jpharmsci.org/article/S0022354915309333/pdf` (403 to curl with browser headers and to WebFetch, although Europe PMC lists it as a free Unpaywall PDF), CDX for that path and for `onlinelibrary.wiley.com/doi/10.1002/jps.23567/pdf` (no captures). Full text **N**.
2. **Our claim:** fig5lit — "Human OCT2 protein is PT-wide with no segment resolution".
3. **Verdict: SUPPORTED (abstract)**, with one tightening.
   - (Abstract) "the expression and distribution of apical MATE1 and MATE2-K, and basolateral OAT1, OAT3, and OCT2 were compared using serial sections of human kidney cortex."
   - (Abstract) "The results quantitatively indicated that these transporters were localized to proximal tubules in the renal cortex."
   - (Abstract) "MATE1 and MATE2-K protein were detected in proximal epithelial cells in which OCT2 protein was expressed at the basolateral membranes."
   **Correction:** write "cortical proximal tubule, without S1/S2/S3 resolution" rather than "PT-wide" — only cortex was sampled (laser microdissection for mRNA plus serial-section IHC), so outer-stripe S3, exactly where a late human OCT2 would sit, was never looked at. OAT2 is **not** among the transporters in this paper.
4. **Axial wording we can borrow:** "serial sections of human kidney cortex"; "localized to proximal tubules in the renal cortex".

#### Shen 2015 — OAT2 (SLC22A7): creatinine transport and species-dependent tubular expression (PMID 25904762)
1. **Access:** **Wayback of the archived publisher PDF**, `http://web.archive.org/web/20190430052132if_/http://dmd.aspetjournals.org/content/dmd/43/7/984.full.pdf` (five PDF captures 2017–2020), converted with pdftotext. Full text **Y**. Saved `dl_D/E/shen2015.{pdf,txt}`.
2. **Our claim:** fig5lit — "basolateral+apical in human, apical only in rat".
3. **Verdict: SUPPORTED.**
   - (Results, Localization ... by Immunohistochemistry) "In both human and monkey kidneys, OAT2 was located along the basolateral infoldings of the cell membrane, but an apical localization was also visible ... In contrast, rat Oat2 was restricted to the brush-border (apical) membrane."
   - (Methods, Immunohistochemistry) "The kidneys from human, cynomolgus monkey, and Sprague-Dawley rat were removed and fixed in 10% neutral buffered formalin solution."
   - (Discussion) "Oat2 in rat and mouse kidneys was immunolocalized to the brush-border (apical) membrane of the proximal tubule S3 segment in the outer stripe and medullar rays".
   **Resource caveats:** the IHC carries **no segment resolution** (cortex only, labels are just "GL, glomerulus; PT, proximal tubule"), uses a commercial Sigma polyclonal anti-human OAT2 at 1:200, and reports no number of kidneys and no sex/age for the human or monkey tissue — so it cannot place human OAT2 on the axis at all. The functional arm used one cryopreserved male and one fresh female human proximal-tubule cell lot. The Discussion is also a tidy summary of the rat literature conflict (Kojima's rat Oat2 also in TAL and collecting ducts).
4. **Axial wording we can borrow:** "the brush-border (apical) membrane of the proximal tubule S3 segment in the outer stripe and medullar rays"; "a molecular basis for potential species difference in the transporter handling of creatinine".
