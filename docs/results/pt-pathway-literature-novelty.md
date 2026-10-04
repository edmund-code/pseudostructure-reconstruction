# What is new? Literature check of the PT pathway findings

Checked 2026-10-04 against the confident and robust pathway lists of notebooks 31–33
([pt-pathway-method-selection.md](pt-pathway-method-selection.md)). Private reports, dossiers and
tables: `results/pt_pathway_literature/`.

## How the check was done

- **Scope.** All 82 confident `T_spatial` pathways (8 themes) with up to 8 driver genes each, the 70
  genes whose species difference reverses along the PT, and the 55 pathways called only by
  conventional pseudobulk GSEA (species but not relabeled). Each item came with our per-segment
  pattern for both species, DESeq2 fold changes and robustness labels.
- **Two rounds of parallel literature agents** (11 Sonnet agents each). Round 1 used web search
  until a session-wide search cap stopped it. Round 2 used the Europe PMC API (titles and
  abstracts) to fill gaps and re-verify round-1 claims. Round 2 corrected or withdrew roughly
  25 round-1 attributions (misread abstracts, wrong species, one paper that could not be found).
  Only round-2 statuses are used here. I re-checked the key PMIDs myself.
- **Status meanings.** KNOWN: the same species difference and its PT localization are reported.
  PARTLY KNOWN: one part is reported (for example the mouse segment pattern, or a whole-kidney
  species difference). NOT FOUND: nothing after at least two targeted queries; a candidate, not
  proof of novelty. CONTRADICTS: the literature reports the opposite.
- **Limits.** Evidence is abstract-level. Segment-resolved atlases (Ransick 2019, Chen 2023
  supplements, KPMP) were mostly not machine-readable, so "NOT FOUND" means "not in titles and
  abstracts". No search was run with any number or result from our data.

## Bottom line

1. **No published study compares human and mouse PT pathways along the S1→S3 axis.** Every
   pathway-level result is NOT FOUND at the level of a species-by-position comparison. The
   novelty claim is the positional comparison itself, plus a short list of gene programs where
   even the gene-level contrast is unreported.
2. **Our map reproduces a substantial body of known biology.** Examples are male-mouse Cyp2e1 in S1–S2
   with little human renal CYP2E1, the mouse S3 markers Cyp7b1/Slc22a7/Slc5a8/Slc23a1, the
   S1–S2 restriction of mouse B0AT1, the human-only renal UGT1A9 and CYP3A5, and the rodent-only
   apoB. This is a validation section for the paper.
3. **Sex is the dominant confounder for mouse-high late-PT programs.** All specimens are male. Of
   454 driver genes, 58 are mouse-high **and** male-biased in mouse PT (Xiong et al. 2023 tables).
   The literature places the main mouse PT sex differences in S2/S3. Pathways are therefore
   reported with and without sex-biased genes (51 of 61 robust pathways survive their removal).
4. **A handful of "PT" signals are neighbouring-segment or ambient RNA.** These are Egf, Muc1, Klk1,
   Slc26a7, Car12 (human) and Akr1c18, plus hepatocyte-type transcripts in human. Remove them from
   the story.
5. **The conventional-only (average) differences are mostly technical.** They are dominated by
   ubiquitin, proteasome and V-ATPase genes that are uniformly lower in human, including genes with
   equal probe counts. This is consistent with the relabeling result that average-level screens
   are not specific.

## Tier 1 — candidate novel findings with the strongest support

Every row survives re-anchoring of the human coordinate. "Region check" compares human S3 with
**cortical** mouse S3 only (262 tubules no deeper than mouse S1). "Pairings" gives the sign of
human − mouse in S1, S2, S3 for all four human–mouse specimen pairings (+ or − = all four agree).

| Finding | Genes (both species express unless noted) | Our evidence | Literature status |
|---|---|---|---|
| **Creatine synthesis is a zonation shift.** Human AGAT stays high along the whole PT; mouse AGAT is S1-restricted and gone by S3. | Gatm | Pairings +++; S3 difference unchanged in cortical mouse S3; female-biased in mouse, so male mice may understate it. | GATM is a human PT enzyme (Fanconi syndrome, [PMID 29654216](https://pubmed.ncbi.nlm.nih.gov/29654216/)); kidney is the main creatine-synthesis site. Axial species contrast **NOT FOUND** (7 queries). |
| **Mitochondrial β-oxidation is zoned differently.** Mouse rises along the PT; human is flat. Acadm and Acaa2 reverse direction. | Acadm, Hadh, Acaa2, Acads, Acadl | Pairings +−− (Acadm), −++ (Acaa2), −−− (Hadh); not sex-biased in mouse PT; cortical check holds. | **NOT FOUND** after ≥2 queries per gene; ACAA2/ACSM3 appear only in injury papers. |
| **Glutathione synthesis rises toward S3 in mouse only.** | Gclc, Gclm, Gss | Pairings −−; 9–10 core pathways; not sex-biased; cortical check holds. | **PARTLY KNOWN**: rabbit PT GSH synthesis is highest in S3 ([PMID 9612330](https://pubmed.ncbi.nlm.nih.gov/9612330/)). The flat human profile is not reported. |
| **Sterol/mevalonate (SREBP2) program rises toward S3 in mouse only.** | Hmgcs1, Cyp51, Hmgcr, Sqle (human low) | Pairings −−−; cortical check holds; Cyp51 female-biased in mouse. | **NOT FOUND** (3–4 queries per gene). Context: Cyp7b1 loss lowers male kidney sterol synthesis ([PMID 10748048](https://pubmed.ncbi.nlm.nih.gov/10748048/)), so a male S3 axis is plausible. |
| **Lead, not a result: the organic-anion uptake module persists into late PT in human; mouse confines it to S1–S2.** Whole-PT shows no difference for OAT1 and NaDC3. | Slc22a6 (OAT1), Slc22a8 (OAT3), Slc13a3 (NaDC3) | Pairings +−+ (OAT1, NaDC3), +++ (OAT3); segment DESeq2 confirms opposite signs. About half the S3 difference is outer-stripe region; a smaller difference persists in cortical mouse S3. | **PARTLY KNOWN**: human OAT1–3 are stronger in S1/S2 than S3, and human outer-stripe S3 is unstained (Breljak 2016), so region explains much of the S3 contrast. The cortical-only species difference is **NOT FOUND**. Translationally relevant (drug secretion) if it holds in matched cortex. |
| **Human retains B0AT1 into late PT; mouse is S1–S2 only.** | Slc6a19 | Pairings +++; difference unchanged in cortical mouse S3; region-sensitive pathway, so treat as a lead. | Mouse S1–S2 restriction KNOWN ([PMID 35979966](https://pubmed.ncbi.nlm.nih.gov/35979966/)); human axial profile **NOT FOUND**. |
| **Direction reversals of shared genes.** Best: Dcxr (mouse falls, human rises; no mouse sex bias, and male-mouse microdissection data confirm the fall). | Dcxr, Pah, Igfbp4, Glyat, Ugt3a1, Acox2, Nt5e, Cyp24a1 | Segment DESeq2 confirms opposite signs; pairings consistent. Glyat and Cyp24a1 are male-biased in mouse, and the Cndp2 S2 hump is likely androgen-driven. | Mouse PT localization known for several (Igfbp4 [PMID 9324049](https://pubmed.ncbi.nlm.nih.gov/9324049/); rat Nt5e PCT [PMID 19333785](https://pubmed.ncbi.nlm.nih.gov/19333785/)); functional PAH in human kidney known ([PMID 10444341](https://pubmed.ncbi.nlm.nih.gov/10444341/)). Species reversals and axial profiles **NOT FOUND**. |

## Tier 2 — genes expressed by one species only, with zonation in that species

These are species-specific expression rather than differently placed shared genes. The
positional novelty is the zonation description in the expressing species. Each needs per-gene
validation (ISH or protein), because a probe gap would look the same.

| Finding | Literature status |
|---|---|
| **Human PT serine synthesis**: Phgdh/Psat1 high in human S1 and falling, near-absent in mouse PT. | Mouse Phgdh is weak in PT and near-negative in the outer stripe ([PMID 17510490](https://pubmed.ncbi.nlm.nih.gov/17510490/)); human zonation **NOT FOUND** (~12 queries). Both are ATF4/stress-responsive, so donor tissue state matters. |
| **Human UGT1A9, EPHX1, GSTP1 decline from S1 to S3**; mouse Ugt1a9 absent. | Human renal UGT1A9 KNOWN ([PMID 25650382](https://pubmed.ncbi.nlm.nih.gov/25650382/)); axial gradient **NOT FOUND**. The mouse Ugt1a cluster shares exons, so the mouse zero may be a mapping issue. |
| **Human RBP4 rises toward late PT**; mouse near zero. | Only a 1989 rat ISH placing RBP mRNA in outer-stripe S3 ([PMID 2469758](https://pubmed.ncbi.nlm.nih.gov/2469758/)); species contrast **NOT FOUND**. |
| **Human AOX1 rises to S3** (mouse Aox1 near zero); **human AOC1** in PT. | Human PT AOX1 by IHC ([PMID 17992631](https://pubmed.ncbi.nlm.nih.gov/17992631/)); rodents have four AOX genes. AOC1 human PT localization is unresolved. |
| **Mouse-only Inmt (S2 peak), Me1, Lpl, Mogat1, CoA synthesis (Pank1, Coasy)**. | Mouse renal Inmt and Lpl are known; human contrasts **NOT FOUND**. Pank1, Coasy and Mogat1 are male-biased in mouse PT; the CoA cluster has a PPARα hypothesis. Human whole-kidney LPL exists ([PMID 24371263](https://pubmed.ncbi.nlm.nih.gov/24371263/)), so near-zero human PT Lpl needs ISH. |

## Known biology the map reproduces (validation)

| Gene or program | What is known | Reference |
|---|---|---|
| Cyp2e1 mouse S1–S2, little in human kidney | Male-mouse, androgen-induced; no p-nitrophenol oxidase activity in human kidney (one 2023 report finds human PT CYP2E1) | [PMID 7839370](https://pubmed.ncbi.nlm.nih.gov/7839370/), [PMID 38042273](https://pubmed.ncbi.nlm.nih.gov/38042273/) |
| Mouse male S2/S3 programs: Cyp7b1, Slc27a2, Acox1–3, Acsm2/3, Hsd11b1, Hsd17b2 | Main mouse PT sex differences in S2/S3 | [PMID 36758122](https://pubmed.ncbi.nlm.nih.gov/36758122/), [PMID 37673062](https://pubmed.ncbi.nlm.nih.gov/37673062/) |
| Mouse peroxisomal FAO sexual dimorphism | Androgen-linked | [PMID 34651140](https://pubmed.ncbi.nlm.nih.gov/34651140/) |
| Slc22a7 (OAT2) S3 in rodents; human OAT2 on both membranes | Species difference in membrane localization | [PMID 16885152](https://pubmed.ncbi.nlm.nih.gov/16885152/), [PMID 25904762](https://pubmed.ncbi.nlm.nih.gov/25904762/) |
| Slc5a8 and Slc23a1 mostly in mouse S3 | — | [PMID 17692818](https://pubmed.ncbi.nlm.nih.gov/17692818/), [PMID 18614995](https://pubmed.ncbi.nlm.nih.gov/18614995/) |
| Mouse apoB made in kidney; human kidney makes apoE, not apoB | Matches our directions | [PMID 20103594](https://pubmed.ncbi.nlm.nih.gov/20103594/), [PMID 6572003](https://pubmed.ncbi.nlm.nih.gov/6572003/) |
| Adh1 androgen-induced in mouse kidney | Explains mouse-high Adh1 | [PMID 8973327](https://pubmed.ncbi.nlm.nih.gov/8973327/) |
| Gcnt1 (core 2 GlcNAc-T) in mouse S3 | Mouse-specific Gsl5 regulatory element | [PMID 16278214](https://pubmed.ncbi.nlm.nih.gov/16278214/) |
| Pck1 rising S1→S3 | Rat microdissection and human isolated segments agree | [PMID 11939717](https://pubmed.ncbi.nlm.nih.gov/11939717/), [PMID 11716765](https://pubmed.ncbi.nlm.nih.gov/11716765/) |
| Rodent kidney DMET proteins mostly above human | Cross-species proteomics | [PMID 38711199](https://pubmed.ncbi.nlm.nih.gov/38711199/), [PMID 31123036](https://pubmed.ncbi.nlm.nih.gov/31123036/) |

## Contradictions and non-PT signals (remove or caveat)

| Signal | Literature | Likely explanation |
|---|---|---|
| Egf in mouse S3 | TAL/DCT only; PT negative in mouse and human | Neighbouring-segment RNA |
| Muc1 rising in human S3 | Negative in normal human PT | Distal-nephron spill |
| Klk1 (mouse), Slc26a7 | Connecting tubule; intercalated cell | Not PT |
| Car12 high in human PT | Human CA XII weak in PCT | Neighbour RNA or mRNA–protein difference |
| Akr1c18 in male mouse S3 | 20α-HSD mRNA in DCT, female-high | Probe cross-hybridization or boundary |
| Aqp1 flat in human | Human AQP1 highest in straight tubule ([PMID 9013443](https://pubmed.ncbi.nlm.nih.gov/9013443/)) | Saturation (99% detection) or coordinate mismatch |
| Hsd11b1 near zero in human | 11β-HSD1 protein in human PT ([PMID 18261751](https://pubmed.ncbi.nlm.nih.gov/18261751/)) | mRNA level versus protein |
| Alb, Ttr and other hepatocyte-type transcripts in human | Not PT genes | Ambient RNA; phototransduction, sensory perception and phosphorylation sets are driven by these |
| Slc22a7 rising toward S3 in human | Human OAT2 protein is S1/S2 > S3 | mRNA versus protein, or human late PT mis-defined |
| Fgf1 in mouse S3 | Conflicting rat reports (S3 versus absent) | Unresolved |

## Confounders that apply across findings

- **Sex.** All four specimens are male. 1,211 genes are sex-biased in mouse PT (511 male, 700
  female; Xiong et al. 2023 public tables). Removing them keeps 51 of 61 testable robust
  pathways, including all five flagships, but drug-metabolism sets (Aspirin and Paracetamol ADME,
  Drug metabolism) depend on sex-biased genes. Female mice are needed to separate species from sex.
- **Region.** Mouse S3 is mostly outer stripe (83–85%); human sections are cortex. The cortical-S3
  check keeps most late-PT differences (Gatm, FAO, glutathione, sterol, Slc6a19), weakens the
  OAT module by about half, and removes Car2.
- **One-species detection.** 14 of 66 robust pathways depend on genes one species barely detects
  (probe gap or true absence). The Tier 2 genes belong here.
- **"Human retains an S1/S2 program in late PT."** Slc6a19, Slc4a4, Slc9a3 and the OAT module
  share this pattern. It could be biology, or a partly mis-defined human late PT (human S1 spans
  more of the coordinate, and human S3 is cortical). It is one hypothesis, not several findings.
- **Technical offset.** Ubiquitin, proteasome and V-ATPase genes are uniformly lower in human,
  including at equal probe counts. This affects average comparisons, not `T_spatial`.
- **Ortholog families** that are not one-to-one: Ugt1a/2b, Cyp3a, Cyp4a, Akr1c, Aox, Ces1, Adh1.

## Suggested framing for the paper

Lead with the method result (only position-dependent screens are specific), then show the map
reproduces known zonation and species biology. Present the Tier 1 programs as the new positional
findings, each with its replicate, region and sex checks: creatine synthesis zonation,
mitochondrial and peroxisomal FAO, glutathione and sterol synthesis in mouse late PT, and the
OAT module. Present Tier 2 as human PT zonation of species-specific genes (serine synthesis,
UGT1A9, RBP4, AOX1) for validation by ISH. Keep the male-only and cortex-only caveats beside every
late-PT result.
