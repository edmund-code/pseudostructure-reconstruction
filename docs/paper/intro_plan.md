# Introduction plan (Nature Methods Article)

Compiled 2026-11-03 and revised the same day. This is a plan, not prose.
- Outline: GASTON (*Nat. Methods* 2025).
- Substance and tone: the Fan lab's Introductions (median length 477 words).
- Reading notes: `docs/paper/citation_checks/B_intro_refs_22_40.md` and `F_fanlab_intros.md`.

**Target: 500 words in 5 paragraphs, about 35 references.**

The Introduction is methods-first. It opens on a limitation of spatial transcriptomics, then brings in the kidney as
the model example. Each paragraph has one job. Anything a reader does not need in order to accept the gap goes to
Methods, Results or the Discussion.

## P1. The limitation of spatial transcriptomics (~90 words)
- Spatial transcriptomics measures expression in place.
- Many organs are built from repeated structures whose function changes along their length.
- Where such a structure winds out of the section plane, a section cuts each copy many times:
  - each cross-section comes from an unknown position;
  - neighbouring cross-sections belong to different copies;
  - so position in the tissue does not give position along the structure.
- Tracing through serial sections, or profiling in 3D, recovers position but does not scale to whole organs.
- Choose examples that fit the claim: nephrons, and seminiferous tubules. Liver lobules are a weaker example,
  because their axis is in the section plane, marked by the portal and central veins. Use them, if at all, for
  contrast: there, a landmark gives the axis.

## P2. The kidney as the model example (~90 words)
- The nephron runs through ordered segments. The proximal tubule's S1, S2 and S3 differ in transport, metabolism and
  injury susceptibility.
- Function changes within segments as well as between them; for example, glucose uptake passes from SGLT2 to SGLT1.
- A single section holds thousands of PT cross-sections, yet atlases resolve the PT only as three segments.
- The methods point: the kidney is a good test case, because the known segment order gives a check that needs no
  tissue coordinates.

## P3. Prior methods in two categories (~130 words)
- **From the tissue plane:**
  - clustering and spatial domains average within segments;
  - gradient and alignment methods need tissue coordinates, or structure shared in the plane.
- **From expression:**
  - trajectory inference needs a start cell;
  - repeated-structure ordering needs a change in cell-type composition, a known cycle, or whole structures captured
    in 3D.
- **Close:** the PT has one cell type, no landmarks, and partial sampling at unknown positions.

## P4. The concept (~80 words)
- Pseudostructure treats each cross-section as a sample of one shared structure and places it at a position t,
  using only its expression (Fig. 1).
- t carries both segment order and position within a segment.
- As pseudotime orders cells along a process, t orders cross-sections along a structure. It preserves order, not
  length.
- One sentence on the assumption: the copies are exchangeable and are sampled along one unbranched axis.

## P5. The tool and the results (~110 words)
- "Here, we developed pseudostructure reconstruction, a label-free method that …". Name the elastic principal curve,
  the NB gene models and the relabeling calibration.
- Data: 4 mouse sections (2 control, 2 AKI) and 2 human sections.
- Results in this order, method first:
  - recovers segment order without labels, and stays stable under resampling;
  - resolves within-segment genes;
  - as applications: reversals between human and mouse, with human transporter switch points later (**genes**, not
    pathways), and AKI switch points moving upstream.

## Working draft (2026-11-03)

**P1 (agreed by the user).** Spatial transcriptomic (ST) technologies enable investigations into tissue organization
and cellular function by capturing the spatial location and gene expression profiles of cells within a
two-dimensional tissue slice. However, many organs are built from structures, such as kidney nephrons, seminiferous
tubules and the epididymal duct, that wind back and forth through the tissue, with successive parts of each structure
performing different functions. As such, a two-dimensional slice captures each structure only as scattered
cross-sections, losing their positions along the structure. This lack of positional information hinders the
characterization of gene expression variation along repeated structures. Recovering position requires serial-section
tracing or three-dimensional profiling, which have so far captured morphology alone or transcriptomes over
sub-millimetre depths.
- Still to do: find and read in full a source for seminiferous tubules and one for the epididymal duct.

**P2 (draft v4, why the proximal tubule).** In the kidney, each nephron passes through a known sequence of
segments, from the proximal tubule through the loop of Henle and the distal convoluted tubule to the connecting
tubule, which drains into the collecting duct. These segments perform distinct functions, from bulk reabsorption in
the proximal tubule to salt reabsorption in the thick ascending limb and final adjustment of salt, potassium and water
in the distal convoluted tubule, connecting tubule and collecting duct. Of these, the proximal tubule reabsorbs
approximately 60–70% of filtered water and salt and nearly all filtered nutrients. This reabsorption requires high
metabolic activity, which makes the proximal tubule the segment most susceptible to injury. The proximal tubule is
further divided into three segments, S1, S2 and S3, that differ in solute transport and metabolism. For example,
glucose is reabsorbed mainly by SGLT2 in S1 and S2 and by SGLT1 in S3. As such, the known segment order provides an
opportunity to evaluate whether a method recovers position along a structure. However, transcriptomic atlases resolve
the nephron only into these discrete segments, such that how expression changes within each segment remains an open
question.
- The segments match the order used in the Fig. 3 evaluation: PT-S1 < PT-S2 < PT-S3 < ascending limb < DCT < CNT <
  collecting duct.
- The two new sentences justify the proximal-tubule focus by its function and injury, not by "most studied", which no
  source we have read supports. The explicit focus goes in P5 (sentence 5.3: "genes that vary within proximal-tubule
  segments") and the Fig. 4 Results opening.
- Sources (cite later), both read in full:
  - Curthoys & Moe 2014: "reabsorbing approximately 60%–70% of the water and NaCl, a greater proportion of the NaHCO3,
    and nearly all of the nutrients in the ultrafiltrate"; the energy comes from the basolateral Na+/K+-ATPase.
  - Kirita 2020 (PNAS): "proximal tubule (PT) epithelia (the segment most susceptible to injury)"; "this segment
    suffers the most injury due to high metabolic activity".
  - Optional primary source, not yet read: Chevalier 2016 (Am J Physiol Renal Physiol), "The proximal tubule is the
    primary target of injury and progression of kidney disease".
- "Susceptibility to injury" was dropped from the S1–S3 sentence to avoid repeating it. Ransick's "the S3 region is
  particularly vulnerable to injury" can go in the AKI Results instead.
- "Specifically" became "For example" in the SGLT sentence, so the paragraph does not repeat a connective.
- Other sources: Curthoys & Moe for S1–S3 and SGLT2/SGLT1; Chen 2021 and Lake 2023 for the atlases. The functions of
  the distal segments still need a physiology review, read in full.

**P3 (draft).** Attempts at such analysis can be performed using clustering or spatial-domain approaches that
assign each spot or cross-section to a discrete group, assuming that mean expression is constant within each group.
Other approaches model continuous gradients from tissue coordinates or compare slices after aligning shared tissue
structures, and therefore rely on tissue position, which does not reflect position along a nephron. Alternatively,
trajectory inference orders cells along processes from their expression, and related approaches have ordered
structures such as intestinal crypts, seminiferous tubules and hair follicles using changes in cell-type composition,
a known cyclic topology or three-dimensional profiles of whole structures. However, along the proximal tubule, a
single epithelial cell type spans an open, non-cyclic axis, and a cross-section carries no landmark of its position
other than its expression. As such, an approach that orders cross-sections from their expression, without tissue
coordinates or segment labels, provides an alternative strategy for such structures.
- Sources (cite later): BayesSpace, SpaGCN, GASTON (constant-mean assumption); GASTON, CoPro; STalign, STcompare;
  Monocle, DPT, Saelens; SpaceFold (crypts, composition); Chakravorty (testis, cycle); Asami (follicles, 3D).

**P4 (draft).** Thus, we sought to treat cross-sections from the many copies of a structure, such as the nephrons
of a kidney, as samples of a pseudostructure, defined as one shared structure along which each cross-section can be
placed by its expression. The position of each cross-section along the pseudostructure, t, captures both the order of
segments and the position within a segment (Fig. 1). In the context of single-cell data, pseudotime orders cells
along a process from their expression. In the context of tissue structures, a pseudostructure analogously orders
cross-sections along a structure, preserving their order but not physical distance.
- Optional (STdeconvolve's "Briefly, … (i)–(iii)"): Briefly, this ordering leverages (i) the many copies of a
  structure in each slice, (ii) one expression program shared along an unbranched axis, and (iii) slices that
  together span the axis even when each covers it only partially.
- Sources (cite later): Monocle, DPT for pseudotime.

## Sentence plan

Key to the labels:
- **COPY**: a Fan-lab sentence used verbatim, apart from citations, the tool name and the URL.
- **ADAPT**: her sentence word for word, with our nouns swapped in.
- **OURS**: our own sentence, written in her pattern.

Citations are given as keys because the reference numbers will change. Nature runs a similarity check on submissions.
One copied sentence is normal; the adapted ones differ enough to pass. Total: about 510 words (cut 1.4 to reach about 500).

### P1. The limitation of spatial transcriptomics
- **1.1 COPY**, STcompare sentence 1:
  Spatial transcriptomic (ST) technologies enable the investigation of how tissue organization relates to cellular function by capturing the spatial location of cells within a tissue and the gene expression profiles associated with those locations [Ståhl 2016; Moses 2022].
- **1.2 OURS:**
  Many organs are organized as repeated structures, such as nephrons in the kidney and seminiferous tubules in the testis [Zhai 2006; Chakravorty 2026].
- **1.3 OURS**, the "However" pattern:
  However, a two-dimensional section cuts each copy many times, and neighbouring cross-sections often belong to different copies, so tissue position does not reveal position along the structure.
- **1.4 ADAPT**, STdeconvolve ("This lack of single-cell resolution hinders the characterization of cell-type-specific spatial organization and gene expression variation"):
  This lack of positional information hinders the characterization of gene expression variation along repeated structures.
- **1.5 OURS:**
  Recovering position requires serial-section tracing [Zhai 2006] or three-dimensional profiling [Schott 2024], which have so far captured morphology alone or transcriptomes over sub-millimetre depths.

### P2. The kidney as the model example
- **2.1 OURS**, CRAWDAD's "For example" pattern:
  In the kidney, for example, the proximal tubule of each nephron runs through three segments, S1, S2 and S3, that differ in solute transport and metabolism [Curthoys & Moe 2014] and in susceptibility to injury [Ransick 2019].
- **2.2 OURS**, STdeconvolve's "Specifically" pattern:
  Specifically, glucose is reabsorbed mainly by SGLT2 in S1 and S2 and by SGLT1 in S3 [Curthoys & Moe 2014].
- **2.3 ADAPT**, MERINGUE ("Thus, how these subpopulations of cells are organized in space … remains an open question in many systems"):
  Thus, as transcriptomic atlases resolve the proximal tubule into these three segments [Ransick 2019; Chen 2021; Lake 2023], how gene expression varies within segments remains an open question.
- **2.4 OURS**, CRAWDAD's "provide an opportunity" pattern:
  As its segment order is known, the kidney further provides an opportunity to evaluate recovered position.

### P3. Prior methods
- **3.1 ADAPT**, STcompare ("Attempts at such comparative analysis of ST datasets can be performed using traditional … approaches such as …"):
  Attempts at such analysis can be performed using clustering or spatial-domain approaches that assign each spot or cross-section to a discrete group [Zhao 2021; Hu 2021], assuming that mean expression is constant within each group [Chitra 2025].
- **3.2 ADAPT**, STdeconvolve ("As such, these approaches rely on the availability of a suitable single-cell reference"):
  Other approaches rely on tissue coordinates, modelling continuous gradients from tissue position [Chitra 2025; Miao 2026] or comparing sections after aligning shared tissue structures [Clifton 2023; Clifton 2026].
- **3.3 OURS**, MERINGUE's "Alternative approaches" pattern:
  Alternatively, trajectory inference orders cells along processes from their expression [Trapnell 2014; Haghverdi 2016; Saelens 2019], and related approaches have ordered repeated structures by cell-type composition [Niec 2022], a known cycle [Chakravorty 2026] or whole structures profiled in three dimensions [Asami 2026].
- **3.4 OURS**, the "However" pattern:
  However, along the proximal tubule, a single epithelial cell type spans the axis, and no landmark marks where a cross-section lies.
- **3.5 ADAPT**, STdeconvolve ("As such, a reference-free deconvolution approach provides an alternative strategy for deconvolving cell types when an appropriate reference is not available"):
  As such, an approach that orders cross-sections without tissue coordinates or segment labels provides an alternative strategy for such structures.

### P4. The concept
- **4.1 ADAPT**, CRAWDAD ("Thus, we sought to consider …"), with MERINGUE's inline "defined as":
  Thus, we sought to place every cross-section along a pseudostructure, defined as one shared structure of which each cross-section is a sample, at a position t that captures both segment order and position within a segment (Fig. 1).
- **4.2 ADAPT**, STdeconvolve's parallel pair ("In the context of natural language processing, … In the context of ST data, …"):
  In the context of single cells, pseudotime orders cells along a process from their expression [Trapnell 2014; Haghverdi 2016]. In the context of repeated structures, pseudostructure analogously orders cross-sections along a structure, preserving order but not physical length.

### P5. The tool and results
- **5.1 ADAPT**, STalign ("To address these challenges, we present an approach called STalign that builds on recent developments in Large Deformation Diffeomorphic Metric Mapping to align ST datasets using image varifolds"):
  To address these challenges, we present an approach called pseudostructure reconstruction that builds on elastic principal curves [Albergante 2020] to place cross-sections along a shared axis without segment labels or tissue coordinates.
- **5.2 ADAPT**, CRAWDAD ("We demonstrate CRAWDAD on simulated as well as real SRO datasets …, though CRAWDAD is amenable to any SRO data for which cell positions and cell-type annotations can be obtained"):
  We demonstrate pseudostructure reconstruction on Visium HD data [Oliveira 2025] from four mouse and two human kidney sections, though pseudostructure reconstruction is amenable to any ST data in which repeated structures can be segmented.
- **5.3 ADAPT**, MERINGUE ("we show that MERINGUE is able to identify …"):
  We show that pseudostructure reconstruction is able to recover segment order without labels, resolve genes that vary within segments, and identify genes and transporter switch points that differ between species or shift after acute kidney injury.
- **5.4 COPY**, STalign (tool name and URL changed):
  Pseudostructure reconstruction is available as an open-source Python package at [URL].

### Sentence habits these follow
- Every sentence opens with a connective that states its role:
  - a limit: *However*;
  - a consequence: *As such*, *Thus*;
  - an instance: *For example*, *Specifically*;
  - a parallel: *Alternatively*;
  - the turn to the tool: *To address these challenges*.
- One claim per sentence, with the citation at the end of the clause it supports.
- The same term every time: "cross-section", "repeated structure", "pseudostructure reconstruction".
- A new term is defined inline with "defined as".
- Prior methods are described by what they do and what they rely on. The tool's name is the subject of the sentences
  that say what it does.

## Cut from the current draft (787 words)
| Material | Where it goes |
|---|---|
| Visium HD specifications | Methods |
| "distinguished by transporter expression, not ultrastructure" | Results |
| Descriptions of CoPro, SpaceFold, Chakravorty and Asami | Discussion (keep one citation clause each here) |
| The segmentation pipeline | Methods (keep half a clause here) |

## Style
- Short sentences, one claim each.
- No "first", "novel" or "state-of-the-art".
- Criticise the input a category relies on, not a named tool.
- Copy Fan-lab sentences that apply directly (the user, 2026-11-03), such as STcompare's first sentence; see the
  sentence plan.

## For the Discussion
- How we differ from Asami.
- A shared axis is not new by itself (CoPro and Asami already compare samples on one).
- t is not linear in physical length.
- Axis error is an alternative explanation for shifts in t.
- The design has two mice per condition and one human donor.

## References for the Introduction (22)

Every reference below was read in full (Moses 2022 and CoPro as preprints). The numbers in brackets are the current
draft's reference numbers.

| Paragraph | Reference | What it supports |
|---|---|---|
| P1 | Ståhl 2016 [1] | Spatial transcriptomics measures expression in place |
| P1 | Moses & Pachter 2022 [3] | The field and its technologies |
| P1 | Zhai 2006 [17] | Nephrons wind through the tissue; tracing them takes serial sections |
| P1 | Schott 2024, Open-ST [21] | 3D spatial transcriptomics, over limited depth |
| P1, P3 | Chakravorty 2026 [38] | Seminiferous tubules; ordering around a known cycle |
| P2 | Curthoys & Moe 2014 [9] | S1–S3 transport; SGLT2 in S1/S2 and SGLT1 in S3 |
| P2 | Chen 2021 [10] | mRNA map across the 14 mouse tubule segments |
| P2 | Ransick 2019 [8] | Mouse PT regional states and injury susceptibility |
| P2 | Lake 2023 [13] | Human healthy and injured states, including spatial data |
| P3 | Zhao 2021, BayesSpace [24]; Hu 2021, SpaGCN [25] | Spatial domains |
| P3 | Chitra 2025, GASTON [27] | Gradients built from tissue x/y |
| P3 | Miao 2026, CoPro [29] | Continuous axis; was given segment order in its kidney run |
| P3 | Clifton 2023, STalign [30]; Clifton 2026, STcompare [31] | Comparison needs structure shared in the plane |
| P3, P4 | Trapnell 2014 [32]; Haghverdi 2016 [33] | Pseudotime; needs a start cell |
| P3 | Saelens 2019 [36] | Trajectory-inference benchmark |
| P3 | Niec 2022, SpaceFold [37] | Ordering by cell-type composition |
| P3 | Asami 2026 [40] | Whole structures captured in 3D, ordered by developmental age |
| P5 | Albergante 2020, ElPiGraph [47] | Elastic principal curve |
| P5 | Oliveira 2025 [4] | Visium HD |

### Moved out of the Introduction
- **To Methods:** Leiden and Seurat; the foundation-model references; Harmony; Hastie & Stuetzle; Gorban &
  Zinovyev; McCarthy; tradeSeq.
- **To Results:** Zhai 2003; Slingshot, PAGA and scFates (as baselines); Kirita 2020 and Dixon 2022 (for AKI);
  Muto 2021 (for human).
- **To the Discussion:**
  - scPrisma;
  - Belayer;
  - Halpern 2017 and Moor 2018: liver and villus have landmarks in the plane, which is the contrast case.
- **Dropped:** the Method of the Year editorial, Lee 2015, Melo Ferreira 2021, CODA, STARmap, Vickovic 2022, GraphST.

### Links (DOIs checked against Crossref, 2026-11-03)
- Ståhl 2016: https://doi.org/10.1126/science.aaf2403
- Moses & Pachter 2022: https://doi.org/10.1038/s41592-022-01409-2
- Zhai 2006: https://doi.org/10.1681/ASN.2005080796
- Schott 2024 (Open-ST): https://doi.org/10.1016/j.cell.2024.05.055
- Chakravorty 2026: https://doi.org/10.1016/j.cell.2026.04.036
- Curthoys & Moe 2014: https://doi.org/10.2215/CJN.10391012
- Chen 2021: https://doi.org/10.1681/ASN.2020101406
- Ransick 2019: https://doi.org/10.1016/j.devcel.2019.10.005
- Lake 2023: https://doi.org/10.1038/s41586-023-05769-3
- BayesSpace: https://doi.org/10.1038/s41587-021-00935-2
- SpaGCN: https://doi.org/10.1038/s41592-021-01255-8
- GASTON: https://doi.org/10.1038/s41592-024-02503-3
- CoPro (preprint): https://doi.org/10.64898/2026.04.17.719309
- STalign: https://doi.org/10.1038/s41467-023-43915-7
- STcompare: https://doi.org/10.1093/bioinformatics/btag644
- Monocle (Trapnell 2014): https://doi.org/10.1038/nbt.2859
- DPT (Haghverdi 2016): https://doi.org/10.1038/nmeth.3971
- Saelens 2019: https://doi.org/10.1038/s41587-019-0071-9
- SpaceFold (Niec 2022): https://doi.org/10.1016/j.stem.2022.05.007
- Asami 2026: https://doi.org/10.1016/j.cell.2026.06.014
- ElPiGraph (Albergante 2020): https://doi.org/10.3390/e22030296
- Visium HD (Oliveira 2025): https://doi.org/10.1038/s41588-025-02193-3
