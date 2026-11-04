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

## Cut from the current draft (787 words)
| Material | Where it goes |
|---|---|
| Visium HD specifications | Methods |
| "distinguished by transporter expression, not ultrastructure" | Results |
| Descriptions of CoPro, SpaceFold, Chakravorty and Asami | Discussion (keep one citation clause each here) |
| The segmentation pipeline | Methods (keep half a clause here) |
| The URL | Code availability |

## Style
- Short sentences, one claim each.
- No "first", "novel" or "state-of-the-art".
- Criticise the input a category relies on, not a named tool.
- The first sentence must not echo STalign's ("Spatial transcriptomics technologies enable…").

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
