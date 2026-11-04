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
