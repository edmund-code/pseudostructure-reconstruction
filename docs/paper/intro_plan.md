# Introduction plan (Nature Methods Article)

Compiled 2026-11-03, revised the same day for a shorter target. This is a plan, not prose.
- Outline: GASTON (*Nat. Methods* 2025).
- Substance and tone: the Fan lab's Introductions, whose median is 477 words.
- Reading notes: `docs/paper/citation_checks/B_intro_refs_22_40.md` and `F_fanlab_intros.md`.

**Target: 500 words, 5 paragraphs, about 35 references.** Each paragraph has one job. Anything a reader doesn't need
in order to accept the gap goes to Methods, Results or the Discussion.

## P1. Biology and the gap (~90 words)
- Organs are built from repeated units whose function changes along one axis: lobule, villus, nephron.
- The proximal tubule's S1–S3 differ in transport, metabolism and injury. Function changes at segment boundaries and
  within segments; for example, glucose uptake passes from SGLT2 to SGLT1.
- The gap: atlases resolve the PT only as three segments.

## P2. The data limit (~80 words)
- Spatial transcriptomics measures expression in place. But a section cuts each nephron many times, and neighbouring
  cross-sections belong to different nephrons, so tissue position does not give position along the structure.
- Tracing nephrons through serial sections, or profiling in 3D, recovers position but does not scale to whole organs.

## P3. Prior methods: two categories (~140 words)
- **From the tissue plane.**
  - Clustering and spatial domains average within segments, so they miss within-segment change.
  - Gradient and alignment methods need tissue coordinates or shared structure in the plane.
- **From expression.**
  - Trajectory inference needs a start cell.
  - Repeated-structure ordering needs a change in cell-type composition, a known cycle, or whole structures captured
    in 3D.
- **Close:** the PT has one cell type, no landmarks and partial sampling at unknown positions, so none of these
  applies.

## P4. The concept (~80 words)
- Pseudostructure treats each cross-section as a sample of one shared structure and places it by its expression at a
  position t that carries both segment order and position within a segment (Fig. 1).
- As pseudotime orders cells along a process, t orders cross-sections along a structure. It preserves order, not
  length.
- One sentence on the assumption: the copies are exchangeable and sampled along one unbranched axis.

## P5. The tool and results (~110 words)
- "Here, we developed pseudostructure reconstruction, a label-free method that …". Name the elastic principal curve,
  the NB gene models and the relabeling calibration. Segmentation detail goes to Methods.
- The data, in one clause: 4 mouse sections (2 control, 2 AKI) and 2 human sections.
- One clause per result:
  - recovers segment order without labels;
  - resolves within-segment genes;
  - finds human–mouse reversals and later human transporter switch points (**genes**, not pathways);
  - finds AKI switch points moving upstream.

## Cut from the current draft (787 words)
- Visium HD grid specifications. These go to Methods.
- "Distinguished by transporter expression, not ultrastructure". This goes to Results.
- Named descriptions of CoPro, SpaceFold, Chakravorty and Asami. Keep one citation clause each, and put the
  distinctions in the Discussion.
- The segmentation pipeline sentence. Keep half a clause.
- The URL. It goes in Code availability.

## Style rules
- Short sentences, one claim each. No "first", "novel" or "state-of-the-art".
- Criticise the input a category relies on, not a named tool.
- The first sentence must not echo STalign's ("Spatial transcriptomics technologies enable…").

## Left for the Discussion
- How we differ from Asami.
- A shared axis is not new by itself (CoPro, Asami).
- t is not linear in physical length.
- Axis error as an alternative explanation for shifts in t.
- The two-mouse, one-donor design.
