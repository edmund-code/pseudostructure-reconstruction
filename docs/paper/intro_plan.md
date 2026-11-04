# Introduction plan (Nature Methods Article)

Compiled 2026-11-03. This is a plan in bullet points, not draft prose. It gives the outline, taken from GASTON
(Chitra et al., *Nat. Methods* 2025), and the substance and tone, taken from the Fan lab's Introductions.

Sources, all read in full text:
- GASTON and the other structure examples: `docs/paper/citation_checks/B_intro_refs_22_40.md`.
- Nine Fan-lab Introductions, plus STalign and CRAWDAD: `docs/paper/citation_checks/F_fanlab_intros.md`.

## Size

- Nature Methods sets no length for the Introduction. It has no heading and counts toward the 3,000-word main text.
- GASTON runs 743 words in 5 paragraphs and PASTE 756 words in 5. Fan-lab Introductions run 320–650 words in 2–4
  paragraphs; their Nature Methods paper (PAGODA) has 477.
- Target: **5 paragraphs, 700–800 words**, about 45 references.
- The current draft has 787 words in 5 paragraphs (137 / 98 / 259 / 125 / 168). The length is right; the problem is
  the order:
  - P3 carries every prior method at once.
  - P4 fuses the concept with the pipeline.

## The outline (GASTON's five moves)

| GASTON | Ours | Words |
|---|---|---|
| P1 Biology: expression varies both as sharp boundaries and as continuous gradients | P1 Repeated structures and what is known along the PT | ~130 |
| P2 The technology and its data limit (sparsity) | P2 Spatial transcriptomics and its data limit (position along the structure is not observed) | ~120 |
| P3 Prior methods in two categories; neither models X | P3 Methods that use the tissue plane, and methods that order from expression | ~220 |
| P4 New concept (topography, isodepth) explained by an analogy | P4 Pseudostructure: the concept, the coordinate t, the assumptions | ~120 |
| P5 Tool plus a results summary | P5 Pseudostructure reconstruction: the pipeline, the data, the results in figure order | ~170 |

## P1. Repeated structures and what is known along the PT (~130 words)

Cover:
- Many organs are built from repeated units whose function changes continuously along one axis: the liver lobule,
  the intestinal villus, the kidney nephron [5–8].
- Along the nephron, the proximal tubule runs through S1, S2 and S3. These differ in transport and metabolism
  [9,10] and in injury susceptibility [8]. In the mouse they are told apart by transporter expression, not
  ultrastructure [11].
- The GASTON move: function changes both at segment boundaries and within segments. One clause, with a
  textbook example: glucose uptake passes from SGLT2 early to SGLT1 late (cite the transporter literature, not our
  data).
- Single-cell and spatial atlases map segment states in mouse, human and injury [8,12–16].
- End on the gap, in the Fan lab's phrasing ("remains an open question", "is not well understood"): expression along
  the tubule is known mostly at the resolution of three segments.

Avoid:
- Opening on "Spatial transcriptomics technologies enable…". It echoes STalign's first sentence (author note F).
  Open on the biology, as GASTON does.

## P2. The data limit: position along the structure is not observed (~120 words)

The Fan-lab funnel: "*X technologies enable Y. However, …*". GASTON gives the data limit its own paragraph, before
any method is criticised. Ours is that the variable we need is missing.

Cover:
- Spatial transcriptomics measures expression in place [1–3]. Visium HD covers the transcriptome on a 2-µm grid over
  a 6.5 × 6.5-mm capture area [4].
- The limit: a two-dimensional section cuts each nephron many times.
  - A section holds thousands of tubule cross-sections, each from an unknown position.
  - Neighbouring cross-sections often belong to different nephrons.
  - So tissue x/y does not give position along the structure.
- Recovering position physically needs tracing through hundreds to thousands of serial sections [11,17,18] or 3D
  profiling [19–21]. These have captured morphology, targeted panels or sub-millimetre depths, not whole organs.
  Use a "remains challenging" register, not "laborious".

## P3. Two categories of prior methods (~220 words)

This is the paragraph to fix. GASTON uses two categories with citation blocks. The Fan lab criticises the **input a
category relies on**, closes each category with "As such" or "Thus", and adds a concession. Name only the closest
competitors.

**Category 1: methods that work in the tissue plane.** The shared input is tissue coordinates or shared structure in
the plane.
- Clustering assigns cross-sections to discrete groups that are annotated as segments [22,23]. Spatial-domain methods
  smooth across neighbours [24–26]. Conditions are then compared on group means.
  - What they lack: a group mean "assumes mean expression is constant within" the group (GASTON's own phrase, citable
    to [27]). It cannot represent a within-segment change, a boundary peak or a handoff.
  - Their prior misleads here: a neighbour in the plane is often another nephron.
- Continuous-gradient methods (GASTON, Belayer, CoPro) build the coordinate as a **smooth function of tissue x/y**
  [27–29]. The objection is smoothness, not contiguity: GASTON allows non-contiguous domains. In its kidney run,
  CoPro was given the segment order.
- Alignment and comparison methods (STalign, STcompare) need sections that share structure in the plane [30,31].
  The concession clause: they resolve the kidney at compartment scale (cortex and medulla), which is what they were
  built for. Be gentle; these are the lab's own tools.

**Category 2: methods that order from expression.** The shared input is a known process, cycle or whole structure.
- Trajectory inference orders cells along processes from expression. It is usually given a start cell or cluster
  [32–36]. Concession: it is the closest analogue, and we borrow its logic.
- Repeated-structure ordering:
  - crypt positions from deconvolved cell-type composition (SpaceFold) [37];
  - testis tubules around the known seminiferous cycle (Chakravorty, with the scPrisma template) [38,39];
  - hair follicles by developmental age, each captured whole in 3D (Asami) [40].
  - "As such", they rely on a change in composition, a known cyclic topology or whole-structure capture.

**Close with the conditions the PT meets and none of these handles:**
- one epithelial cell type along the axis;
- an open axis with no landmark marking where a cross-section lies;
- every section sampling partial, unknown positions.
- Fan-lab ending: "*As such, an approach that … provides an alternative strategy when …*" (STdeconvolve).

Avoid:
- "first", "novel" and "state-of-the-art". None of the nine Fan-lab Introductions in report F uses them.
- Claiming that ordering cross-sections is new (author note E). The new part is an axis built with no tissue
  coordinates, labels, reference program or whole-structure capture.

## P4. The concept: pseudostructure (~120 words)

GASTON names its concept ("gene expression topography", the isodepth) one paragraph before the tool. Fan-lab Nature
Methods and Nature Communications papers explain the model by an analogy (PAGODA: GSEA; STdeconvolve: LDA's
words → documents → topics). The current draft fuses this paragraph with P5.

Cover:
- **Name and define.** Pseudostructure treats every cross-section of a repeated structure as a sample of one
  shared structure and places it by its expression.
- **The coordinate, defined the way GASTON defines the isodepth:** one position t ∈ [0, 1] that carries both the
  order of segments and the position within a segment.
- **The analogy.** As pseudotime orders cells along a process, pseudostructure orders cross-sections along a structure.
  - DPT's accepted text: "pseudotimes are proxies for universal time".
  - Likewise, t is a proxy for position: it preserves order, not physical length.
  - Asami already used a cryo-EM analogy for follicles. Don't reuse it without citing.
- **Assumptions as properties of the data**, in the STdeconvolve style ("these unique aspects … include (i)…(iv)"):
  1. The copies are exchangeable samples of one shared expression program. This is a modelling assumption, as
     the workflow doc says.
  2. Each section holds many cross-sections, at positions spread along the axis, with partial coverage allowed.
  3. Within the scope analysed, the program changes along one unbranched axis.
  - The things it does not identify (orientation, physical length) go to the Discussion's limitations paragraph.
- **Generality**, claimed by listing structures, not by adjectives: lobules, villi, crypts, nephrons. Use the
  CRAWDAD form ("amenable to any … for which … can be obtained"): any tissue where repeated structures can be
  segmented and profiled.
- A **Fig. 1** callout on the concept. STdeconvolve puts Fig. 1 on its tool sentence, and CRAWDAD cites Fig. 1a,b
  for its concept.

## P5. The tool and results summary (~170 words)

The Fan-lab tool sentence: "*Here, we developed X, a [property] method for [task] using [mechanism] (Fig. 1)*". Both
GASTON and PASTE then give a results summary against comparators. Fan-lab Introductions that summarise results do
it in one or two "*we show / demonstrate that…*" sentences.

Cover:
- The tool sentence: a label-free method for placing cross-sections of a repeated structure along their shared axis,
  using an elastic principal curve.
- The pipeline in one sentence, as GASTON does:
  - H&E segmentation (frozen pathology foundation model plus a trained decoder) [41–43];
  - sum the Visium HD bins in each cross-section;
  - harmonise [44];
  - one elastic principal curve [45–47];
  - negative-binomial spline gene models, each gene on an axis built without it [48,49];
  - comparisons calibrated by balanced specimen relabelings.
- The data: Visium HD of 4 mouse kidney sections (2 control, 2 after ischaemia–reperfusion) and 2 human cortex
  sections.
- The results, in figure order, one clause each:
  - Without segment labels it recovers the segment order in every section, stays stable under downsampling and
    half-gene axes, and orders segments as well as the best tree-based method while fitting one unbranched curve.
  - It resolves within-segment and boundary genes that clustering obscures.
  - Species: genes whose PT gradient reverses, and transporter switch points that sit later along the human PT
    (glucose, amino acids). Say **genes**, not pathways (author note D2).
  - AKI: switch points move upstream as injured S3 loses late transporters; injury creates or erases within-segment
    gradients.
- The URL can follow the tool sentence (STdeconvolve style) or close the paragraph. GASTON and PASTE leave it to
  Code availability. One clause either way.

## Style rules taken from the reading

- Limitations are capabilities the methods lack, or inputs they rely on: "rely on", "do not model", "cannot
  represent", "are limited to". Hedge with "may" where the failure is conditional.
- Give one shared limitation per category, introduced with "As such" or "Thus", plus one concession for each
  category.
- The tool is "an alternative strategy when …", never "the first".
- Claim generality by listing data and structures.
- Assumptions go in the Introduction as data properties. Limitations go in a Discussion paragraph: "*Several
  limitations of … should be considered*" (VeloViz, SEraster).
- The Discussion's first paragraph restates the need, then "*Here, we presented …*", then the results in order
  (Fan-lab habit).

## Left for the Discussion, not the Introduction

- Asami 2026 orders whole structures that develop asynchronously, by developmental age. We order cross-sections of
  mature structures by position, from 2D sections in which no structure is whole (author note F).
- CoPro and Asami already compare samples on a shared axis. A shared axis is not new by itself.
- t is not linear in physical length. Asami's line is the model: "may not have a linear relationship with real time".
- A shift in t between groups could reflect axis error. Our answers are the held-out-gene axes and the relabeling
  null, and STcompare's alignment-error caveat is the model sentence.
- The two-mouse, one-donor design limits generalisation.
