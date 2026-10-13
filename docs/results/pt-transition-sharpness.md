# Sharp or gradual: how segment-level pathways change along the elastic curve (notebook 85)

Frozen on 2026-10-12, before any transition width, binned profile or label width of notebook 85 was computed. The user
asked, after notebook 84: "can we further classify into sharp transitions or gradual transitions?"

**What I had seen before freezing.**
- Notebook 84's results (thirds inside segments) and its called pathways.
- Notebook 83's figures (71's spline curves).
- Notebook 81's descriptors: "no sharp switch (width ≤ 0.5) is confirmed in either species".
- 72's `transitions.csv` (label transition positions per specimen) and `noise_budget.csv`.
- The number of PT units per specimen: Ctrl1A2 2,644 and Ctrl1A4 2,650.

## Two findings from reading the code (no fit), which set the design

1. **The spline may not be able to show a step.**
   - 71's gene curves, which 81, 83 and 84 read, are NB GLMs on a 6-df cubic B-spline. Its three interior knots sit at
     the position quartiles.
   - A spline this stiff may draw a true step as a graded change. If so, 81's "sharp" class (≤ 0.5 segment units) was
     out of reach by construction.
   - Q1 below measures this. Notebook 85 measures widths without the spline.
2. **Pooled segment units misplace the mouse.**
   - Notebooks 79–84 read positions in 72's pooled units: the mean of all four specimens' label transitions (0.473,
     0.605).
   - The mouse's own S1→S2 label transition is at curve position 0.343 and 0.340. In pooled units that is S1 0.72.
   - So a mouse "switch inside S1 at 0.69–0.72" (notebook 83: NaPi-IIc, U-M3) sits at the mouse's own S1/S2 boundary.
   - Pooled S2 [1, 2) covers only the second half of the mouse's own S2. That was the "S2" of 84's thirds.
   - Notebook 85 reads mouse positions in **mouse units**: the mean of the two mice's own logistic transitions. S0
     converts the positions that 81, 83 and 84 quoted.

## Question

Chen 2021's segment-level analysis calls pathways "enriched in PTS1, PTS2 or PTS3". At that resolution a step at a
segment boundary and a gradient across segments look alike. Along the elastic curve:
- Is each pathway's change **sharp** (done within half a segment) or **gradual** (spread over a segment or more)?
- **Where** is it: at a mouse label boundary or inside a segment?
- Does the curve resolve that difference at that place?

## Labels

- Notebook: `analysis/notebooks/85_pt_transition_sharpness.ipynb`
- Logic: `85.transition_sharpness.1`
- Module: `pseudospace/transition_sharpness.py`, tested in `tests/test_transition_sharpness.py`
- Results: `results/pt_transition_sharpness/` (`tables/`, `figures/`, `run_record.json`)

## Inputs (read only)

| Source | What 85 reads |
|---|---|
| 70, 71 | `package_inputs.hvm_inputs` (mt-free, as 71, 81, 84): counts, library, curve position, labels, universe. 71's stage cache: the condition model's convergence and separation flags and its spline knots |
| 72 | `transitions.csv` (the per-specimen label transitions) |
| 84 | `tables/segment_pathways.csv` (the Chen-called pathways, their cores and representatives) and `run_record.json` (logic 84.within_segment_pathways.1) |
| Chen 2021 | the GSE150338 TPM file (sha256 `daf8b6f4…5488`), only for the segment means drawn beside the curves |

**Guards (not cached).**
- As 84: G1 is the Chen sha256; G2 checks that `hvm_inputs` reproduces 71's universe and units; G3 checks that the
  logistic transitions reproduce 72's.
- G4: 84's pathways and cores are re-derived from its `segment_pathways.csv`. Every core gene must be usable in both
  mice by 81's rule.

## Coordinate and units

- **Coordinate:** the package's elastic curve (stage 04 `pseudospace`), the only coordinate (user instruction).
- **Mouse units:** the mean of the two mice's own logistic label transitions (72's method), mapped piecewise
  linearly: [0, t12) → [0, 1), [t12, t23) → [1, 2), [t23, 1] → [2, 3].
- **Human:** not analysed. The human label transitions cross on the curve (HUK1_MED1: S1→S2 0.667 after S2→S3 0.590;
  issue 6), so human segment units are undefined.

## P · Profiles (no spline)

- **P1 · Bins.** Per mouse, PT units are sorted by curve position and cut into 40 equal-count bins of about 66 units.
  A bin's position is the median position of its units, in mouse units.
- **P2 · Gene profile.** Per bin: log2((Σ counts + 0.5) / Σ library). The library is `hvm_inputs`' library (mt-free),
  the same offset as 71.
- **P3 · Centring.** Each gene's profile minus its mean over bins, then minus the mouse's reference profile. The
  reference is the median over usable genes of their centred profiles (81's deviation 1; `pf.reference_curve`).
- **P4 · Pathway profile.** The median over a pathway's core genes of their P3 profiles.

Usable means 81's rule: 71's curve converged and not separated, and detection ≥ 5% in that mouse.

## E · Edges and their width

- **E1 · Models.** Both are fitted by least squares on the 40 bins, with x in mouse units:
  - M1, one edge: y = a + b·σ((x − m)/s);
  - M2, two edges: y = a + b1·σ((x − m1)/s1) + b2·σ((x − m2)/s2), with m1 < m2.

  σ is the logistic function. The width of an edge is its 10–90% distance, w = ln(81)·s, bounded to [0.05, 6] mouse
  units. A midpoint m must lie between the second and the second-to-last bin positions. The fit is an exhaustive grid
  over (m, w), linear in the amplitudes, then refined with bounded least squares from the best grid point.
- **E2 · Model choice.** M2 when BIC(M2) ≤ BIC(M1) − 6 and both M2 edges reach |b| ≥ 0.25 log2; otherwise M1.
- **E3 · Counted edges.** An edge counts when |b| ≥ 0.25 log2 (84's flat threshold). A profile without a counted edge
  is **flat**.
- **E4 · Class** of an edge, from its width (81's thresholds):

  | Class | Width (mouse units) |
  |---|---|
  | sharp | w ≤ 0.5 |
  | graded | 0.5 < w < 1 |
  | gradual | w ≥ 1 |

- **E5 · Location** from the midpoint: "at S1/S2" when |m − 1| ≤ 0.15, "at S2/S3" when |m − 2| ≤ 0.15, otherwise
  "inside" the segment that holds m.

## R · Replication (the two mice)

- **R1 · Agreement.** An edge in one mouse agrees with the other mouse's edge of the same direction whose midpoint is
  nearest, if they lie within 0.25 mouse units. Matching is one-to-one, smallest distance first.
- **R2 · Agreed class.** The class of the two mice's mean width and the location of their mean midpoint. When one
  mouse is sharp and the other gradual, the edge is "mice disagree on width" and gets no class.
- **R3 · Main edge** of a pathway: its agreed edge with the largest mean |b|. A pathway with no agreed edge is "no
  agreed edge" (flat in both mice is "flat").

## Q · What the curve can resolve

A gradual edge may be a step that the coordinate blurred. Blur only widens edges, so a sharp agreed edge needs no
reference. A gradual or graded edge does.
- **Q1 · Spline step response.**
  - Input: a noiseless unit-level step at x0 ∈ {0.5, 1, 1.5, 2, 2.5} mouse units, on each mouse's real unit positions.
  - Fit: least squares on 71's spline basis with 71's stored knots, evaluated on 71's grid.
  - Output: the 10–90% width. It shows whether 81, 83 and 84 could have shown a sharp switch. This is an
    approximation: 71 fits on the log link of NB counts.
- **Q2 · Binning step response.** The same synthetic step through P1 and E1 (M1). This is the floor of 85's own method.
- **Q3 · Label step.**
  - Per mouse and boundary, a bin's label fraction: S2 or S3 units for S1/S2, S3 units for S2/S3.
  - The fraction is fitted with M1, and its width says how sharply a discrete segment label switches on this curve.
  - This is the empirical step response at the two boundaries.
- **Q4 · Interior resolution.**
  - Agreed gene edges (G1) whose mean midpoint lies inside a segment (outside the ±0.15 boundary windows).
  - The segment's interior reference is the 5th percentile of their mean widths. It needs ≥ 20 such edges; otherwise
    the interior is "unresolved".
- **Rule.** A graded or gradual edge is **resolved** when the reference at its location (Q3 at a boundary, Q4 inside)
  is ≤ 0.5. Otherwise it reads "graded or gradual, or a blurred step (unresolved here)".

## G · Gene level (descriptive)

- **G1 · Edges for every gene** usable in both mice: P2–P3, E1–E5 and R1–R2, as for pathways.
- **G2 · Composition of a pathway edge.**
  - Take the core genes with an agreed edge of the same direction whose mean midpoint lies within 0.5 of the pathway
    edge's.
  - Record the share of them that is sharp, graded and gradual, and the interquartile range of their midpoints.
  - Label for a gradual or graded pathway edge:
    - "staggered gene switches" when ≥ 50% of those genes are sharp and the interquartile range is ≥ 0.25;
    - "gene gradient" when ≥ 50% are gradual;
    - otherwise "mixed".
  - With fewer than 3 such genes, it is "not evaluable".

## S0 · Earlier positions re-read in mouse units

A table converts the mouse positions quoted by 81, 83 and 84 from pooled to mouse units:
- NaPi-IIa at S3 0.37 and NaPi-IIc at S1 0.69 (83);
- S-M1 at S3 0.13 (83) and S3 0.18 (81);
- U-M3 at S1 0.72 (83) and S1 0.80 (81);
- 84's pooled S1 and S2.

This is arithmetic on the cuts, not a fit.

## Pathways

84's Chen design, unchanged:
- every called pathway with ≥ 5 core genes: S1 9, S2 53, S3 1;
- the representatives for the figures: 5, 30, 1;
- each pathway's core is its 84 core.

Every one of them is classified. There is no between-pathway test.

## Bounds

- **C14.** Mouse fold noise ÷ within-segment spread is 4.2 in S1, 2.0 in S2 and 0.98 in S3 (pooled units). Inside a
  segment, Q4 decides whether a width is read.
- **Cohort.** Two male control mice. Edges are agreed between them, not tested.
- **Descriptive.** Every class is a description. There are no p-values.

## Figures

- **85a** (discrete beside continuous):
  - Left: Chen's PTS1, PTS2 and PTS3 means of each representative's core.
  - Right: the representative's binned profile, the mean of the two mice in mouse units.
  - A mark on each agreed edge: filled for sharp, half-filled for graded, open for gradual. Unresolved edges get a
    grey mark.
  - Rows are grouped by Chen segment, then by main-edge class and position.
- **85b** (examples): within one Chen segment class, a sharp pathway and a gradual one. Each shows both mice's binned
  profiles with the fitted edges, and the label fraction at the nearest boundary in grey.
- **85c** (resolution):
  - Q1, the spline and binning step responses;
  - Q3, the label fractions with their fits;
  - Q4, agreed gene edge widths by location.

## What would make us say the classification adds nothing

- No pathway has an agreed edge: there is nothing to classify.
- Every agreed main edge has the same class: the curve does not split Chen's segment classes by sharpness.
- Q3 is above 0.5 at both boundaries and Q4 is unresolved in every segment: on this curve, sharp and gradual cannot be
  told apart, except where an edge is already sharp.
