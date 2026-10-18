# Transporter handoffs along the PT: a family-wide screen (notebook 101)

**The question.** The outline's Results section 3 claims "micro-gradients in solute handling". Notebook 92 mapped one
handoff along the continuous axis: glucose uptake passing from SGLT2 to SGLT1. Is that a special case?

Within each solute-carrier (SLC) and ABC family, do early (S1-high) and late (S3-high) members hand off more often
than chance? Where along the axis do the handoffs sit, and are they conserved between mouse and human? This is a
screen over every transporter family. Nobody can see it with S1/S2/S3 labels, because a handoff's position and width
are continuous quantities.

**What the literature fixes (notebook 92's reads, full text where reachable):**
- Rodent SGLT2 sits in S1/S2 and is absent from S3 (Sabolic 2012, PMID 22262063; Vallon 2011, PMID 20616166).
- Human SGLT1 is strongest in S3, while mouse SGLT1 is reported higher in S2 than S3 (Vallon 2020, PMID 32144488).
- GLUT2 is strongest in S1/S2, and GLUT1 runs along the PT, highest distally (Vallon 2020).
- Human SGLT2 and SGLT1 have similar glucose affinity (Hummel 2011, PMID 20980548). A handoff is therefore a matter of
  position and abundance, and no affinity design principle is assumed here.

Protocol frozen on 2026-10-18, before any transporter pair was formed or profiled.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/101_pt_transporter_handoffs.ipynb`; logic `101.handoffs.1`; seed 101.
- Results `results/pt_transporter_handoffs/`.
- New helpers go in `pseudospace/positional_decomposition.py` or a small module `pseudospace/handoffs.py`, with tests.
- **Gradients:** notebook 90's per-specimen end contrasts E (`d_end_contrasts.csv.gz`; Ctrl1A2, Ctrl1A4, HUK1_COR1,
  HUK1_MED1; with its fold-protected flag). E is the mean of the last 30% minus the first 30% of each specimen's own
  curve, in log2.
- **Curves:** notebook 71's human vs mouse species gene model, read from its cache. These are the specimen curves on
  the common grid 0.05–0.95 (natural-log rate).
- **Families:** HGNC group branches 752 (solute carriers) and 417 (ABC transporters) under
  `data/external/hgnc/`. The family is the `Group name` (for example "Solute carrier family 5"). Human symbols map to
  mouse through the one-to-one ortholog map.
- **External:**
  - Chen 2021 microdissection (`data/external/segment_microdissection/GSE150338_Mouse_TPM_Replicates.txt.gz`; male
    mice; PTS1, PTS2, PTS3 replicates);
  - Census human PT donors (`data/external/census_pt_segments/`; convoluted PT against segment 3, donors with ≥ 50
    cells in each).

## 2 · Definitions (fixed now)

- **Eligible genes:** family members with E in notebook 90's table and finite curves in both specimens of the species.
- **Early in a species:** E ≤ −0.5 in both of its specimens. **Late:** E ≥ +0.5 in both.
- **Handoff pair in a species:** one early and one late member of the same family. Each pair is reported in both
  species.
- **Profiles:** each specimen curve converted to log2 and min–max scaled over the grid points within that specimen's
  range. The species profile is the mean of its two specimens.
- **Crossing x\*:** the first grid position where the early profile minus the late profile changes sign from positive
  to negative.
- **Handoff width:** x(late first reaches 0.5) − x(early last falls below 0.5). A negative width means the halves
  overlap (notebook 92's question for glucose).
- **Segment cuts:** notebook 71's zone cuts for human vs mouse (`method_b_zones.csv`). A crossing is *inside a
  segment* when it lies more than 0.05 from both cuts.

## 3 · Questions and rules

**Q1 · Are handoffs enriched within families?**
- Statistic: the number of handoff pairs within families, per species.
- Null: family labels permuted among the eligible transporter genes of that species, family sizes kept, 1,000 times.
  The empirical p is one-sided.
- **Holds** in a species when p ≤ 0.05.
- Also reported: the number of families with at least one handoff, and the share of early and late genes among
  transporters against all eligible genes.

**Q2 · Does microdissection or the atlas confirm the direction?**
- Mouse: for each mouse handoff pair, Chen 2021's log2 (mean PTS1 TPM + 1) / (mean PTS3 TPM + 1). The pair is
  **confirmed** when the early member has a ratio > 0 and the late member < 0.
- Human: the same with Census donor-mean log2 CPM, convoluted against segment 3.
- **Holds** when ≥ 2/3 of a species' pairs are confirmed.

**Q3 · Conservation (descriptive).**
- Pairs that are handoffs in both species with the same orientation, and pairs whose orientation reverses.
- Against the expected count, the product of the two species' pair rates times the number of within-family pairs.
- Crossing positions, mouse against human, for conserved pairs.

**Q4 · Where the handoffs sit (continuous against discrete; descriptive).**
- The share of crossings inside a segment, and the width distribution.
- Against the same for random early × late gene pairs from different families. This asks whether transporter handoffs
  sit in particular places.

**Recovery checks (named now).**
- R1: Slc2a2 (GLUT2) early and Slc2a1 (GLUT1) late, in both species.
- R2: Slc5a2 (SGLT2) early and Slc5a1 (SGLT1) late, in human.
- R3 (mouse; descriptive, since the literature conflicts): Slc5a2 against Slc5a1.

A miss is reported, not tuned. Slc5a2 helped build the coordinate (notebook 92), so R2 and R3 are partly circular.

## 4 · Outputs and figures

- **Tables:**
  - `eligible.csv` (gene, family, E per specimen, early or late per species);
  - `pairs.csv` (family, early, late, species, x\*, width, Chen and Census ratios, confirmed);
  - `q1_enrichment.csv`, `q2_confirmation.csv`, `q3_conservation.csv`, `q4_positions.csv`;
  - `recovery_checks.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 101a: a handoff ladder per species, each pair a segment from the early member's 50% fall to the late member's 50%
    rise, with the cuts;
  - 101b: profiles of the recovered and the strongest pairs, mouse and human;
  - 101c: crossing x\*, mouse against human, for pairs in both species;
  - 101d: the Q1 null histograms.

## 5 · Caveats

- mRNA, not protein. Sabolic 2012 shows Sglt2 mRNA and SGLT2 protein can move in opposite directions between sexes.
- One human donor; human early and late calls are replicated across two sections only.
- Within-family pairs share sequence and sometimes function, but not necessarily substrate. A pair is a positional
  handoff, not a proven functional one.
- Genes that built the coordinate make their own positions partly circular. Notebook 90's fold-protected flag is
  reported per gene.

## Notebook 101 results

`analysis/notebooks/101_pt_transporter_handoffs.ipynb` (logic `101.handoffs.1`) writes to
`results/pt_transporter_handoffs/`. It followed the protocol frozen in 66afdb6 and ran in 5 s. G1 passed: 71's curves
for 11,055 genes on 0.050–0.950, with cuts 0.304 and 0.620. Module `pseudospace/handoffs.py` with
`tests/test_handoffs.py` (2 tests).

**How the protocol was read.**
- A transporter in two HGNC groups keeps its first group.
- Chen's ratio uses the mean TPM over the replicates.
- Census uses the 6 donors with ≥ 50 cells of both types.
- The random comparison pairs are non-transporter genes.

**Eligible transporters.**
- 479 transporters have a one-to-one ortholog, in 77 families. 237 are eligible in each species.
- Mouse: 50 early and 46 late. Human: 31 early and 24 late.
- Transporters are graded far more often than genes overall: 0.41 against 0.17 in mouse, 0.23 against 0.08 in human.

| Question | Mouse | Human | Holds |
|---|---|---|---|
| Q1 within-family handoff pairs (null mean, p) | 52 (78.1, 0.95) | 20 (24.9, 0.72) | **no** |
| Q2 direction confirmed (Chen; Census) | 50 / 52 = 0.96 | 16 / 20 = 0.80 | yes |
| Q3 conserved pairs (expected by chance) | 11 (0.55); reversed 0 | | descriptive |
| Q4 crossings inside a segment (random pairs) | 0.60 (0.67) | 0.50 (0.72) | descriptive |
| Q4 overlapping handoffs, width < 0 (random) | 0.64 (median width −0.07; random 0.00) | 0.95 (−0.15; random −0.03) | descriptive |

- **Q1 fails.** Within a family, early and late members do not co-occur more than chance. They co-occur *less*, so
  families tend to share one direction. Handoffs are not a family-wide design principle of the PT.
- **Specific handoffs are real and conserved.**
  - Microdissection confirms 96% of the mouse pairs, and the atlas 80% of the human pairs.
  - 11 pairs are handoffs in both species, against 0.55 by chance, with none reversed:
    - amino acids: B0AT1 → B0AT3 (Slc6a19 → Slc6a18), and b0,+AT, LAT2 and y+LAT1 → AGT1 (Slc7a9, Slc7a8 and Slc7a7
      → Slc7a13);
    - organic anions: OAT3 → OAT2, URAT1 and OAT5/Slc22a13 (Slc22a8 → Slc22a7, Slc22a12, Slc22a13), and Slc22a23 → the
      same three;
    - sugars: SGLT2 → SGLT5 (Slc5a2 → Slc5a10).
- **The handoffs overlap and sit inside S2.**
  - Mouse crossings run from x 0.30 to 0.66, almost all between the two cuts. Median x\* is 0.50 in mouse and 0.56 in
    human.
  - Most handoffs overlap: the late member rises before the early one falls (64% in mouse, 95% in human).
  - A segment view puts both members' changes into "S2" and cannot show where, or how wide, the handoff is.
- **Recovery checks missed.**
  - R1 GLUT2 → GLUT1: Slc2a2 is early in both species, but Slc2a1 is not late in either. Its gradient does not
    replicate in our mice (E −0.44, 0.83) and is flat in human.
  - R2 and R3 SGLT2 → SGLT1: Slc5a1 is not in notebook 90's eligible table, so it was not testable here. Notebook 92
    studied it separately.
  - Instead, SGLT2 hands off to SGLT5 (Slc5a10) in both species.

**What this answers.**
1. **Transporter gradients are common, but handoffs within a family are not a general rule.** Transporters are 2–3
   times more often graded than other genes, and family members tend to share a direction. The design principle "an
   early paralog hands off to a late one" is not supported genome-wide.
2. **A small set of conserved, microdissection-confirmed handoffs exists,** in amino-acid, organic-anion and sugar
   transport. They overlap rather than switch, and they sit inside S2. Only a continuous axis can show where and how
   wide each handoff is, so these are concrete "micro-gradients in solute handling" for the paper. The claim should be
   these pairs, not a general principle.
3. **As a method:** a family-wide screen with a label-permutation null separated a general claim (rejected) from
   specific findings (11 conserved pairs against 0.55 by chance), with external confirmation.
