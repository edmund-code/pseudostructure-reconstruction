# Question bank: what COAST can answer about nephron gradients (protocol 2026-11-10, then results)

The question list comes from an external brief: "High-value scientific questions enabled by fine-resolution nephron
gene-expression gradients in healthy human versus mouse and mouse AKI". It is a question list, not a
literature source; none of its claims are cited from it.

Everything here is **post hoc and descriptive**. No readout feeds the Fig. 5 or Fig. 6 calls. The readouts below were
fixed before any of them was computed. Code: `analysis/paper/question_bank.py`.

## Frozen protocol (2026-11-10, before any readout)

### Data and axes
- **PT questions** reuse the Fig. 5 and Fig. 6 screen curves (`fig{5,6}/screen/all_curves.npz`).
  - These are NB curves on the gene-held-out axis, with x the control (mouse) quantile along the PT on a 201-point
    grid.
  - b12 and b23 are the medians over the gene patterns in `axes.csv`.
- **Distal and whole-nephron questions** use the package pseudostructure t of each run (`04_pseudospace`), over all
  labelled units.
  - Per section, 60 equal-width bins of t on [0, 1]. A bin needs at least 15 units.
  - Each bin's rate is log2(10^4 · Σcounts / Σlibrary + 0.5), where the library is non-mitochondrial counts.
  - Smoothing is a 3-bin running mean.
  - Ceiling: the axis is not gene-held-out, so a gene that helped build the axis is placed by it. Both groups
    share that axis, and the relabelings carry the same bias.
- **Runs.**
  - Human vs mouse: the `hvm` run, sections Ctrl1A2 and Ctrl1A4 against HUK1_COR1 and HUK1_MED1.
  - AKI: the `aki` run, sections Ctrl1A2 and Ctrl1A4 against IR2A2 and IR2A4.

### Readouts (per section)
- **Scaled profile.** Min–max over the question's window.
- **Peak.** Position of the maximum of the profile.
- **Domain.** First and last positions where the scaled profile is at least 0.5; width is the difference.
- **Crossing.** The switch point of two scaled profiles. As in Fig. 5: the crossing nearest the window centre where
  the first profile goes from above to below the second.
- **Overlap.** The share of the window where both scaled profiles are at least 0.5.
- **AKI response.** IR curve minus the mean control curve, in log2, at each position.
  - The response peak is the position of the largest |response|.
  - The front is the most proximal position where the response reaches half of its peak, in the peak's sign.

### Contrast and false-call control
- For any per-section readout Q, the group contrast is the mean of the two test sections minus the mean of the two
  reference sections.
- The two balanced relabelings are those of Figs. 5–6.
- **Call** a difference when:
  - all four cross pairs share the contrast's sign; and
  - |contrast| exceeds both relabelings' |contrast|.
- Positions are reported in units of the window: the PT for PT questions, and the nephron (t in [0, 1]) for distal
  ones.
- **Ranking criterion, fixed now.** Called before suggestive; within that, ranked by how directly the result answers
  a question clustering cannot.

### Questions not answerable here, stated in advance
- **Human medullary questions** (medullary TAL, OMCD/IMCD osmoadaptation). Both human sections are cortex.
- **Human-specific genes.** OAT4 (SLC22A11) and MATE2-K (SLC47A2) have no mouse ortholog in the panel.
- **DCT1 vs DCT2 by label.** The hvm run has only DCT1 and the aki run only DCT2. Positions along t are read
  instead of labels.
- **Repair fate** (adaptive vs failed repair arising from pre-injury positions). It needs lineage tracing or a later
  time point; our IR mice sit at about 12 h–2 d. Only where the markers sit now is reported.
- **Inverted or duplicated identity after repair.** Same reason, no repaired time point. The order of segment
  markers in IR is reported.

## Changes after the first run (made before reading any values; 2026-11-10)

- **Label registration (distal readouts).**
  - Problem: human distal labels sit at different t (cTAL median 0.43 against 0.48 in mouse), which would make
    every distal gene look shifted.
  - Fix: each section's t is mapped onto the mouse controls' t by piecewise-linear interpolation between its own
    segment-label medians.
  - Effect: distal positions are now read relative to each section's morphological segments.
- **PT end trim.** The screen splines overshoot at the axis ends (mouse SGLT2 rose again at x = 1), so PT readouts use
  x in [0.05, 0.95], as `fig6_interpret.py` does.
- **Curve floor.** The screen curves are unfloored log2, and a gene absent in control ran toward −∞ (KRT20 gave
  +34,000 log2). They are now floored at +0.5 per 10k, as the distal profiles are.
- **Crossing.** This is Fig. 5's rule, the first crossing, and not "nearest the window centre" as first written.
- **A stricter call ("strong").** The frozen rule (four pairs agree and both relabelings are smaller) passes about
  one exchangeable 2-vs-2 null in six, and it called most readouts. Each row also carries the Figs. 5–6 rule: both
  relabelings ≤ 0.25 × the contrast. The ranking uses "strong".
- **Post-hoc additions.**
  - `injury_origin.csv`: where injury-induced genes sit relative to the healthy nephron.
  - Their S3 − S1 slopes.

## Answers

"Strong" means the four cross pairs agree in sign and both relabelings are ≤ 0.25 × the contrast. "Called" means
only the frozen rule passed. "Suggestive" means neither rule passed but the pattern is clear. Positions are fractions
of the PT (PT questions) or of the registered nephron axis (distal questions). The human values are from one donor
throughout.

### Healthy human vs mouse

| Question (brief) | Answer | Evidence |
|---|---|---|
| SGLT2-dominant zone (Q "does mouse overestimate it?") | **No, the opposite.** The human glucose-entry zone is about twice as long. SGLT2 stays above half its range to 0.66–0.71 of the human PT against 0.38–0.39 in mouse. SGLT1 starts at 0.70–0.76 against 0.35. GLUT2's boundary does not move (0.55–0.58 in both). | Strong (SGLT2 end +0.30, relabelings 0.02; SGLT1 start +0.38; switch +0.34). SGLT2 and GLUT2 are also "steeper in mouse" beyond the genome-wide flattening in the Fig. 5 screen (0 relabeling calls). |
| Drug-secretion domains | **Reorganised.** OAT1 peaks mid-PT in mouse (0.46) but at the start in human (0.13–0.16). OAT3 peaks at 0.24–0.28 against 0.06. OCT2 starts later in human (0.51–0.57 against 0.32–0.36). MRP2 starts at the PT origin in human but at 0.55–0.59 in mouse. In mouse, OAT1 and OAT3 vanish from S3; in human they persist. MATE1 is flat in both. | Strong: OAT1 peak, OAT3 peak, OCT2 start and MRP2 start. OAT1 and OAT3 are steeper in mouse beyond flattening (Fig. 5 screen). MRP4 extends further in human (end 0.71–0.94 against 0.22–0.26, strong), but the screen marked MRP4 not assessable, so it is not used. |
| Phosphate handoff (SLC34A1 → SLC34A3; the brief's #1) | **Suggestive.** In both species NaPi-IIc is an early-S1 gene. Its half-range domain ends at 0.26–0.30 in mouse and at 0.44–0.60 in human, so in human it overlaps NaPi-IIa over the first half of the PT. NaPi-IIa is flat along the human PT and falls in mouse S3. There is no clean switch point in human, because NaPi-IIa never falls. | Called, not strong (relabelings 0.10 and 0.06 against +0.24). SLC34A3 is not assessable in the Fig. 5 screen (low detection). |
| Urate | **URAT1 shifts late.** Its peak moves from mid-PT (0.54) to late PT (0.84). It keeps its amplitude in human while most genes flatten. PDZK1 is flat in both species, so its peak position means little. GLUT9 and ABCG2 are not in the PT curve set. | Strong (URAT1 peak +0.29, relabelings 0.00). This agrees with the Fig. 5 handoff SLC22A23 → URAT1 at +0.25. |
| Nutrient → secretion transition | **Not resolved.** The module crossing moves later in human (0.33 to 0.42 and 0.72), but the two human sections disagree. | Called, not strong. |
| TAL claudins and the TAL transition | **Conserved in cortex.** CLDN10 peaks before CLDN16/14 in both species, and their crossing is at 0.45 vs 0.49. UMOD and CASR occupy a wider share of human TAL, but human TAL is cortex only, so the corticomedullary question cannot be asked. | Strong for the UMOD and CASR start, but read as sampling: the human sections contain no medulla. |
| NCC–ENaC overlap | **Suggestive: no overlap zone in human.** NCC and ENaC-β/γ overlap over 12–19% of the DCT→CD window in mouse and 0% in both human sections. ENaC-β/γ start 0.05–0.07 later in human. This contradicts published human DCT2 overlap, and the feature is only 2–3 bins wide. | Strong by rule, but the axis is not gene-held-out and the bins are 0.017 wide. Needs a check before any claim. |
| DCT1 → DCT2, the Ca/Mg programs | **Mg and Ca programs extend further distally in human.** The PVALB domain ends at 0.76–0.78 against 0.64–0.68; the CALB1 domain at 0.92–0.94 against 0.81. The PVALB→CALB1 crossing does not move (0.62–0.64). | Strong for the PVALB peak and the CALB1 end, with the same axis caveat. |
| AQP2 / vasopressin onset | **Conserved.** AQP2 starts at 0.76–0.78 in both species, and the NCC→AQP2 crossing is at 0.74 in both. AVPR2 starts later in human (+0.08). | Strong for AVPR2, with the axis caveat. |
| Intercalated-cell order | **Inverted, but untrustworthy.** In mouse, pendrin peaks before AE1 (0.78 vs 0.85); in human the order is reversed (0.94 vs 0.72–0.74). Human CD is cortical only and IC cells are mixed into units. | Strong by rule; not interpretable. |
| Boundary genes | PVALB, CALB1 and SLC8A1 are segment-wide domains, not narrow transition-zone genes. | Descriptive. |
| Medullary osmoadaptation; OAT4, MATE2-K | **Not answerable.** There is no human medulla, and these genes are not in the ortholog panel. | — |

### Healthy vs AKI mouse (IR about 12 h–2 d)

| Question (brief) | Answer | Evidence |
|---|---|---|
| Does injury reuse the healthy coordinate or create new ones? (the brief's conceptual test) | **Mostly new.** 360 genes are induced in a PT label in all four IR–control pairs, and 0 under either relabeling. 252 (70%) exceed their level in every healthy segment, including KIM-1, KRT20, KRT8, HMOX1, LGALS3, SPP1 and CYP24A1. 89 (25%) reach a level seen in some healthy distal segment; their IMCD share (63%) is near background (53%). 19 (5%) stay within the healthy PT range. Of 278 repressed genes, 149 fall only to a level found elsewhere in the healthy PT; 129 fall below every healthy PT label. | Calibrated: 0 calls in both relabelings. Distal markers *fall* in injured PT units (UMOD −0.9 to −1.6, EGF −2.8 to −3.8), so contamination from neighbouring distal tubules does not explain the reused group. |
| Do injury genes have their own geometry? | **Yes: new S3 gradients.** Of the 329 induced genes that are flat in the healthy PT, 169 gain an S3 > S1 gradient in both IR mice and 10 an S1 > S3 one. Of the 31 induced genes with a healthy gradient, 26 keep its direction while being amplified (SPP1: a healthy early-S1 gradient, about 16× higher, still S1-tilted). | Post hoc. Both mice agree by construction of the count (each ≥ 0.5 log2). |
| Is there an ordered injury front (S3 → S2 → S1)? | **Yes, as nested domains rather than one moving front.** KIM-1 is induced along the whole PT and rises S1 3.8 → S2 5.3 → S3 6.4 log2. Half-range onsets in the two IR mice: KIM-1 0.46/0.46 → SOX9 0.50/0.57 → ICAM1 0.55/0.56 → VCAM1 0.58/0.65 → PCNA 0.62/0.67 → TOP2A, KRT8, KRT20, MKI67 0.64–0.72 → LCN2 0.70/0.74 → KRT18 0.80/0.84. The proliferation and stress-keratin onset sits at the control S2/S3 boundary (0.67). | Each gene's response passes its relabelings. The order agrees between the two IR mice within about 0.07 (IL34 differs by 0.13). Onsets are descriptive. |
| Is transporter loss centred where injury markers peak? | **No, it is offset.** Each segment loses its own specialised transport: SGLT2 −1.0 in S1 (0 in S3); NaPi-IIa −1.5/−1.8 in S1/S2 (0.1 in S3); OAT3 −1.8 in S1; AGT1 (SLC7A13) −3.6 in S3. KIM-1 peaks in S3, so S1 transport is lost where KIM-1 is lowest. Megalin and ACSM2 fall evenly (−1.1 to −1.7). | Each response passes its relabelings. Caveat: this is the gradient flattening of Fig. 6b, and blurred IR positions would produce the same pattern; genes whose level is unchanged also flatten (S2 retention 0.66–0.68). |
| In what order is mature transport lost? | **By segment, not by program.** There is no global order such as secretion before nutrients. The biggest loss is S3 amino-acid transport (AGT1 −3.6), then early nutrient and phosphate carriers in S1–S2, with endocytosis lost evenly. | As above. |
| Does metabolic failure precede, coincide with or extend beyond injury? | **It coincides with the proliferative domain and sits inside the KIM-1 field.** The ACADM loss onset is 0.59 (S3 −2.5 against S1 −0.6) and ACOX1's is 0.78. PPARA and CPT1A are *not* lost at this time point (PPARA +0.35 in S2). | Called (PPARA not called). |
| VCAM1: a distal extension of the KIM-1 domain? | **No.** VCAM1 and ICAM1 sit inside the KIM-1 field (VCAM1 peak 0.76 against 0.82). At 12 h–2 d there is no separate failed-repair domain. | Descriptive. |
| SPP1: ectopic? | **No, an amplified healthy gradient.** Healthy SPP1 falls along S1 (3.3 → 0.5). In IR it is high along the whole PT (7.6 → 6.2) and keeps the S1 tilt. This is the opposite end from KIM-1. | Strong response. |
| Common injury program across segments? | **Yes, weighted by segment.** LCN2, SPP1, CLU and KRT8 rise in every segment. LCN2 and CLU are largest in TAL and CD (+5 to +7.6) and smallest in PT (+1.8 to +3.4); KIM-1 is PT-weighted. | Strong in most segments. Distal KIM-1 may be PT admixture. |
| Is TAL injury position-specific? | **Yes.** NKCC2 falls −4.1 in the medullary third of TAL and −1.9 in the cortical third; ROMK −2.1 against −0.8. Both mice agree. LCN2 and SPP1 rise evenly (+5.3 to +6.0). Transport loss is localised; the injury markers are not. | Strong (relabeling peaks −0.4 and −0.2 against −4.7). |
| Distal compensation (NCC/ENaC up)? | **No.** NCC and ENaC are unchanged. DCT Ca/Mg genes fall (CALB1 −2.8, PVALB −1.9 in DCT2), and AQP2 and AVPR2 fall in CCD (−0.9, −1.2). | Strong for CALB1 and AVPR2. |
| DCT1→DCT2 boundary, ENaC/AQP2 expansion | **No shift beyond one bin.** The PVALB→CALB1 crossing moves −0.02 and the AQP2 onset −0.017 (bin 0.017). | Called; within resolution. |
| Repair-fate origin; inverted identity after repair | **Not answerable.** Single time point (about 12 h–2 d), no lineage tracing. | — |

### A correction for the Fig. 6 text

The Fig. 6 handoff OAT1 → OCT2 moves upstream (−0.14) on scaled curves. In absolute terms, OAT1 is *retained*
further into S3 after injury (+1.8 log2 in S3), while OCT2 loses its S3 peak (−1.2). The crossing moves because OCT2
collapses, not because OAT1 retreats. The Results sentence "switch points move upstream" should say which partner
moves.

## Ranking (most significant first)

The criterion was fixed in the protocol: strong before suggestive; then how directly the result answers a question
that segment clusters cannot; then translational weight.

1. **The human glucose-reabsorption zone is about twice as long as the mouse's.** SGLT2 extends to about 0.7 of the
   human PT against about 0.4 in mouse. SGLT1 begins correspondingly later. GLUT2's boundary does not move.
   - Strong and calibrated.
   - The gliflozin target, and the collaborators' SGLT2-inhibitor interest.
   - This answers the brief's question with the opposite of its premise.
2. **Injury builds new geometry rather than reusing healthy expression.** 70% of 360 induced genes exceed every
   healthy segment (0 under relabelings). Genes that are flat in health gain S3 > S1 gradients 17 times more often
   than the reverse. This is the brief's "warped versus new coordinates" test, answered genome-wide.
3. **Injury programs form nested, reproducible domains along the PT.**
   - KIM-1 marks a whole-PT field rising toward S3.
   - SOX9, ICAM1 and VCAM1 begin in mid-S2.
   - Proliferation and stress keratins switch on at the S2/S3 boundary.
   - KRT18 begins later still.
   - Both IR mice agree within about 0.07.
4. **The human drug-secretion apparatus is reorganised.**
   - Anion uptake (OAT1, OAT3) sits at the PT start and persists along the whole human PT.
   - In mouse, OAT1 peaks mid-PT, and both OATs vanish from S3.
   - OCT2 is confined to the later human PT; MRP2 starts at the human PT origin.
   - Strong. Relevant to drug clearance and nephrotoxicity.
5. **Transporter loss is offset from the injury markers.** Each segment loses its own specialised transport, so S1
   nutrient carriers are lost where KIM-1 is lowest. Caveat: this cannot be separated from axis blur.
6. **TAL injury is position-specific.** NKCC2 and ROMK loss is about twice as deep in medullary TAL, while LCN2 and
   SPP1 rise evenly. The loop-diuretic target.
7. **Metabolic failure (ACADM, ACOX1) coincides with the proliferative S3 domain, without PPARα or CPT1A loss.**
8. **URAT1 moves to the late human PT and keeps its gradient** (strong). The phosphate result, NaPi-IIc extending
   into the human mid-PT, is *suggestive* only. It would be the anatomical correlate of the brief's #1, but it does
   not pass the strict rule.
9. **SPP1 is an amplified healthy S1 gradient, not an ectopic domain**, and it sits opposite KIM-1.
10. **Distal nephron.** No NCC/ENaC compensation after AKI; Ca/Mg and water genes fall. In human, the Mg and Ca
    programs extend further distally and NCC–ENaC do not overlap. These rest on a non-held-out axis at 2–3 bin
    resolution, so they need a check before any claim.

**Not answerable here:**
- medullary osmoadaptation and the human corticomedullary TAL (no human medulla);
- OAT4 and MATE2-K (no ortholog);
- repair-fate origins and post-repair identity (one early time point);
- intercalated-cell order (cortical CD only, mixed units).
