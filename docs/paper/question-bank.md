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
