# Borrowing time: which injury phases our IR mice carry, and where along the PT each acts (notebook 99)

**The question.** Our IR mice were sampled at one unknown time after injury (`pt-package-issues.md`, item 10). Notebook
94 aligned them best with 12 h of Kirita's time course (ρ 0.52), but the cell-cycle program pointed to 2 d. Kirita
2020 (main text, read in full, PMC7355049) describes three phases of the PT response:
- **acute:** "'Injured S1/2' and 'injured S3' were primarily composed of cells from the 4- and 12-h time points";
- **repair:** a repairing, mitotic PT state arising at 2 d (Dataset S2 description, as recorded in notebook 80);
- **late:** failed-repair PT that "emerged … at 2 d after IRI and then expanded and remained … at late time points".

Ischemic injury "is most apparent in the S3 segment" (Bonventre & Yang 2011, read in full).

This notebook uses Kirita's time course as a ruler. It defines acute, repair, late and sustained gene programs from
Kirita's mice alone, then asks:
1. which of these phases our IR mice carry;
2. where along the PT each phase acts.

A single-time-point spatial dataset cannot show this by itself.

Protocol frozen on 2026-10-18, before any phase program was scored on our units. Notebook 94's external gene table was
known, since it is its result. No phase program had been defined or intersected with our data.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/99_pt_injury_phases.ipynb`; logic `99.phases.1`; seed 99.
- Results `results/pt_injury_phases/`.
- Module: `pseudospace/tissue_location.py` (`program_score`, `position_expectation`), plus one new helper there,
  `gene_effects`, with tests added to `tests/test_tissue_location.py`.
- **Kirita statistics:** notebook 94's `external_genes.csv.gz`. These are the centred whole-PT level coefficients
  (log2) and q at 4 h, 12 h, 2 d, 14 d and 6 wk, over notebook 94's E.
- **Our units:** the AKI export through `support_positions` (notebook 96's construction: 14,048 common-support units,
  positions x ∈ [0, 1], package segment labels). The controls' segment cuts c12 and c23 come from `segment_cuts`.

## 2 · Phase programs (from Kirita's data only; fixed now)

Per gene and time t, *up_t* means q ≤ 0.05 and coefficient ≥ 1, and *down_t* means q ≤ 0.05 and coefficient ≤ −1.
"Quiet at t" means coefficient < 0.5 for up programs and > −0.5 for down programs.

| Program | Up version | Down version |
|---|---|---|
| ACUTE | up at 4 h or 12 h; quiet at 2 d and 14 d | the same with down |
| REPAIR | up at 2 d; quiet at 4 h and 12 h | the same with down |
| LATE | up at 14 d; quiet at 4 h, 12 h and 2 d | the same with down |
| SUSTAINED | up at 12 h, 2 d and 14 d | the same with down |

- Programs keep members detected in ≥ 2% of our PT units in both control mice, or in both IR mice (91's rule).
- They are scored only with ≥ 10 such members.
- 6 wk is not used, because Kirita has three level changes at 6 wk.

## 3 · Per-gene injury effects (our units)

- **Gene z per unit:** (log1p CP10k − control mean) / max(control SD, 0.25), with the mean and SD over the pooled
  control units. This is notebook 96's program score with one gene.
- **Gene effect:** z minus the control expectation at the unit's position. The expectation is the binned mean, 0.05
  wide, interpolated; control units are cross-fitted from the other control.
- **Per gene and mouse:**
  - *level* = the mean effect over the mouse's units;
  - *S3−S1* = the mean effect in units at x ≥ c23 minus the mean in units at x < c12;
  - the binned effect profile, 20 bins of x.
- A program's statistics are the means of its members' statistics. They are linear, so they equal the statistics of
  the program score.

## 4 · Analyses and decision rules

**A1 · Which phases do our IR mice carry?**
- For each program and IR mouse: level z = (program level − null mean) / null SD. The null is 1,000 random programs of
  the same size, drawn from the scored genes and matched member-for-member on control-expression decile.
- **A phase is present** when the up and the down versions both have z ≥ 3 in the expected direction (up positive,
  down negative) in both IR mice. Each version is also reported alone.
- The same statistics in the two control mice (cross-fitted) give the reference. They should be about 0.
- **Time mixture (descriptive):**
  - non-negative least squares of our gene L (notebook 91, condition) on Kirita's five level coefficients, over the
    genes in both;
  - the weights, normalised to sum to 1, and the fit's R².

**A2 · Where along the PT does each present phase act?**
- For each present program and IR mouse: S3−S1 z against the same matched null, and the peak bin of its profile.
- **H1 (acute injury is S3-weighted):** S3−S1 z ≥ 3 in both IR mice for ACUTE-up, and ≤ −3 for ACUTE-down.
- **H2 (repair is S3-weighted):** the same for REPAIR-up and REPAIR-down.
- **H3 (phases differ in location; descriptive):** the S3−S1 contrast divided by the level, per program, so that a
  larger program does not look more located. ACUTE against REPAIR against LATE.
- A program that is not present is reported, not judged.

## 5 · Outputs and figures

- **Tables:**
  - `phase_programs.csv` (members, Kirita coefficients);
  - `gene_effects.csv.gz`;
  - `a1_presence.csv`, `a1_time_mixture.csv`, `a2_location.csv`, `profiles.csv`;
  - `verdicts.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 99a: Kirita's phase programs over the time course, as mean coefficient against time (the definitions);
  - 99b: our IR mice's level z per program, against controls;
  - 99c: each present program's effect profile along x, per IR mouse, with the S1/S2/S3 cuts;
  - 99d: location per level by phase.

## 6 · Caveats

- Kirita profiles nuclei from bilateral 18-min IRI in male C57BL/6J mice. Our surgery is unknown. Phases are defined
  by their timing in that model.
- A phase program can contain genes that respond at several times. The quiet rules make the programs specific, at the
  cost of size.
- The deepest injured S3 units are under-represented in our export (notebook 76). Notebook 96's placed-unit arm did
  not change its location verdicts, and this notebook does not repeat it.

## Notebook 99 results

`analysis/notebooks/99_pt_injury_phases.ipynb` (logic `99.phases.1`) writes to `results/pt_injury_phases/`. It
followed the protocol frozen in 1a92c73. G1 passed: 14,048 units, with control cuts 0.343 and 0.647. New helpers are
`tissue_location.gene_effects` and `external_replication.matched_means`, each with a test. The `gene_effects` test
checks it against notebook 96's per-program path.

**How the protocol was read.**
- Kirita's coefficients and q are notebook 94's.
- The time mixture has no intercept.
- The deciles come from the controls' mean log1p CP10k.
- The profile peak is the bin with the largest signed effect.
- One post hoc section (5.1).

**Phase programs (Kirita only; members scored among 10,683 genes):**

| Program | Members | Mean log2 FC (4 h, 12 h, 2 d, 14 d, 6 wk) |
|---|---|---|
| ACUTE-up | 320 | 1.16, 1.27, −0.22, −0.08, −0.07 |
| REPAIR-up | 223 | −0.36, −0.30, 1.97, 0.62, 0.45 |
| LATE-up | 16 | −0.45, −0.34, −0.01, 2.02, 0.77 |
| SUSTAINED-up | 257 | 2.51, 3.06, 2.79, 2.81, 1.54 |
| ACUTE-down | 473 | −1.16, −1.63, 0.23, 0.04, 0.10 |
| REPAIR-down | 19 | 0.33, 0.16, −2.05, −0.74, −0.79 |
| SUSTAINED-down | 48 | −2.36, −3.04, −2.45, −1.72, −0.70 |

LATE-down (1 gene) is not scored.

**A1 · Which phases our IR mice carry** (level z against matched random programs; IR2A2, IR2A4 | controls):

| Program | IR | Controls |
|---|---|---|
| ACUTE-up / -down | 5.7, 4.2 / −9.0, −7.4 | 0.1, 0.0 / ±4.7 |
| REPAIR-up / -down | 2.4, 3.6 / −2.9, −3.2 | ±1.0 / ±0.8 |
| LATE-up | −1.0, −1.1 | ±3.3 |
| SUSTAINED-up / -down | 26.1, 24.4 / −8.6, −8.7 | ±10.1 / ±1.7 |

- **Present:** ACUTE and SUSTAINED.
- **Not present:** REPAIR (IR2A2's up version, at 2.4, misses the threshold of 3) and LATE.
- **Time mixture** (non-negative least squares of our L on Kirita's five time points; R² 0.38 over 10,017 genes):
  12 h 0.57, 2 d 0.43, nothing at 4 h, 14 d or 6 wk.
- **Our IR mice sit between 12 h and 2 d after IRI.** They carry the acute and sustained responses fully, the 2-day
  repair response partly, and none of the late or failed-repair response.

**A2 · Where along the PT each phase acts** (S3−S1 z; IR2A2, IR2A4):

| Program | S3−S1 z | Per unit of level | Peak x |
|---|---|---|---|
| ACUTE-up | **7.9, 4.5** | 1.75 (IR2A4) | 0.93, 0.88 |
| ACUTE-down | **−6.2, −3.5** | – | 0.98, 0.98 |
| SUSTAINED-up | 20.4, 18.8 | 0.98 | 0.88, 0.88 |
| SUSTAINED-down | −6.0, −7.1 | 0.61 | 0.83, 0.88 |
| REPAIR-up (not present) | 2.7, 6.9 | 2.16 (IR2A4) | 0.78, 0.88 |
| REPAIR-down (not present) | −2.4, −3.3 | – | 0.78, 0.78 |

- **H1 holds.** Acute injury is S3-weighted in both IR mice, in both directions.
- **H2 is not judged,** because repair is not present. Descriptively, the partial repair response is the most
  S3-weighted per unit of level (IR2A4 2.16, against 1.75 for acute).
- **The up programs share the shape of notebook 91's induced program,** an observation not specified in advance:
  - raised in early S1;
  - a dip in late S1 around x ≈ 0.3, the same dip 91 found for its induced program;
  - a rise through S2 into late S3.

**Post hoc 5.1: IR against the control mice's own spread.**
- The cross-fitted controls are not at 0, because of the A2/A4 factor: SUSTAINED-up level control z ±10.1, ACUTE S3−S1
  control z ±3.7.
- Against the larger control |z|, these still put both IR mice beyond the controls:
  - ACUTE (up and down), level and S3−S1;
  - SUSTAINED (up and down), level and S3−S1;
  - REPAIR, level only.
- H1 survives this calibration. IR2A4's ACUTE-down S3−S1 is the narrowest margin: −3.54 against 3.45.

**What this answers.**
1. **Time point.** Our IR mice match 12 h to 2 d after an 18-min bilateral IRI: the acute and sustained programs, part
   of the 2-day repair program, nothing later. This narrows `pt-package-issues.md` item 10 from "unknown" to an early
   repair window. The cell-cycle program's 2-day signature (notebook 94, K2) fits the repair share.
2. **Acute injury is S3-weighted along the continuous axis,** in both directions and both mice, beyond the controls'
   spread. This is the textbook S3 vulnerability (Bonventre & Yang), measured gene-program by gene-program on a
   continuous coordinate, with the phase defined in an independent experiment.
3. **The late-S1 dip is a feature of every injury phase,** not of one program. Late S1 is the least responsive
   stretch of the PT in both directions, as notebook 91 found.
4. **As a method,** an external time course used as a ruler turns a single-time-point spatial dataset into a phase
   decomposition, with location per phase. It needs the external experiment and the matched null, and it should be
   read against the controls' own spread.
