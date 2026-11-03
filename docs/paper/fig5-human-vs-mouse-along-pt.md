# Fig. 5: mouse against human along the proximal tubule — protocol and results

Fig. 5 is the paper's summary of the human vs mouse comparison. The abstract claims "pathways and paired solute
transporters whose proximal tubule positions differ between mouse and human". This page freezes the screen that
backs that sentence on the paper's machinery (the Fig. 4 and Fig. 6 axes), and it reuses what earlier rounds already
established instead of re-deriving it:

- notebook 101's transporter pairs (`docs/results/pt-transporter-handoffs.md`), frozen before any species comparison;
- the species confounds measured in `human-mouse-paper-ready.md` and `hvm-continuous-discovery-round.md`: probe
  panels differ, so levels are never compared; the human arm is flatter (median amplitude ratio 0.63, cause not
  separable); human S1/S2 labels are unreliable (most human S2-zone units are labelled S1); human within-zone order is
  not trustworthy (±0.06); the human tissue is one donor, two cortex sections;
- the per-gene literature reviews and artifact checks of the 2026-10-21 to 10-27 rounds (gene agents, position-only
  report, pathway ground truth, lactate), used to interpret calls after the screen. New literature is read only for
  calls those rounds did not cover.

## Frozen protocol (2026-11-03, before any fit)

### 1 · Data and axis

- **Fit.** The human vs mouse package fit (`runs/hvm`, notebook 70), units with a reference label PT-S1, PT-S2 or
  PT-S3 in Ctrl1A2, Ctrl1A4, HUK1_COR1 and HUK1_MED1. Raw counts from `02_qc/units_qc.h5ad` (1:1 ortholog space,
  mouse symbols); library = total over non-mitochondrial genes.
- **Axis.** As Fig. 4 and Fig. 6: each gene's position for a unit is the mean of the five Fig. 3d gene-half axes of
  this fit built without that gene (`fig3/robustness/hvm/reruns`); a unit needs 3 of 5.
- **Scale x: mouse quantile.** A mouse unit's x is its rank among the PT units of both mice; a human unit's mean
  position is mapped through the mouse units' empirical distribution and clipped to [0, 1]. The clipped share of
  each human section is reported. The mouse sections are the Fig. 4 tubules, so x means the same thing in Figs 4–6.
- **Boundaries and cores** from the mouse units only (Fig. 4's minimum-misclassification rule; cores > 0.06 from
  b12 and b23). Human labels are not used to place boundaries, and no human segment is named; the segment band in
  every panel is the mouse's.
- **Knots** at mouse x quantiles; the same basis in all four sections.

### 2 · Contrasts

For any per-section quantity Q:
- **species:** mean(Q of HUK1_COR1, HUK1_MED1) − mean(Q of Ctrl1A2, Ctrl1A4);
- **the four pairs:** Q(human section) − Q(mouse section);
- **relabeling 1:** mean(COR1, Ctrl1A2) − mean(MED1, Ctrl1A4); **relabeling 2:** mean(COR1, Ctrl1A4) −
  mean(MED1, Ctrl1A2).
- **Specific:** both relabelings' |value| ≤ 0.25 × the species |value|.
- The relabelings measure section-to-section variation. Both human sections come from one donor, so they cannot
  measure donor-to-donor variation; every call is a statement about this donor.

### 3 · Genes and models

- Mitochondrial genes are left out. **Screened:** detected in ≥ 10% of PT units in each of the four sections (a
  shape comparison needs the gene in both species).
- **Neighbour filter:** Fig. 4's rule (dropped if AL, DCT or CD is more than twice the PT rate), within each
  species' sections; failing in either species drops the gene.
- **Models** per gene and section, as Fig. 6: NB GLMs with a log(library) offset; curve (6-df spline), clusters
  (label factor; descriptive only, human labels are unreliable) and bulk; dispersion under clusters + curve; curves
  on the 201-point grid in log2 counts per 10k; 20 equal-count bin means.
- **Unsupported:** a gene with no count at all in the first or last third of the S1 or S3 core, in any section, is
  not assessed for gradients (Fig. 6's post-hoc rule, adopted here before any fit).

### 4 · Gradients and their comparison

- **E (per section):** mean of the curve over the S3 core minus its mean over the S1 core (log2).
- **Raw R (per section):** log2 rate (summed counts + 0.5, per 10k) over S3-labelled units inside the S3 core
  minus the same over S1-labelled units inside the S1 core. No spline, no registration inside a core.
- **Support:** each of those two unit sets holds ≥ 25 units in every section.

**Category R — the gradient reverses between species.** All of:
1. both mice share the sign of E and both human sections have the opposite sign; |E| ≥ 0.5 in all four;
2. raw R has the sign of that section's E in all four sections;
3. support, and not unsupported.

**Category H — steeper in human** (attenuation cannot produce it). All of:
1. E has the same sign in all four sections; |E| ≥ 1 in both human sections;
2. |E(human section)| − |E(mouse section)| ≥ 0.5 for all four pairs;
3. raw R has that sign in all four sections and |R(human)| > |R(mouse)| for all four pairs;
4. support, and not unsupported.

The mirror category (steeper in mouse) is only counted: it cannot be told apart from the human arm's general
flattening.

- **Calibration.** Each category is re-applied with the species grouping replaced by each relabeling's grouping
  (arm 1 = the plus sections, arm 2 = the minus sections). False-call share = mean relabeling calls / species calls.
  Above 0.10, the category is reported as not calibrated and the figure shows no list from it.
- **Examples** are ranked by the smallest |E| over the four sections (R) or the smallest pair difference (H).

**Beyond flattening (for sets).** For each human × mouse pair, OLS of E(human section) on E(mouse section) over the
screened, supported genes gives the genome-wide flattening line; a gene's pair residual is its distance from that
line. The species residual d is the mean of the four pair residuals. Relabelings get the same construction on their
own four cross pairs. The four slopes are reported (the flattening itself).

### 5 · Pathways

- Sets: `pt_package_aki/tables/gene_sets.json` (Hallmark, Reactome, KEGG mouse), as Fig. 6; ≥ 10 screened,
  supported members.
- **Statistic:** the members' mean d. Negative: the set sits relatively earlier along the human PT than the mouse
  gradient predicts after flattening; positive: later.
- **Null:** 1,000 random size-matched sets of screened, supported genes; two-sided empirical p; BH within each
  split.
- **Call:** q ≤ 0.05, all four pair means share the sign, specific.
- **Calibration:** mean relabeling calls / species calls ≤ 0.10, else not calibrated.
- **Display (no claim):** per section, the members' mean curve after centring each gene on its own mean (levels are
  not comparable across probe panels).

### 6 · Transporter handoffs

- **Pairs (frozen elsewhere):** notebook 101's mouse pairs confirmed by Chen 2021 (50) and human pairs confirmed by
  Census (16), deduplicated (57), plus the poster's lactate (Slc5a12 → Slc5a8) and glucose (Slc5a2 → Slc5a1) pairs,
  marked as previously shown. Members outside the screened universe are fitted the same way for this purpose only.
- **Profiles and crossing x\*:** as Fig. 6 (`pseudospace.handoffs`; each curve min–max scaled over the grid points
  inside that section's 2.5–97.5% range of x). Min–max scaling makes x\* invariant to a per-gene, per-species offset
  and scale of the log curve, so neither probe efficiency nor flattening moves it.
- **Shifted:** a crossing in all four sections; all four pair Δx\* share a sign; |Δx\*| ≥ 0.05; specific; the 95%
  interval of the species Δx\* from 200 bootstraps (units resampled within section, curves refitted at fixed
  dispersion) excludes 0.
- **Mouse-only:** crossings in both mice and in neither human section. **Human-only:** the reverse.
- **Warp check (descriptive).** The shared axis registers human onto mouse through integration, so a smooth warp
  between species would move every crossing. For each section, the switch point of every screened gene with
  |E| ≥ 0.5 is where its scaled curve first crosses 0.5 in its gradient's direction; a pair's x\* is re-expressed
  as its quantile among those switch points. A shifted pair whose quantile shift has the same sign in all four pairs
  is reported as moving beyond the common warp.
- **Levels are not shown** (probe panels differ).

### 7 · Gradient conservation (descriptive)

Fig. 6's retention, across species: genes with |E| ≥ 0.75 in one mouse; the median signed E of the other mouse
and of each human section over the median selected |E|. The same selected in one human section and read in the
other human section and each mouse. Selecting in one section and reading another avoids regression to the mean.

### 8 · Named checks for anything shown in the main figure (after the screen)

- **Neighbour admixture:** the gene's neighbour ratio in each species, and the correlation of its bin means with a
  distal-identity score (AL, DCT, CD marker rate) along x in each section; admixture predicts a positive one.
- **Human cell state:** human units carry ~6× the mouse mitochondrial fraction and the early human axis tracks it.
  Raw R of each called gene and set recomputed on the lowest-mito tertile of each human section.
- **Literature:** the existing reviews first; new reading only where none exists.

### 9 · Figure plan

Main: reversed gradients; pathway relocations; handoffs in the poster layout (mouse and human side by side).
Extended Data: maps, coverage, gradient conservation, steeper-in-human genes, pathway map, all moving handoffs and
the checks. Which examples are drawn is a presentation choice made after the screen and recorded here.

## Results (2026-11-03)

`analysis/paper/fig5_screen.py` (protocol as frozen in 367e917; about 45 min at 24 jobs), then
`fig5_interpret.py` (post-hoc checks) and `fig5_panels.py`. Outputs in
`results/paper_figures/fig5/{screen,panels}/`; tests `tests/test_paper_fig5.py`.

**How the protocol was read.**
- E is averaged over the core's grid points inside the section's own 2.5–97.5% x range, so no curve is read
  where a section has no tubules (the protocol said "over the core").
- As in Fig. 6, every handoff pair uses the x range of the first handoff gene's axis; the ranges of the 32 axes
  differ by < 0.01.
- Handoff partners outside the screened universe were fitted for the handoffs only (e.g. Slc5a1, Slc5a8, which are
  detected in 6–7% of human PT units).

**Data.**
- PT units: Ctrl1A2 2,644, Ctrl1A4 2,650, HUK1_COR1 3,194, HUK1_MED1 3,459. Median mitochondrial fraction:
  mouse 0.04/0.035, human 0.18/0.19.
- 5,983 genes are detected in ≥ 10% of PT units in all four sections; the neighbour rule drops 2,052 of them (more
  than in the AKI fit) and 59 did not converge, leaving 3,872 screened; 3,869 are assessable.
- Boundaries over the 32 axes: b12 0.318–0.340, b23 0.633–0.660. Up to 3.2% of human units sit before the first
  mouse unit and are clipped to x = 0; ≤ 0.13% at the end.

**The flattening.** OLS of human E on mouse E over the assessable genes has slope 0.09–0.12 in the four pairs
(regression dilution included). Gradient conservation (§7):

| chosen in | genes | other mouse keeps | human sections keep | same sign in human |
|---|---|---|---|---|
| Ctrl1A2 | 516 | 1.18 | 0.05 / 0.09 | 0.59 / 0.62 |
| Ctrl1A4 | 682 | 0.84 | 0.05 / 0.09 | 0.60 / 0.64 |

| chosen in | genes | other human keeps | mice keep | same sign in mouse |
|---|---|---|---|---|
| HUK1_COR1 | 95 | 1.04 | 0.96 / 0.71 | 0.72 / 0.72 |
| HUK1_MED1 | 142 | 0.88 | 0.64 / 0.59 | 0.73 / 0.71 |

Each species reproduces its own gradients. The mice keep most of a human gradient, but the human sections keep
under a tenth of a mouse gradient: human PT gradients are, on this platform and donor, a much flatter version of
the mouse ones. The cause (probe capture, donor, biology) is not separable here (`human-arm-attenuation`), so the
figure compares directions and positions, never amplitudes.

**Category R — reversed gradients: 28 calls, 0 and 0 under the relabelings; calibrated.**
Apoc3, Hpgd, Tmem252, Vnn1, Dcxr, Aldh3b1, Lpar6, Abhd14a, Sh3gl2, Frmd4b, Pitpnc1, Mme, Slc5a11, Slc25a15,
Slc17a3, Sorbs1, Msrb1, Rai14, Lap3, Ctsh, Fuca2, Ugt3a1, Chrna4, Hao2, Lrig1, Odc1, Tex264, Reep6 (by smallest
|E|). Earlier leads not called: Slc10a2 (human MED1 −0.30, under the 0.5 floor), Miox (MED1 +0.49), Acot12
(MED1 −0.14), Pck1 (mouse −0.54/−0.39 is the wrong limb for the floor), Dpp4 (Ctrl1A2 +0.23).

**Category H — steeper in human: 7 calls, 0 and 0; calibrated.** Gpt, Pter, Sel1l3, Slc7a13, Ggact, Fut8, Fmo5.
The mirror category (steeper in mouse) has 131 genes and 0/0 under the relabelings; it is not called.

**Pathways: no call.** 890 sets; 0 species calls and 0/0 under the relabelings. The species p distribution does
not depart from the relabelings' (p ≤ 0.01: 17 sets, against 21 and 20; ED panel f). The best-ranked sets all sit
earlier along the human PT than the trend predicts: Lysosome (KEGG; d −0.22, p 0.001, q 0.22), MHC class II
antigen presentation, Coagulation, Metabolism of vitamins and cofactors. Lysosome and MHC II are led by the
cathepsins and legumain (Lgmn, Ctsl, Ctsh, Ctsb, Ctsf) and survive raw counts, low-mito human units, cortical mouse
S3 and dropping their top three genes (p 0.002, 0.008). This matches the 2026-10-26 relocation screen, but this
screen cannot separate it from section-to-section variation, so it is not a finding.

**Transporter handoffs: 13 of 53 shift, all 13 beyond the common warp; the same rule applied to each relabeling (from the saved bootstraps, after the run) calls 0 and 0.**
59 pairs; switches in both mice 58, in both human sections 54; 2 pairs switch in mouse only (MCT12 → MCT5, MCT12 →
MCT9), none in human only.

| switch later in human | Δx\* (95% interval) | switch earlier in human | Δx\* |
|---|---|---|---|
| SGLT2 → SGLT1 (poster pair) | +0.30 (0.26, 0.32) | MCT1 → MCT9 | −0.49 (−0.53, −0.46) |
| SGLT2 → SMIT2 | +0.26 (0.16, 0.45) | MCT10 → MCT9 | −0.33 (−0.39, −0.17) |
| SLC22A23 → URAT1 | +0.25 (0.08, 0.31) | MCT1 → MCT5 | −0.32 (−0.48, −0.00) |
| y+LAT1 → AGT1 | +0.24 (0.20, 0.26) | OCTN2 → URAT1 | −0.22 (−0.26, −0.02) |
| SLC22A23 → OCT2 | +0.24 (0.11, 0.32) | OCTN2 → OAT2 | −0.16 (−0.33, −0.09) |
| SGLT2 → SGLT5 | +0.21 (0.11, 0.27) | | |
| LAT2 → AGT1 | +0.19 (0.13, 0.22) | | |
| B0AT1 → B0AT3 | +0.18 (0.13, 0.21) | | |

- **What moves them.** In the glucose and amino-acid pairs, the early partner is far flatter in human (SGLT2 E
  −1.4/−1.5 against −4.5/−8.0; B0AT1 −0.6/−1.0 against −4.6/−6.2; y+LAT1 −1.5 against −3.0/−3.8), so it stays on
  through the mouse S2 and the late partner takes over at the mouse S2/S3 boundary. In the MCT pairs the late
  partners reverse: MCT9 and MCT5 rise toward S3 in mouse (E +1.3/+1.5, +1.2/+1.3) and fall in human
  (−0.4/−0.6, −0.2/−0.1), so the human "switch" is two falling curves crossing early. MCT9 is a reversal, not a
  moved handoff.
- **Not shifted:** lactate SMCT2 → SMCT1, +0.40 but its interval reaches −0.01 (human SMCT2 is flat, so the
  crossing is unstable); cystine b0,+AT → AGT1, +0.14, not specific. The poster's "no handoff in human" for lactate
  and B0AT1 does not hold on this axis: both cross, B0AT1 → B0AT3 later.

**Checks (post hoc where marked).**
- **Admixture (reformulated post hoc).** The named ρ with a distal-identity score turned out uninformative: the
  distal signal rises along the PT in both species, so ρ restates each gene's direction. The replacement is a bound.
  The S3 core carries thick-limb signal at 0.14 / 0.05 of a thick-limb tubule's (mouse) and 0.08 / 0.08 (human);
  the S1 core at most 0.05. With a neighbour ratio ≤ 2 (the screen's rule), admixture can add at most 0.22 log2 to
  any called gene's gradient, against a smallest |E| of 0.5, and it can only raise late-PT expression. No call is
  flagged.
- **Human cell state.** Raw R of every call keeps its sign in the lowest-mito tertile of both human sections
  (35/35).
- **Cortex-only human (post hoc).** Restricting mouse S3 to cortical units (8–13% of S3-core units) keeps the sign
  for 32/35 calls; Lrig1, Hao2 (reversed) and Sel1l3 (steeper) lose it.

**Abstract wording.** "We identify pathways and paired solute transporters whose proximal tubule positions differ"
is not supported for pathways: the calibrated set screen calls none. Supported: "genes and paired solute
transporters".
