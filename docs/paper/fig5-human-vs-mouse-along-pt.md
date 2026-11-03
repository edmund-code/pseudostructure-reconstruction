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
