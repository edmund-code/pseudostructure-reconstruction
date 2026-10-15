# Six ways to use the continuous PT axis, from the literature review (notebook 90)

Frozen on 2026-10-14, before any statistic of notebook 90 was computed. The user asked, after the literature review
(`docs/research/continuous-analysis-literature.md`): "run all the experiments and send me a single figure for each.
the goal is get the best results to highlight in the paper for why pseudostructure is useful. this doesn't have to be
statistical, can be descriptive. but can also be statistical. we want to fully use the continuous nature of the
pseudostructure, expand upon known biology or discover things that might be missed by current analyses."

"All the experiments" are the review's recommendations 1–5 and the IR part of recommendation 4. Recommendation 6
(MEFISTO) is not run: `mofapy2` is not installed, it would add a dependency for an exploratory method with no test and
no count model, and notebook 89's E1 already learned smooth programs and replicated them between the mice.

**What I had seen before freezing.**
- The results of notebooks 74, 76, 84, 85 and 89; in particular C14 (the curve's fold noise exceeds the
  within-segment spread in mouse S1 and S2), 85's resolution (a label switch is about 0.3 mouse units wide on the
  curve; a 6-df spline draws a step at S2/S3 about 1 unit wide), 89's two opposite pathway gradients and its absence
  of boundary programs, and 76's AKI verdicts.
- Design facts only, from 89's loaders: mouse U_PT units in the 0.1-unit bins on [0, 2.5) are 1,737 (Ctrl1A2) and
  1,844 (Ctrl1A4), 22–152 per bin; 9,487 usable PT-expressed mouse genes, 590 of them among the 1,105 genes that built
  the coordinate or the labels; 1,404 of 71's sets have ≥ 10 usable members (1,356 with ≥ 10 outside the built
  genes); the human sections hold 2,625 and 2,636 units, and HUK1_MED1's transitions cross (0.667, 0.590).

## Labels

- Notebook: `analysis/notebooks/90_pt_continuous_axis_experiments.ipynb`
- Logic version: `90.continuous_axis.1`
- Results: `results/pt_continuous_axis_experiments/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`)
- New module: `pseudospace/continuous_axis.py`, with `tests/test_continuous_axis.py`.

## Common inputs (89's, unchanged)

- The hvm export and 71's inputs (G1), 72's transitions (G2), the stage-01 native panels (G3), 86's mouse NB
  dispersions (G5), the AKI export (G7). The guards are 89's.
- **Units:** U_PT (package `primary_identity` Proximal_Tubule), so the thin-limb tail of late S3 (85, PH2) is out.
- **Mouse axis:** mouse segment units (the mean of the two mice's logistic transitions), 25 fixed 0.1-unit bins on
  [0, 2.5).
- **Profiles:** per mouse, each gene's binned log2 rate minus its mean, minus the reference (the median usable gene
  per bin): 89's deviation profile. Human: 25 equal bins over each section's 1–99% range (89), descriptive.
- **Genes:** usable PT-expressed genes (detection ≥ 5% in both mice; 86's rule).
- **Fold protection:** genes outside the 1,105 built genes (integration features of both passes and the 15 label and
  orientation genes; 74's rule). Every readout is reported for all usable genes and for the fold-protected genes;
  any statement that rests on the axis organising genes rests on the fold-protected genes.
- **Thresholds:** effect 0.25 log2 (the project's), BIC margin 6 (85).

## A · Gene shape classes, with segment steps as a competing model (review recommendation 1)

**Question.** For each gene, is its profile along the PT better described by steps at the segment labels (what an
S1/S2/S3 analysis assumes) or by a continuous shape: a gradient, a switch at a free position, or a transient bump?

**Models.** Weighted least squares on the two mice's deviation profiles together (50 points; weights = units per
bin, scaled to mean 1; one intercept per mouse):

| Model | Shape | Parameters |
|---|---|---|
| F | flat | 2 |
| L | linear gradient | 3 |
| S | logistic switch: midpoint m ∈ [0.15, 2.35], 10–90% width w ∈ [0.1, 4] | 5 |
| B | Gaussian bump or trough: centre c ∈ [0.15, 2.35], FWHM ∈ [0.2, 1.2] | 5 |
| D | discrete: constant within each label; the expected bin value is Σ_l f_bl·β_l, with f_bl the mouse's share of its units in bin b labelled l | 4 |
| DS | discrete plus a slope within each label | 7 |

Nonlinear parameters are fitted on a grid (m and c 45 values, w and FWHM 25 geometric values), amplitudes by least
squares. BIC = n·ln(WRSS/n) + k·ln n.

**Classes.** The model with the lowest BIC; "flat" when that is F or when the fitted curve's range is < 0.25 log2
(mean of the mice). S is split by location (within 0.15 of a boundary: "at S1/S2", "at S2/S3"; else "inside Sk") and
width (sharp ≤ 0.5, graded, gradual ≥ 1; 85's). B by the segment of its peak and its sign.

**Continuous against discrete.** ΔBIC = BIC(D) − min(BIC(L), BIC(S), BIC(B)), on non-flat genes. Above 6: continuous
favoured; below −6: discrete favoured; else undecided.

**Replication.** Each mouse fitted alone (25 points, one intercept); the share of non-flat genes with the same class in
both mice, against the agreement expected from the two mice's class frequencies.

**Semi-synthetic calibration (the comparator for the continuous share).** 400 non-flat fold-protected genes
(seed 90). Per gene and mouse, NB counts on the real units with the gene's 86 dispersion and the real libraries
(89's `simulate_counts`) under (i) its fitted D model, units taking their own label's level, and (ii) its best
continuous model at the unit's position. Binned, classified as the data. Readouts: the share of D-truth genes called
continuous (ΔBIC > 6) and of continuous-truth genes called discrete (ΔBIC < −6).

**References (descriptive).** Ransick 2019's segment markers (Slc5a2, Spp2 S1; Slc34a1 S1–S2; Cyp2e1 S2; Slc22a7,
Cyp7b1 S3; the male S2/S3 genes form "a tight, non-overlapping boundary" by RNAscope). Several are label or
orientation genes, so they show what the classifier does, not evidence.

**Highlight criteria (descriptive).**
- Continuous share: the share of non-flat fold-protected genes with ΔBIC > 6, if it exceeds the D-truth
  false-continuous share by at least 2-fold.
- Bumps: genes classified B in the joint fit and in both single-mouse fits, with peaks within 0.3 units of each
  other.
- Switch positions: where switches sit, inside segments against at boundaries.

**Figure 90a.** Class counts (with the replicated share); the genes in class blocks sorted by position; ΔBIC for the
data against the two simulated truths; one exemplar gene per class (both mice, best fit, D fit).

## B · Where pathways sit, and pathways a median curve flattens (recommendation 3)

**Set positions.** Per usable gene and mouse, the centre of mass CoM = Σ_b u_b·2^{d_b} / Σ_b 2^{d_b} (u_b bin centres,
d the deviation profile). Per set (≥ 10 usable members), the mean member CoM; null: 2,000 random sets of the same
size matched on expression decile (89's `strata_permutations`). z and two-sided p per mouse; p_set = the larger p of
the two mice, same direction required; BH over sets. Gene-sampling nulls are anti-conservative for correlated
members (review, PSEA, CAMERA), so the calls are descriptive.

**Timeline.** For the top 24 positioned redundancy-group representatives (89's Jaccard 0.5 groups) by mean |z|: the
density of member change positions (A's midpoint for S; peak for B; the half-crossing of the fitted curve for L, D
and DS), ordered by median position.

**Sub-trends and a signed summary.** Per set:
- A_med: the range of the median member profile (89's program profile), mean of the mice.
- A_al (cross-validated aligned amplitude): member signs from the first principal component of the standardised member
  profiles in one mouse; in the other mouse, the range of the median of sign-aligned member profiles; the mean of the
  two directions.
- Null: A_al of 500 matched random sets.
- **Hidden by the median:** A_al above its null 95th percentile and ≥ 0.25, while A_med < 0.25; each sign holds ≥ 3
  members. The two sign groups are the sub-trends.

**Fold protection** (sensitivity): positions and A_al on fold-protected members only.

**Figure 90b.** The timeline; set z in mouse A against mouse B; A_med against A_al with the hidden sets marked; one
hidden set's median curve and its two sub-trends in both mice.

## C · Positional information: decoding position from a pathway's genes (recommendation 2)

**Decoder.** Train on one mouse, test on the other (both directions). Per gene, the training mouse's binned rate
r_gb = (counts + 0.5) / library in each of the 25 bins. For a test unit with counts c and library L, the NB
log-likelihood of each bin with the gene's 86 dispersion; the posterior over bins with the training mouse's bin
occupancy as prior. The decoded bin is the posterior mode.

**Comparators.**
- Labels: the posterior is the training mouse's p(bin | label).
- Labels + genes: p(bin | label) × likelihood.
- Random gene sets: fold-protected genes, sizes 10, 15, 20, 30, 50, 75, 100, 150, 200, 300; 200 sets per size,
  matched to the expression-decile composition of the fold-protected universe.
- All fold-protected genes.
- On/off code (Dubuis): detection with Bernoulli per-bin detection rates, for the top sets and all genes.

**Information.** I = MI(true bin; decoded bin) in bits, plug-in, minus the mean of 20 shuffles of the decoded bins,
averaged over the two directions. A lower bound on the information the genes carry about the curve position. Also
the information within labels: Σ_l p(l)·MI(true bin; decoded bin | label). The median absolute error of the posterior
mean (segment units) and σ_x(x), the RMS error per true bin.

**Gene sets.** 71's sets, fold-protected members only (≥ 10). The target is the reconstructed coordinate, not
anatomy; fold protection keeps the genes that built it out.

**Highlight criteria (descriptive).** The number of single pathways whose information exceeds the labels'; the
information within segments that all fold-protected genes carry and the labels cannot; sets above their size's
random 95th percentile.

**Figure 90c.** Information per set against its size, with the random band and the labels' and all genes' lines; the
total split into between-label and within-label parts (labels, top sets, all genes); σ_x(x) along the PT; decoded
against true position in the test mouse for the top set.

## D · Human against mouse without levels, against within-species baselines (recommendation 4, species)

**Genes.** 71's 1:1 universe, usable in both species.

**Summary.** Per specimen and gene, the end contrast E = mean deviation over the last 8 bins − over the first 8
(mouse units for the mice; each section's 25 equal bins for human). A within-specimen log ratio, so a probe-panel
offset cancels. Also the Spearman of the profile with bin order.

**Agreement.** Pearson r of E across genes: Ctrl1A2 vs Ctrl1A4; HUK1_COR1 vs HUK1_MED1; the four human × mouse
pairs; the species means. Conservation corrected for replicate noise: r(means) / √(R_h·R_m), with R = 2r/(1 + r) the
reliability of a two-specimen mean. The human "replicates" are two sections of one donor, so donor variation counts as
species difference and the corrected value is a lower bound.

**Categories.** Graded in a species: |E| ≥ 0.5 in both its specimens with the same sign. Conserved: graded in both,
same sign. Reversed: graded in both, opposite signs. Mouse-only / human-only: graded in one species, |E| < 0.2 in
both specimens of the other. Relabeling control: the same counts for the pseudo-species {Ctrl1A2, HUK1_COR1} vs
{Ctrl1A4, HUK1_MED1} and {Ctrl1A2, HUK1_MED1} vs {Ctrl1A4, HUK1_COR1}.

**Highlight criteria (descriptive).** The corrected conservation; reversed and species-only genes beyond both
relabelings; named anchor genes.

**Figure 90d.** E in mouse A against B; human section against section; human against mouse with the categories and
anchor genes; the agreement and category counts against the relabelings; profiles of conserved and reversed genes in
all four specimens.

## E · Where along the PT the transcriptome changes (recommendation 5)

**Local change.** Per usable gene, mouse and interior position p_j (j = 2..23, p = 0.2..2.3):
Δ_j = mean(d over bins j, j+1) − mean(d over bins j−2, j−1). A gene changes at p_j when |Δ_j| ≥ 0.25 in both mice with
the same sign.

**Comparators.**
- Chance agreement: each gene's mouse-A Δ paired with a random usable gene of the same expression decile in mouse B
  (100 permutations).
- **The discrete transcriptome:** per gene and mouse, the profile a constant-within-label gene would give (each bin
  Σ_l f_bl·R_l, R_l the label pseudobulk rate), run through the same profile, reference and Δ steps.

**Readouts.** Agreed changes per position (rising, falling) minus chance; the same for the discrete transcriptome;
the share of excess change more than 0.3 units from both boundaries ("inside segments"), observed against discrete.

**Pathways by window.** Genes rising or falling at each position, grouped in windows of 0.3 units; one-sided Fisher
of 71's sets (≥ 10 usable members) against the usable genes; BH within window. Descriptive.

**Figure 90e.** The change-density track (observed, chance, discrete) with the label transitions; the pathway ×
window heatmap; the inside-segment shares.

## F · Where along the PT injury acts (recommendation 4, IR)

**Data.** 89's E4 inputs: the AKI export, U_PT units of the four mice, AKI-control mouse units, PT-expressed genes
(≥ 2% in each mouse), 0.1-unit bins. Δ(x) = IR mean − control mean of the binned log2 rates, centred by the median
usable gene (89's `a_delta`); the four IR × control pairs; the two relabelings.

**Responsive.** A gene responds at a bin when |Δ| ≥ 0.5, all four pairs share the sign, and both relabelings are
≤ 0.5·|Δ| there; responsive when this holds in ≥ 2 adjacent bins. Per responsive gene: peak position, extent (bins
responding), effect gradient (Spearman of |Δ| with position over valid bins).

**Position-matched against segment-level fold changes.** Per label, the IR − control log2 FC of label pseudobulks
(what a segment analysis reports) and the position-matched FC (the control bins reweighted to the IR units' bin
distribution within the label). Readouts: the within-label position shift (mean units, IR − control) per label; the
genes whose segment-level |FC| ≥ 0.5 while the matched |FC| < 0.25, or the reverse.

**Highlight criteria (descriptive).** Where injury acts (peak positions, graded vs local effects); genes whose
segment-level change is a within-segment composition shift.

**Figure 90f.** Unit density along the PT per mouse; the effect-location heatmap of responsive genes sorted by peak;
the median |Δ| of up and down genes along the PT with the relabeling band; segment-level against position-matched FC.

## Reading rules

- Two control mice, two IR mice, one donor in two sections. Everything is descriptive unless stated; no p-value here
  is confirmatory.
- C14: statements inside mouse S1 and S2 are descriptive, and position differences below about 0.3 units are below the
  curve's resolution.
- 85's resolution: a step at a label boundary appears about 0.3 units wide on the curve; A's D model includes that
  blur by construction.
- Figures carry plots only (no captions, step boxes or footnotes inside).
- Anything added after the run is labelled post hoc.
