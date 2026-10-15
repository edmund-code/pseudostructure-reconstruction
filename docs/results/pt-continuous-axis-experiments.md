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

## Notebook 90 results

`analysis/notebooks/90_pt_continuous_axis_experiments.ipynb` (logic `90.continuous_axis.1`) writes to
`results/pt_continuous_axis_experiments/`. It followed the protocol frozen in 3f7d89d. Guards G1, G2, G3, G5 and G7
passed. Wall time about 3.5 min cold and 1 min cached; a cold run and a cached run give identical readouts. Module
`pseudospace/continuous_axis.py` with `tests/test_continuous_axis.py` (12 tests).

**How the protocol was read (also in `run_record.json`).** Two readings were decided after the first run, and two
post hoc analyses were added. Everything else is as frozen.
1. **B, p values.** The p per mouse comes from the permutation z (normal), not the empirical p. With 2,000
   permutations the empirical p cannot fall below 5 × 10⁻⁴, and BH over about 1,400 sets then cannot pass. The first
   run gave 18 calls, and 0 on fold-protected members. The empirical p is kept as a column.
2. **B, timeline.** Gradient members are left out, because the half-crossing of a straight line is always the middle
   of the axis.
3. **E, post hoc.** A reproducible change energy, Σ_g Δ_A·Δ_B. Chance agreement makes up 73–95% of the frozen
   agreed-change counts at every position, and their power depends on how many units a bin holds. The product of the
   two mice's changes has the true change as its expectation whatever the noise, so it needs no threshold.
4. **F, post hoc.** The responsive rule applied to the two relabelings as contrasts: the specificity check the
   protocol did not state.

### A · Gene shape classes against segment steps

| Class | All usable | Fold-protected | Same class in both single-mouse fits |
|---|---|---|---|
| discrete (steps at the labels) | 1,171 | 1,106 | 208 |
| discrete + slope within segments | 544 | 462 | 209 |
| gradient | 2,162 | 2,045 | 623 |
| switch | 824 | 748 | 171 |
| bump | 376 | 351 | 63 |
| trough | 633 | 547 | 88 |
| flat | 3,777 | 3,638 | – |

- **Continuous against discrete (fold-protected, non-flat genes):**

  | Reading | Data | Step truth (simulated) | Continuous truth (simulated) |
  |---|---|---|---|
  | continuous favoured (ΔBIC > 6) | 0.281 | 0.028 | 0.459 |
  | steps favoured (ΔBIC < −6) | 0.027 | 0.219 | 0.003 |
  | undecided | 0.692 | 0.753 | 0.538 |

  Shares of the non-flat genes in each column. The data call continuous 10 times as often as a step world would, and
  steps 8 times less often. Post hoc, as a two-component mixture of the simulated readings, the continuous share is
  0.59 (from the continuous calls) to 0.89 (from the step calls). This rests on the simulations matching the data's
  noise.
- **Replication.** 0.283 of non-flat genes have the same class in both mice alone, against 0.172 by chance. Single-mouse
  fits have half the points, so class agreement is a floor.
- **Switches** concentrate in late S1 (0.8–0.9 mouse units), just before the S1/S2 label transition, and at the end of
  the axis. Of the joint switches, 280 sit inside S1, 154 inside S2, 150 inside S3, 169 at S1/S2 and 71 at S2/S3.
- **Bumps and troughs:** 1,009 in the joint fit, 141 replicated (same class in both mice, peaks within 0.3 units;
  115 fold-protected). They sit inside S2 (57), at S1/S2 (45) and inside S1 (34). Among expressed genes (log2 CP10k
  > 0), a mid-S2 peak (1.35–1.55 units) of lipid and peroxisomal enzymes: *Degs2, Hsd17b11, Abhd3, Aldh3a2, Acox3,
  Dhrs4, Gmpr, Kyat1, Fhl1*. A segment analysis calls these "S2"; the curve shows a transient peak that falls before
  S3. 89 found single-gene boundary peaks no more often than in shuffled bins, so the 45 "at S1/S2" are not read as
  boundary programs.
- **References.** *Cyp7b1* switches at S2/S3 (graded) in both mice, matching Ransick's sharp S2/S3 boundary at the
  curve's resolution. *Cyp2e1* is discrete plus within-segment slope, with steps favoured (ΔBIC −26). *Slc34a1*
  and *Slc5a2* have continuous shapes favoured (ΔBIC 86, 42). *Cyp7b1, Cyp2e1, Slc5a2* and *Slc22a7* built the labels,
  so they illustrate the classifier and are not evidence.
- Figure 90a.

### B · Where pathways sit

- **Positioned sets** (mean member centre of mass against expression-matched random sets, both mice, BH q ≤ 0.05):
  30 of 1,404 (15 early, 15 late); 15 of 1,356 on fold-protected members.
  - Early (S1-high): basigin interactions, mineral absorption, inorganic cation/anion and amino-acid transport,
    hypoxia, oxidative phosphorylation, amino acids regulating mTORC1.
  - Late (S3-high): peroxisome, peroxisomal lipid metabolism, PPAR signalling, bile acid metabolism, steroid
    metabolism, vitamins and cofactors, pantothenate and CoA, cysteine and methionine.
  - These match the known S1 transport and oxidative profile and the S3 peroxisomal and organic-anion profile. Set z
    replicates between the mice (figure 90b, panel b).
- **Sub-trends a median curve flattens: none beyond chance.** 19 sets meet the rule, against about 17 expected (5% of
  the 343 sets with a flat median and members of both signs). Random gene sets also hold rising and falling genes
  whose median cancels. A_al exceeds its null 95th percentile in 12.4% of all sets, almost all with a large median
  amplitude too.
- Figure 90b.

### C · Positional information

| Decoder (fold-protected genes, cross-mouse) | Bits | Within segments | Median error (units) |
|---|---|---|---|
| S1/S2/S3 labels | 1.123 | 0 | 0.177 |
| all fold-protected genes (8,897) | 1.505 | 0.541 | 0.127 |
| labels + all genes | 1.512 | – | – |
| Fatty acid metabolism (Hallmark, 123 genes) | 1.367 | 0.458 | 0.120 |
| Glutathione metabolism (KEGG, 29 genes) | 1.366 | 0.468 | 0.123 |
| Xenobiotic metabolism (Hallmark, 129 genes) | 1.360 | 0.469 | 0.121 |

- **26 single pathways localize a unit better than its segment label.** Glutathione metabolism does so with 29 genes
  that did not build the coordinate. 132 of 1,356 sets exceed the 95th percentile of random genes of their size.
- **The labels are worst at the boundaries.** Their position error peaks at 0.45 units at S1/S2 and S2/S3. The gene
  decoders stay at 0.2–0.3 units there.
- **Counts against an on/off code (Dubuis):** detection alone gives 0.25–0.51 bits for the top sets and 0.46 for all
  genes, against 1.36–1.51 from counts.
- The information is about the reconstructed coordinate, not about anatomy; fold protection keeps the genes that
  built it out of the decoders.
- Figure 90c.

### D · Human against mouse without levels

| r of the end contrast (7,971 genes; 7,516 fold-protected) | All | Fold-protected |
|---|---|---|
| mouse Ctrl1A2 vs Ctrl1A4 | 0.891 | 0.881 |
| human HUK1_COR1 vs HUK1_MED1 (one donor) | 0.633 | 0.576 |
| human × mouse pairs | 0.189–0.205 | 0.137–0.154 |
| species means | 0.225 | 0.168 |
| corrected for replicate noise | 0.264 | 0.203 |

- **About a quarter of each species' reproducible PT gradient is shared.** This is a lower bound: donor variation
  counts as species difference.
- **Categories** (graded = |E| ≥ 0.5 in both specimens of a species, same sign):

  | | Species | A2+COR1 \| A4+MED1 | A2+MED1 \| A4+COR1 |
  |---|---|---|---|
  | conserved | 189 | 189 | 189 |
  | reversed | 74 (65 fold-protected) | 0 | 0 |
  | mouse only | 229 | 6 | 8 |
  | human only | 70 | 5 | 8 |

  "Conserved" is the same under every split by construction (all four specimens graded with one sign).
- **Conserved:** *Slc7a13, Slc6a18, Slc7a8, Slc5a2, Slc39a8, Prodh2, Slc22a7, Slc22a8, Psat1, Car2* among others.
- **Reversed** (strongest first): *Apoc3, Gcat, Sh3gl2, Pdk4, Vnn1, Hpgd, Socs2, Oat, Lpar6, Ddit4, Abhd14a, Slc17a3,
  Mme, Prom1, Lap3, Ugt3a1, Kcnk1, Notch1*. For example, *Apoc3* rises toward S3 in mouse (E ≈ +2.9) and falls in human
  (E ≈ −3.1).
- **Mouse only:** *Gldc, Alpl, Cyp51, Spp1, Ehhadh, Acsf2*. **Human only:** *Sel1l3, Serping1, Dpp4, Htra1, Pck1*.
- Readings stay descriptive: one donor, and species registration is unstable (72's C8), so the end contrasts compare
  the first and last 30% of each specimen's own curve.
- Figure 90d.

### E · Where along the PT the transcriptome changes

- **Frozen counts:** agreed changes beyond chance peak at 0.8 (435 genes) and 0.9 (376) and at 2.1–2.2 (337–365).
  At the label transitions themselves (1.0 and 2.0) they are 124 and 277. A step transcriptome would put its changes
  at 1.0 (1,110) and 2.0 (1,081). 28.7% of the excess change sits more than 0.3 units from both boundaries, against
  1.2% for steps; 45.5% of positions are that far from a boundary.
- **Post hoc, reproducible change energy:**
  - 28.3% inside segments (fold-protected 26.3%), against 1.3% for steps.
  - A step transcriptome accounts for 76% of the observed change energy.
  - The energy peaks in late S1 (0.8–0.9) and early S3 (2.1–2.2), flanking the label transitions where steps would
    put it (1.0, 2.0).
  - The late-S1 peak sits next to the curve's densest bin (0.7–0.8, 1.8–2 times the units of its neighbours), so its
    exact place is not resolved. The early-S3 peak is in sparse bins.
- **Pathways:** windowed enrichment (frozen) found 4 sets: EMT falling at 0.8–1.1; SLC transport and inorganic ion
  transport falling at 1.7–2.0; non-homologous end joining rising at 0.5–0.8. The post hoc set energies (figure 90e,
  panel d) show three change zones:
  - late S1: amino-acid, ion, Na⁺/Cl⁻ and hexose transport, scavenging;
  - late S2, 1.6–1.8: SLC transport, lipoprotein assembly and remodelling, LDL clearance, fat and vitamin uptake,
    cholesterol metabolism;
  - early S3: drug and xenobiotic metabolism, CYP450, bile salt transport, steroid biosynthesis, folate.
- C14 applies inside S1 and S2: the zones are descriptive.
- Figure 90e.

### F · Where along the PT injury acts

- **Responsive genes:** 4,456 (2,244 down, 2,212 up). The same rule with the relabelings as contrasts gives 335 and
  281, so the rule is specific about 13-fold (post hoc).
- **Location:** peaks fall in S3 for 1,205 down and 1,044 up genes, against 541/624 in S1 and 498/544 in S2. S3 is
  20% of the axis. |Δ| rises with position (Spearman > 0.5) for 30% of responsive genes and falls for 5%.
- **The injury effect is graded, not stepped.** The median Δ of the up and down genes grows steadily from mid-S2 to
  the end of S3 (about 0.55 to 0.85 log2), with no jump at S2/S3. The relabelings stay near 0. This matches IRI's
  outer-stripe S3 injury (76's priors) and shows the effect building through S2.
- **What segment pseudobulks confuse.** Within S2 and S3, IR units sit 0.12–0.13 units later than control units
  (S1 −0.05). In S2 51 genes and in S3 46 genes reach |log2 FC| ≥ 0.5 in the segment pseudobulk but < 0.25 at
  matched position; 9 and 35 genes show the reverse. The two fold changes agree closely overall (r 0.98–0.99).
- **Limitation:** the package drops much of the injured S3 (`pt-package-issues.md`, item 7), so S3 effects are
  measured on the units it kept.
- Figure 90f.

### What to highlight (descriptive)

1. **Most zonated PT genes change continuously.** A continuous shape beats segment steps for 28% of them, 10 times a
   step world's rate; steps win for 3%, 8 times less than in a step world. Steps still capture about three quarters of
   reproducible change energy, so segments are a good first approximation and the curve adds the rest (A, E).
2. **Position is information the labels discard.** Single pathways not used to build the coordinate localize a unit
   better than its S1/S2/S3 label, and genes carry half a bit inside segments that labels cannot (C).
3. **A level-free species comparison** finds 74 genes with reversed PT gradients (none under either relabeling)
   and shows that about a quarter of each species' reproducible gradient is shared (D).
4. **Injury builds continuously from mid-S2 into S3,** and segment pseudobulks mistake within-segment position shifts
   for expression change in about 100 genes (F).
5. **Null results:** pathways whose median hides opposite sub-trends occur at chance (B), and windowed enrichment of
   change is sparse (E).
