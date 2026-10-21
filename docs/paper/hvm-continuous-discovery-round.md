# Continuous-axis human-vs-mouse discovery: what survived, and what it costs the thesis

Written 2026-10-20. Companion to `human-mouse-paper-ready.md`. This page records one round: six
continuous-axis analyses of the human-vs-mouse PT contrast, each attacked three ways, on **our own four
specimens only**. No external dataset was opened anywhere in the round — no microdissection comparison, no
public atlas, no cross-cohort arbitration. Literature was used only to say what a gene does and to judge
novelty, never to adjudicate a measurement. No p-value gates any claim; the gates are artefact filters,
within-species sign replication, the two committed balanced-relabelling refits, a model-free binned
estimator, and a pure-compartment contrast.

**The round's own verdict goes against the paper's second title claim.** Of the six analyses, four
continuous-only claim families failed, and every surviving gene-level finding is one a discrete
S1-versus-S3 contrast makes more cheaply. The honest version of the methods claim is a screening-and-QC
advantage, not a detection advantage. That is §4, and it is the most important section on this page.

## 0 · What I verified myself before relaying any of it

The round's conclusions overturn statements I made earlier the same day, so each load-bearing number was
recomputed independently against the committed inputs rather than taken from an agent.

| check | measured | what it settles |
|---|---|---|
| gene-order alignment | `max&#124;gs[('condition','T_spatial')] − Z['stat&#124;T_spatial']&#124;` = 9.09e-13 | the 11,055 curves and the statistics table are the same objects in the same order |
| COM is an amplitude readout | Spearman(&#124;COM − 0.5&#124;, curve amplitude) = **0.891** mouse, **0.956** human | "where the mass sits" is numerically almost "how much dynamic range the fit has" |
| the registration is shrinkage | affine COM_h = 0.2719 + 0.4567·COM_m, fixed point **0.5005** against a grid midpoint of **0.49999**; reverse slope 1.0003 | a genuine 2.19× compression requires a reverse slope of 2.1895; no compression is fully consistent |
| amplitude manufactures movers | mouse shape rescaled to each gene's observed human amplitude, **zero re-timing injected**: max&#124;resid&#124; **0.162**, 18 of 95 above the 0.06 mover cut, 9 above 0.08 | the quoted "null floor 0.0434" understates the artefact floor by ~3.7× |
| which movers are artefacts | 7 of 29 reproduced at ≥100% same-sign by amplitude alone: Snta1 1.28, Lgmn 1.21, Slc36a1 1.19, Pter 1.14, Proc 1.12, Fbxo21 1.08, Snx8 1.07 (Gatm 0.80) | four of these were the showcase examples for "registration mattered" |
| the within-species null | trend-sign flips: **0 of 85** mouse-vs-mouse and **0 of 93** human-vs-human among genes with a real trend (&#124;ρ&#124; > 0.5 in both), against **15 of 95** cross-species; the single mouse-mouse flip is Txnip at ρ −0.244 / +0.034, a flat gene | the reversals rest on a null with power, and mouse is two different animals |
| the pure-compartment contrast | 14 of 14 called reversals flip sign with 4-of-4 specimen agreement; stable at purity 0.90 / 0.95 / 0.99 | the surviving evidence needs no registration, no spline and no centre of mass |
| the discreteness recount | of the 29 movers, **16 have opposite-sign end contrasts**, 12 more share a sign but differ in magnitude by 1.5–7.0×, and **exactly 1** (Hao2, +0.652 vs +0.619) is genuinely invisible at the ends | my earlier "21 of 29 are invisible to the endpoint contrast" was wrong |
| the attenuation is already ours | `criteria_hvm.csv` **C12, named "A3: no attenuation flag"**, fraction 0.829, percentile 0.943, verdict **not met** | the committed coordinate validation recorded this before the round began |

Inputs: `results/pt_pathways_package_pseudospace/stage_cache/gene_model_human_vs_mouse_species__06ba95871e922d1a.npz`
(11,055 × 101 fitted log-rate curves, grid 0.0503–0.9497, knots 0.240 / 0.410 / 0.649, mouse = reference),
`tables/gene_statistics_human_vs_mouse.csv`, `tables/s10_mechanics_human_vs_mouse.npz` (20-bin observed rates,
model-free), `results/pt_continuous_axis_experiments/tables/d_end_contrasts.csv.gz`,
`results/pt_package_coordinate_validation/tables/{criteria_hvm,composition,noise_budget}.csv`, and the two
committed relabelling refits `..._swap_HUK1_COR1_Ctrl1A2__3838c468cc28c97a.npz` and
`..._swap_HUK1_COR1_Ctrl1A4__999b3aa3d4a15fe3.npz`.

## 1 · What survived

### F1 · Five genes run their PT gradient in opposite directions in human and mouse
This is the only positive gene-level biology that survived, and it survived in one specific form: a **sign
flip of a pure-S1-versus-pure-S3 contrast**, not a shift, displacement or re-timing of a position. Three of
the round's six analyses collapse into this one claim.

The surviving statistic is registration-free. For each specimen, take only the axis bins whose composition is
≥ 90% PT-S1 or ≥ 90% PT-S3, and contrast mean log1p(observed rate per 10k) between them. No spline, no
knots, no cross-species registration, no centre of mass. (It does use the axis to define the bins — the
honest description is "compositionally pure compartments of the axis", not "segment labels only".)

| gene | mouse S3−S1 (Ctrl1A2 / Ctrl1A4) | human S3−S1 (COR1 / MED1) | detection m/h | amplitude ratio h/m |
|---|---|---|---|---|
| Ugt3a1 / UGT3A1 | +1.168 / +1.289 | −0.789 / −0.832 | 0.754 / 0.497 | 0.755 |
| Dcxr / DCXR | −0.536 / −0.565 | +1.064 / +1.194 | 0.868 / 0.778 | 1.743 |
| Slc17a3 / SLC17A3 | +0.605 / +0.714 | −0.886 / −0.907 | 0.931 / 0.700 | 1.310 |
| Abhd14a / ABHD14A | −0.198 / −0.429 | +1.220 / +1.294 | 0.415 / 0.261 | 2.656 |
| Abhd17b / ABHD17B | +0.208 / +0.214 | −0.341 / −0.371 | 0.354 / 0.352 | 1.178 |

All five carry balanced probes (3 mouse, 3 human), no structural zero, negative edge enrichment (edge-*depleted*)
in all four specimens, and are neighbour-marker-negative for both symbols. Four of the five have human
amplitude **larger** than mouse, so the attenuation artefact that killed most of this round cannot produce
them; for Ugt3a1 the amplitude-only surrogate runs the **wrong sign** (+0.029 against a real −0.155), as it
does for Dcxr, Abhd14a and Slc17a3. All five also reverse in the committed notebook-90 per-specimen end
contrasts, a two-bin statistic computed in a different notebook: Ugt3a1 +2.30/+2.47 → −1.07/−0.85,
Dcxr −0.70/−0.78 → +1.32/+1.38, Abhd14a −1.02/−1.42 → +2.57/+2.11. Leave-one-specimen-out keeps 13 of 14
calls and all five of these.

Ranked by exposure: **Ugt3a1** is the one to stake an experiment on (survived every attack at full strength).
**Dcxr** is the cleanest curve pair in the cohort (mouse monotonicity 0.985, interior-knot ρ −1.000) and the
only gene whose reversal also holds *within* pure S1 (−0.92 vs +0.93) and *within* pure S3 (−0.99 vs +0.84).
**Slc17a3** has the cleanest neighbour profile in the whole core but holds as a direction reversal only — its
mouse within-S3 trend is flat (+0.02) and it fails the shape half of its own call model-free.
**Abhd14a** is the largest effect anywhere and the most exposed (detection_human 0.261 against a 0.25 floor,
24–30% of its human mass in the terminal grid band): report its direction, never its magnitude.
**Abhd17b** is the most novel and the thinnest — 47 and 56 raw human S3 counts, reversal confined to the
distal half, and the only survivor outside the union of all 83 called pathways, so our own pathway screen is
structurally blind to it.

### F2 · The human arm has systematically shallower gradients than the mouse arm
The largest effect in the dataset, upstream of every per-gene claim, and the source of most of the round's
artefacts. Median human/mouse curve log-amplitude **0.633** on the fitted curves and 0.73–0.79 on the raw
20-bin observed rates; human flatter for **75.8%** of the 95 core genes. 99% of 478 core gene pairs with a
mouse interquartile gap > 0.05 close in human, median closure +0.322. Our own committed validation already
recorded it: C12 "A3: no attenuation flag", fraction 0.829, **not met**.

**Its cause is not resolvable in this cohort, and that is the single most important open item.** Depth,
probe-panel composition, segmentation quality, axis quality and one-donor sampling are perfectly confounded
with species in a 2-vs-2 design, and all of them are present identically in both arms of both relabelling
refits. Internal columns do not pin it down: Spearman with the per-gene human/mouse probe ratio is −0.024 and
with detection_human only −0.19. The asymmetry is present in the raw binned rates, which excludes the spline
basis but not the coordinate.

Worst-attenuated genes, which is why they generated the worst artefacts: Gatm (amplitude 5.267 → 0.751, ratio
0.143), Slc34a1 (0.100), Slc22a6 (0.153), Nt5e (0.187), Miox (0.235), Pm20d1 (0.362). Opposite direction:
Abhd14a (2.656), Pter (2.537), Slc36a1 (2.207), Gpt (2.033).

### F3 · Methodological: the softmax centre of mass is not a position statistic on this cohort
`COM = softmax(log-rate) @ grid` is exactly invariant to an additive per-gene per-species log offset. That is
true, and it is what "probe-panel safe by algebra" rested on. It is also **maximally sensitive to curve
amplitude**, and amplitude is precisely what differs by species here (F2). A flatter curve's mass centre is
dragged toward the grid midpoint, so a contrast-strength difference reads out as a positional difference.

Measured: Spearman(|COM − 0.5|, amplitude) 0.891 mouse / 0.956 human; the affine registration's fixed point
is the grid midpoint to four decimals; the identity b = r · sd_h/sd_m = 0.6759 × 0.6757 = 0.4567 holds
exactly, i.e. the slope is regression dilution times a spread ratio. Two amplitude-free location statistics
built with measured within-species reliability collapse the cross-species map: x50 (within-mouse r 0.738,
within-human 0.660) gives cross-species r of +0.015 / +0.081 / −0.028 / −0.003 across the four pairings, and
nCOM gives +0.246, against softmax COM's +0.676 to +0.683. The sharpest version: among the 39 core genes
strongly monotone in both species, softmax COM correlates across species at +0.800 while x50 gives +0.219 —
**the registration looks best exactly where the statistic least measures position.**

The committed relabelling refits cannot see this failure mode, because each arm mixes one mouse with one
human specimen, so the species-level amplitude asymmetry cancels in both arms by construction. Their zeros
(0 of 95 above threshold, slopes 0.906 and 1.061) are closer to identities than to measurements: the arms'
per-gene COMs correlate at Spearman 0.995 and the maximum attainable arm-to-arm pair gap is 0.081 against the
0.160 a sign flip would need. **The null that carries the reversals is the within-species one** (§0), which
this project had not been running.

### F4 · Negative result: positional movement and co-expression rewiring are independent here
"Movers are ~2.3× more rewired" is arithmetically correct against the wrong comparator. Movers' median
residual-layer rewiring D is 0.3036 against a table-wide 0.1333 = 2.28× — but the 54 **non-movers** that
passed the same core gate sit at 0.3073 = 2.30×, so the mover-minus-non-mover difference is −0.004. The
mechanism is a selection effect that never looks at movement: D tracks neighbourhood reproducibility at
Spearman +0.738 (W_m) and +0.643 (W_h), and the movement-blind filter W_m ≥ 0.7 ∧ W_h ≥ 0.7 alone selects 551
genes at median D 0.3090 = 2.32×, reproducing the whole headline with no reference to movement. The specified
2×2 is degenerate — both bottom-quartile cells are empty. The layer also fails its own committed gate (Q1
reproducibility 0.670 against 0.7).

Two consequences for wording. No gene list should be carried from this layer: the three within-core quadrant
lists are a reproducibility split, with the "moved but not rewired" genes (Slc36a1, Txnip, Lgmn) simply the
three lowest-W members of the group, and Slc5a2 topping the ranking at W_m exactly 1.000 with a 7× layer
disagreement — **do not let Slc5a2 carry an SGLT2 rewiring reading**. And 22 of 95 core genes are absent from
the rewiring table altogether, 10 of them movers including the two largest, so any rewiring statement is
silent on a third of the core and on the strongest signals.

## 2 · Hypotheses for the lab

One row per testable prediction. Every one is scored **anatomically** — glomerulus-attached convoluted
profiles as S1/S2, medullary-ray and outer-stripe straight profiles as S3 — with no marker-gene anchor.
That is deliberate: `Slc34a1`, the panel anchor the round first proposed, has mouse amplitude 3.41 and human
0.31 and marks no position in the arm the panel exists to read. The only core gene graded in both species
that could serve as an anchor is **Slc4a4** (detection 0.866/0.921, amplitude 2.62/1.08, within-species ρ
0.97/0.99, 3/3 probes).

| # | prediction | gene(s) | assay | falsified if | confidence |
|---|---|---|---|---|---|
| H1 | UGT3A1 conjugation capacity is proximal-weighted in human PT, distal-weighted in mouse | Ugt3a1 / UGT3A1 | RNAscope, one mouse kidney + one human cortex/outer-medulla block | human straight ≥ convoluted, or mouse convoluted ≥ straight | strong |
| H2 | DCXR (polyol/uronate cycle, dicarbonyl detox) is distal-weighted in human, proximal-weighted in mouse | Dcxr / DCXR | IHC (abundant protein, antibodies available), same blocks | human straight not the high compartment, mouse convoluted not — or signal proves non-epithelial | strong |
| H3 | **Either** human PT zonation is genuinely shallower than mouse (biology, and then this is the headline) **or** our human arm under-reports amplitude (measurement, and then every magnitude in the project is recoverable) | Nt5e, Miox, Slc34a3 + Lrp2 as positive control | RNAscope, human block, S1-vs-S3 fold change per gene against the mouse values | human fold changes comparable to mouse across unrelated genes ⇒ artefact; ~0.6× across unrelated genes ⇒ real species finding | strong |
| H4 | Our human arm cannot recover a known gradient: megalin/cubilin have essentially no axial gradient in our human sections (amplitudes 0.170 / 0.197, the lowest checked; LRP2's two human sections anti-correlate at ρ −0.168) | Lrp2 / LRP2, Cubn / CUBN | megalin IHC, same human block | a textbook S1-weighted megalin gradient on the slide settles H3 as measurement in one stain | strong |
| H5 | SLC17A3/NPT4 apical urate-and-anion efflux is proximal-weighted in human, distal-weighted in mouse — predicts where luminal urate peaks and where secreted anionic drugs concentrate | Slc17a3 / SLC17A3 | RNAscope, same panel | human straight ≥ convoluted or mouse convoluted ≥ straight | moderate |
| H6 | ABHD14A (orphan α/β-hydrolase, no renal literature) rises toward distal PT in human, flat-to-falling in mouse | Abhd14a / ABHD14A | RNAscope with a PT mask | human distal ≤ proximal, or the human signal is not in PT epithelium | moderate |
| H7 | ABHD17B (protein depalmitoylase family) rises toward S3 in mouse and falls in human, confined to the distal half | Abhd17b / ABHD17B | RNAscope, scored as a distal-half contrast | human straight is not the low compartment | moderate |
| H8 | The steep distal-rising CD73 gradient in mouse PT is **absent** in human — a loss of gradient, not a reversal; the human limb must be reported as "not measured" | Nt5e / NT5E | CD73 IHC with **mouse as the positive control** | human shows a clear monotone S1→S3 gradient of mouse-like magnitude | moderate |
| H9 | Human PT has a distal-rising PCK1 gradient; mouse has no usable longitudinal gradient. **This is a validation of the method, not a discovery** | Pck1 / PCK1 | PEPCK-C IHC, cortical labyrinth vs medullary ray and outer stripe | human medullary-ray staining ≤ cortical labyrinth | moderate |

Notes that belong with the panel. H5 cannot address the dominant human urate pathway at all: ABCG2
(detection_human 0.062, human edge enrichment 1.28–1.40, human neighbour marker), SLC2A9 (detection_mouse
0.080, mouse neighbour marker) and SLC22A12 (moves +0.005–0.013, below its own null, and a
coordinate-construction gene) are all excluded — say so rather than implying coverage. H9's human limb is
strong (detection 0.998/0.974, pure-compartment +0.787/+0.903) and its mouse limb is weak and must be stated
as weak (pure-compartment only −0.313, mouse argmax at interior grid 0.401, trend −0.497 when either mouse is
dropped, committed category literally "human only").

**Literature status, flagged against the standing full-text rule.** Every PMID in this round was read at
title/abstract/metadata level only; **none was read in full text in this pass.** Two of those reads are
load-bearing and must be re-read in full before anything is written: **PMID 11716765** (isolated human PT
segments — S2 and S3 synthesise more glucose from lactate than S1, which is why H9 is a validation and not a
discovery) and **PMID 2535703** (rat CD73 activity highest in S1 — the *opposite* orientation to our mouse
mRNA gradient, which is why H8 must be stained before any CD73 claim is published). Literature was used only
for gene function and novelty classification; no literature value was used to validate, calibrate or veto any
measurement, and where it disagrees with us it became an experiment rather than a verdict.

**Novelty, as classified by the same pass (abstract-level).** Discovery-grade: the UGT3A1 intra-renal axial
reversal (the organ-level tissue-distribution species difference is already published, PMID 18981171, but the
axial reversal is not); the SLC17A3/NPT4 axial direction difference (its apical urate-efflux function and gout
genetics are established, PMID 20810651); DCXR's segmental distribution in either species. Most novel and
least interpretable: ABHD17B and ABHD14A, with no renal literature at all in either species — which is why no
physiological reading should be attached to them yet. Novel but weakened genes, carried for completeness:
REEP6 (no renal or segmental data), LAP3 (its mercapturate / Cys-Gly-conjugate dipeptidase role is
established, so a reversal would be a species difference in where glutathione-conjugate degradation happens),
PM20D1 (enzymology established, PMID 27374330; no renal axial data). Not a discovery: PCK1 (H9).

## 3 · What did not survive

All of it shares one cause: the centre of mass and its quantiles read amplitude as position, and the human
arm is flatter.

- **The compression framing.** "The human PT runs the same genes on a positionally compressed axis, b = 0.457"
  — the regression's fixed point is the grid midpoint, the reverse slope is 1.0003 where a real 2.19×
  compression demands 2.1895, so the errors-in-variables bound on the true slope is [0.457, 0.9997] and no
  compression at all is consistent. Amplitude-matching moves b to 0.73.
- **The null floor.** "0.0434, so effects are 1.9–5.8× the floor" — the amplitude-only artefact floor is
  0.151–0.162 with 9 of 95 genes above 0.08.
- **"15 earlier / 14 later, so the set is not one-sided"** — void twice: OLS residuals are mean-zero by
  construction, and the amplitude artefact is itself balanced 14/15 on the same genes. Delete the sentence.
- **The 19 headline magnitudes.** On the 0.20–0.80 window, the only one our own filter validated, the headline
  falls 19 → 7 and the movers 29 → 14, while signs hold 29/29 at Spearman 0.981. Nine movers carry ≥ 15–20%
  of their mass in the terminal grid points.
- **Nine of the 29 movers are amplitude artefacts**, with Gatm 147%, Snta1 109–128%, Slc36a1 104–122%, Pter
  106–123%, Lgmn 100–121%, Proc 112%, Fbxo21 108%, Snx8 107%; for five of them the amplitude-invariant peak
  moves the **opposite** way to the reported direction. Lgmn, Slc36a1, Snx8 and Gatm were the showcase
  examples for "registration mattered" — the section meant to prove the correction necessary is where the
  correction manufactured its own movers.
- **Nine of the fourteen reversal genes.** Dropped with the reason: Pck1 (mouse mid-axis hump, the shape for
  which Msrb1 was correctly excluded; fails leave-one-out) → demoted to validation H9; Miox and Cox7a1 (both
  `fold-protected = False`, i.e. coordinate-construction genes — only 5.7% of genes are unprotected, and a
  gene used to fit an axis cannot be independent evidence about its own position on it; Miox's whole-axis flip
  is a peak shift, with no within-segment reversal at all); Nt5e (human 54/47 counts, monotonicity 0.156) →
  H8 as a loss of gradient; Pm20d1 (~60 human counts); Rhoq (only core gene whose top non-PT absolute
  abundance exceeds the neighbour filter's 1.0 CP10k floor in **both** species; the only gene sign-discordant
  inside both relabelling refits); Lap3 (94/118 human counts); Slc16a9 (human amplitude 0.504, fails the
  segment-anchored test); Reep6 (human section slopes −0.91 vs −1.62, worst human reproducibility of the 14).
  The word **"monotone"** must go from the reversal wording: only 1 of the 14 has both limbs at monotonicity
  ≥ 0.80 and it is Rhoq.
- **The coupled-pair analysis entirely.** Pck1–Fbp1, its cleanest result: 88% of Pck1's move is predicted by
  pure segment re-allocation, and Fbp1 has no axial gradient in either species (amplitudes 0.404 / 0.347, COM
  pinned 0.022–0.032 from the grid midpoint), so "in mouse they sit on top of each other" is two
  uninformative numbers that happen to be adjacent. The Slc22a7/OAT2 collapse reverses: four of five gap
  ratios sit **at** the cohort-wide contraction slope (0.516) and OAT2's own residual (+0.041) is below its
  own null (0.059) — OAT2 moves *less* than the cohort trend. The method sentence "a gap cancels the global
  compression, so it needs no registration" is false for a rescaling, which is what is present.
- **The order-inversion count.** "151 of 1,932 pairs" moves to 116–131 under landmark re-registration, 125
  depth-renormalised, 145 composition-renormalised, 148 over random monotone warps, 204 under x50, with the
  denominator moving 1,780–2,121. The unit is ~13 genes, not 65: 25 of the 65 "inverting" genes have
  |COM shift| under 0.03, inside both the 0.012–0.016 reliability and the ~0.08 normalisation offset. The
  claimed monotone-invariance is false as mathematics (a mass centre's order is not monotone-invariant; a
  quantile's is), though it happens to hold empirically here. The 12+1 multiplex panel should not be ordered.
- **The onset/offset asymmetry entirely.** The direction of all 15 calls is predicted **15/15** by whether
  that gene's human curve is flatter or sharper than its mouse curve, with no positional input; the amplitude
  surrogate reproduces 51–123% of 14 of the 15 and collapses Spearman(d25, d75) from 0.821 to 0.468;
  re-calling on the residual leaves 2 of 15. Five of the eight "onset" genes sit at the amplitude filter
  floor with human q25 spanning 0.324–0.352 — five nearly identical near-flat curves, not five advanced
  onsets. The OAT2/OAT3 handoff **reverses**: amplitude-free the gap widens by 0.046, in all four pairings
  individually, against a reported closure of 0.264 — and it is unremarkable, sitting at the 6th–9th
  percentile of position-matched pairs.
- **Every sentence naming a human segment.** Human PT-S2 has segment length 0.0609 in HUK1_COR1 and
  **−0.0769** in HUK1_MED1, where the S2→S3 transition (0.5895) occurs *before* S1→S2 (0.6665); human S1
  occupies 0.543–0.666 of the axis against mouse 0.342/0.339. Grid 0.40 is mouse S2 but human S1, and human
  S2 is not an interval. "Peaks in human S2" is not a statement this cohort can make.
- **Every effect size for every surviving gene.** Human fold-refit noise in axis units is 0.163–0.221, so 27
  of the 29 reported shift magnitudes were smaller than the smallest human noise floor. Direction only.
- **"COM is probe/depth-safe by algebra"** as a general reassurance. True for an additive offset, false for
  scale, and the operative species asymmetry is scale.
- **The SGLT2→SGLT1 handoff, in either direction** — which matters because collaborators asked for it.
  Slc5a1 has detection_human 0.061 and is a mouse neighbour marker (atlas ratio 7.12); the nearest testable
  substitute Slc5a2–Slc2a2 is flat (+0.021 against a 0.072 threshold). Any SGLT2i/glucose question needs an
  assay that actually detects SGLT1 in human; this panel does not.
- **The CYP4A story** (one-to-many orthology across a family where mouse 20-HETE comes only from
  male-specific Cyp4a12a — the textbook species difference; and the amplitude surrogate reproduces 100% of
  its move) and **PTER** (106–123% amplitude-reproduced, amplitude-invariant boundary stationary at spread
  0.005, human curve 2.5× sharper) — a shame, because PTER's 2024 identification as the physiological
  N-acetyltaurine hydrolase is fresh biology; it is just not this gene's result here.
- **The megalin/cubilin endocytosis spacing story** — not measurable in our human arm, and that failure is
  now hypothesis H4.
- **The one passing rewiring verdict** (Q2, ρ 0.297) is the external arm and must not be cited under the
  own-data-only scope, leaving the rewiring layer with no independent corroboration inside this cohort.

Bookkeeping, for trust rather than consequence: the step-8 count of 163 reproduces only under an unstated
missing-value policy (naive readings give 165/167/169/171 and a 96-gene core; Dpep1 is the single gene whose
membership turns on it, and it is not a reversal either way); `pair_positions.csv` has 61 unique ordered pairs
with 54 having x\* in both species, not 72 and 65; the atlas ratio "9.0/14.3" is Car12's mouse/human pair, not
one value each for Car12 and Gpx3 (Gpx3's mouse ratio is 50.88); Abcg2's 1.397 is the HUK1_MED1 value and
1.276 the filter-relevant minimum — both exceed 1, exclusion unaffected. Msrb1 is correctly **not** a
reversal (|trend_mouse| 0.471 Pearson, 0.394 Spearman, below 0.5 by every definition), but its exclusion is
one specimen and one statistic deep: it becomes a 15th call if Ctrl1A4 rather than Ctrl1A2 is dropped.

## 4 · Continuous versus discrete, answered as a count

**Of the 5 surviving gene-level findings, 5 of 5 are findable by a discrete pure-S1-versus-pure-S3 contrast
and 0 of 5 require the continuous axis** — because the discrete contrast *is* the evidence that let them
survive. All five also reverse in the committed notebook-90 two-bin end contrasts, computed in a different
notebook. F2 (the attenuation) is also discretely findable: an S1-vs-S3 fold change per gene compared across
species gives the amplitude ratio directly. The remaining two survivors are a negative result and a
methodological result *about* continuous statistics.

**The count the thesis rested on measures 1 of 29, not 21 of 29.** Reading the end contrasts themselves
rather than the category labels: 16 of 29 movers have **opposite-sign** end contrasts — the most
endpoint-visible pattern a two-end statistic can have — and are unlabelled only because `graded_categories`
requires |E| ≥ 0.5 in both specimens of both species; 12 more share a sign but differ by 1.5–7.0×; exactly
one (Hao2, +0.652 vs +0.619) is genuinely invisible. "Conserved" in that function means *same sign with both
species graded*, not *identical between species* — Gatm is labelled conserved with E −5.59 mouse against
−0.80 human. And the sharpest leg, "8 movers the endpoint screen calls identical", collapses to 1 on the
surviving set, mechanically: a gene whose two ends match has a small raw dCOM, so its residual is dominated
by flattening. The most rhetorically useful part of the count was the most contaminated part.

**Four of four continuous-only claim families died, and the continuous statistic was worse than neutral.**
Centre-of-mass re-timing, onset/offset asymmetry, coupled-pair co-location and pairwise order inversion all
failed for the same reason, and the statistic *manufactured* claims a discrete contrast cannot make: 9 of 29
movers are genes whose curve flattened or sharpened, and for five the amplitude-invariant peak moves the
opposite way. A per-gene S1-vs-S3 fold change is immune to that failure mode by construction. Two further
counts point the same way: the inversion count moves 116–204 across equally legitimate handlings while its
sign field is stable (instability without detections), and the 15 boundary calls have their direction
predicted 15/15 by the sign of the amplitude change alone.

**What is genuinely continuous-only: one sub-property of one gene.** Dcxr reverses *within* pure S1 (mouse
−0.92 vs human +0.93) and *within* pure S3 (−0.99 vs +0.84) with both curves monotone — a statement a
three-segment analysis cannot make, and the only one in the round. Abhd17b's distal-half confinement and
Slc17a3's flat within-S3 mouse trend are two further sub-segment shape facts, and both are visible to a
three-segment analysis. Note also that a three-segment comparison is impossible in *our* human arm regardless
(S2 degenerate, transitions inverted) — so the discrete comparator that beats us is **S1-vs-S3, not
S1/S2/S3**.

**The verdict for the paper.** On this cohort, continuous did not beat discrete at detection. Every survivor
is a claim a discrete S1-vs-S3 contrast makes more cheaply and more robustly, and several continuous-only
claims were artefacts of the continuous statistic itself. The defensible, narrower claim is a
**screening-and-QC advantage**: the continuous reconstruction is what supplied the per-gene quality machinery
— within-species curve reproducibility, amplitude floors on the full and trimmed grid, monotonicity and
within-segment shape, edge-mass and endpoint-pinning diagnostics, and the reduction from 11,055 genes to a
95-gene testable core — and that machinery is what tells you which discrete calls to trust. It is also what
exposed, through C12 and the amplitude surrogate, that our human arm under-reports gradient amplitude; a
discrete analysis would have reported the same shrunken fold changes without ever flagging them. That is a
real contribution, it is honest, and it is a different paper from the one in `outline.md`.

This does **not** settle the AKI axis, where the discreteness evidence is independent (notebook 90's shape
classes, 26 pathways localising a unit better than its own segment label, change density beside rather than
at the transitions). It settles the species axis, where the continuous claim was always weakest.

## 5 · Limits

Two mice, and a human arm that is two sections of **one donor** — so every human statement rests on a single
kidney and the human half of all four-way replication is within-donor technical replication. The two
committed relabelling refits put one mouse and one human specimen in *each* arm: they control for "this
particular block is odd" and can say nothing about species-level measurement asymmetry, and they are
additionally near-incapable of producing the statistics they were used to bound (arm COMs correlate at 0.995,
maximum attainable pair gap 0.081 against the 0.160 needed). The within-species null is the one that carries
the reversals.

The human axis is demonstrably the weak half: C1 (canonical order), C4 (human gene-fold 0.488), C6 (P5 human
S1 0.552), C7 (human S1 gene-fold 0.516), C8 (registration gap 0.369 against a 0.15 rule), C12 (attenuation
0.829) and C14 (10 of 12 specimen-by-segment cells noise-dominated, median ratio 1.734) are all recorded
**not met** in `criteria_hvm.csv`. Human PT-S2 has length 0.0609 and −0.0769 with inverted transitions in
HUK1_MED1. Human fold-refit noise in axis units is 0.163–0.221, which is why only direction is reportable.
The human pure-S3 compartment carrying the surviving evidence is thin (366 and 449 units against 1,661 and
1,743 in pure human S1), so several survivors rest on 47–120 raw human counts.

These are lab-testable hypotheses from four specimens, not findings.

## 6 · Open questions this page does not decide

1. Whether the paper's methods claim is rewritten as the screening-and-QC advantage (§4) or whether the
   species axis is dropped from the methods argument and carried as application only. `outline.md` and
   `paper-goal-methods` both currently assert "continuous beats discrete" as a title claim.
2. Whether F2's cause-resolving stain (H3/H4) runs before any magnitude is published anywhere in the
   project — it bears on the AKI axis too, since both arms share the human-side machinery.
3. Whether the 0.20–0.80 window becomes the mandatory reporting window for every curve statistic in the
   package, given that 9 of 29 movers carry ≥ 15% of their mass in the terminal grid points.
4. Whether the within-species trend-sign null is added to the package as a standard control, since it has
   power where the balanced relabelling does not for species contrasts.
