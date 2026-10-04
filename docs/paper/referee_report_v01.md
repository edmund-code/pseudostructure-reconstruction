# Referee report: "Where human and mouse proximal tubule differ along the tubule axis"

Manuscript: `docs/paper/draft.md` (assembled 2026-10-04). Reviewed against `docs/paper/outline.md`,
`docs/results/pt-pathway-method-selection.md`, `pt-zonation-amplitude.md`,
`pt-pathway-literature-novelty.md`, `docs/paper/r1_draft.md`, `r3_draft.md`, `r5_draft.md`,
`zonation_literature.md`, the saved tables under `results/` (as of 19:35 EDT on 2026-10-04;
`pt_zonation_classes/` was regenerated at 19:19–19:22 with logic version `.2`, which adds settings
A6 and A7), the notebook 35–40 source scripts in SCRATCH (`ws1/nb35.py`, `ws1/nb36.py`,
`ws1/plan.md`, `ws2/nb39.py`, `ws2/nb40.py`, `ws3/nb38.py`) and `pseudospace/zonation_classes.py`
and `pseudospace/zonation_amplitude.py`.

The issues already being fixed are excluded: the 34-row robust-pathway table, the Figure 3a
counts, the old reversal labels in the novelty doc and `r5_draft.md`, R1's attenuation bound, the
03/13 partition concordance, the missing lab records and figure assembly. Where a new problem sits
next to a known one, I say so.

Where I quote a number that is not in the manuscript, I give its table. A few numbers come from my
own quick checks on saved tables. They are marked **[reviewer check]**, use raw Spearman
correlations, and should be redone properly by the authors.

---

## Summary judgement

The manuscript asks a good question: does human PT zonate the same genes as mouse PT? It is
unusually candid about its own weaknesses. Several post-hoc decisions are disclosed, the coordinate's
limits are stated, and pathway-level external replication is reported as failing. Most numbers in
the abstract and Results match their source tables. I list the exceptions in section C.

As written, I would not accept the central biological claim: "beyond a conserved core of transport
markers, human and mouse PT zonate largely different genes". The direction of the result is
probably right. Independent public atlases give a similar picture. The size of the effect and its
attribution to species, however, rest on four things the paper has not yet shown:

1. The human S1/S2/S3 labels are the human members of cross-species Harmony clusters and have no
   anatomical validation. The species asymmetry sits almost entirely in the S2 − S1 contrast, which
   is the contrast where those labels are weakest.
2. No same-species, cross-dataset benchmark tells the reader what correlation to expect between two
   datasets of the same species. My quick check suggests the mouse–rat correlation is as low as the
   human–mouse one.
3. "Flat" calls, which create the species-only classes, are never confirmed externally. Human
   flatness comes from two sections of one donor.
4. Several species-level confounders that external human data share are not discussed: tissue
   procurement, age, fasting/diet/circadian state and disease.

The methodological claim, that only position-dependent screens are specific with two specimens per
species, rests on two relabelings. The coordinate on which it holds was chosen after the
pre-registered rule had selected a different coordinate, on which it fails. The split-plot argument
does not by itself predict the result. The pathway layer mostly restates mouse zonation, as the
authors partly acknowledge.

The method contribution and the biology are not separated. The headline biology (R3, R5) uses
neither the coordinate (R1) nor the pathway screen (R2). A reader will ask what the coordinate and
the Visium HD data contribute to the central claim, which the authors also reproduce from public
snRNA alone (external r = 0.23).

**Recommendation: major revision.** Most of the decisive checks are feasible with data already in
hand (section D).

---

## A. Major issues

### M1. Human segment labels are not anatomically validated, and the central asymmetry is concentrated where they are weakest

**Location.** R1 Cohort ("each labelled S1, S2 or S3 by marker review"); R3 throughout; Methods
6.7–6.8; Abstract.

**Problem.** Every R3 contrast uses the reviewed labels. These labels are joint Leiden clusters on a
Harmony embedding in which batch = specimen and is fully confounded with species. θ = 6, which is
strong mixing pressure. A human structure is "S2" because Harmony co-clusters it with mouse S2.
Nothing in the human data independently shows that it is S2. If the human S1/S2 boundary does not
correspond anatomically to the mouse one, human S2 − S1 gradients are diluted or rotated. That
lowers the cross-species correlation and turns truly zonated human genes into "flat" or
indeterminate ones. This is a direct alternative explanation for the main result.

**Evidence.**
- Human S1 makes up 55% of human PT structures (1,890/3,438 and 2,076/3,768), against 33% in mouse
  (Methods 6.8 table). The S1→S2 transition is at 0.35–0.38 in human and 0.26–0.28 in mouse.
- None of the PT-S2 reference-panel markers is zonated in human. Slc22a6 and Slc13a3 are flat in
  human (mouse-only class), and Cyp2e1 is not measurable. Of the six PT-S1 panel markers, only Slc5a2
  is zonated in human: Gatm and Slc5a12 are mouse-only, and Lrp2, Cubn and Slc34a1 are flat or
  neither (`gene_zonation_classes.csv`). The human S1/S2 split is therefore not supported by the
  panel used to name it.
- The asymmetry is in S2 − S1: mouse-only 219 against human-only 44. For the region-matched
  S3c − early contrast it is 48 against 35 (`class_sizes_by_contrast.csv`). The noise-corrected
  cross-species correlation is 0.14 for S2 − S1 but **0.29 (0.23–0.34) for S3c − early**
  (`posthoc_reference_free_amplitude.csv`). The paper reports the first and not the second.
- Spatial label coherence is weaker in human. The same-label neighbour fraction is 0.59–0.61 against
  0.41 expected in human, and 0.71–0.72 against 0.33 in mouse
  (`phase1_diagnostics/physical_neighbour_diagnostics.csv`).
- Integration genes are concentrated in the conserved class. Notebook 13's 611 construction genes
  make up 15.5% of conserved genes, against 2.0% of mouse-only and 6.0% of indeterminate genes. Among
  zonated construction genes, 24 of 35 (69%) are conserved, against 131 of 485 (27%) of other
  zonated genes **[reviewer check: `nb35/gene_sets_and_external_contrasts.csv` ×
  `gene_zonation_classes.csv`]**. The labels are aligned on the genes that are then called
  conserved.

**Fix.**
1. Validate the human labels anatomically on the existing polygons:
   - the fraction of PT polygons touching a glomerulus polygon (S1 by definition at the urinary pole)
     that are labelled S1, per species;
   - distance to the nearest glomerulus by label in both species;
   - medullary-ray location of human S3.
2. Re-derive labels by species-separate clustering with marker review, and repeat R3.
3. Repeat R3 after removing the integration genes.
4. Report S3c − early as a co-primary contrast with its correlation (0.29) and SD ratio (0.84).
5. State in R3 that labels are cross-species cluster memberships.

### M2. No same-species cross-dataset ceiling, and the rat reference implies one is needed

**Location.** R3 "Human PT is not simply a flatter mouse PT"; Abstract; Discussion point 2; R1 table
(rat omitted); R5 rat statements.

**Problem.** A noise-corrected correlation of 0.14–0.23 means "different genes" only relative to
what two datasets of the *same* species give for the same contrast. The only comparator offered is
mouse AKI against control at 0.83, which uses the same lab, platform, pipeline and labels, so it is
not a fair benchmark. Disattenuation corrects sampling noise only. It does not correct differences in
label definition, dissociation or platform.

**Evidence.**
- Within-species, cross-dataset agreement from saved tables:
  - mouse snRNA against microdissection: 0.55 (S2 − S1) and 0.70 (S3 − S1);
  - human (Lake S3 − S1 against the Census human convoluted-versus-S3 contrast): only 0.35
  **[reviewer check, raw Spearman over 6,357 genes detected in ≥ 10% of structures in both
  species]**.
  Human references therefore disagree with each other far more than mouse references do.
- Cross-species, from the same external data: Lake against mouse gives 0.11 (snRNA) and 0.03
  (microdissection) for S2 − S1.
- **Mouse against rat microdissection: 0.10–0.13**, with **55% direction agreement** for strongly
  zonated rat genes (`pt_external_segment_validation/coordinate_check.csv`;
  `nb35/p1_external_concordance.csv`, rat_micro 0.105). The R1 table omits rat.

**Fix.**
1. Report disattenuated within-species, cross-dataset correlations for S2 − S1 and S3 − early:
   - our mouse against Census mouse and against microdissection;
   - Census mouse against microdissection;
   - our human against Lake;
   - Lake against Census human.
2. Express the cross-species value as a fraction of the within-species ceiling.
3. Decide about rat:
   - If the Lee 2015 rat data are valid, rodent–rodent zonation is as divergent as human–mouse.
     The framing must then change from "human is not a flattened mouse" to "PT zonation programs are
     evolutionarily labile".
   - If the rat data are not valid (likely, given RPKM, few replicates and case-insensitive symbol
     matching), remove the rat-based arguments. These include "rat runs the opposite way, so this is
     a mouse program, not a rodent program" for sterol synthesis, and the rat support cited for Dcxr
     and Rbp4.

### M3. The species-only classes depend on unconfirmed and asymmetric "flat" calls

**Location.** R3 "Zonation classes" ("Species-only genes ... the other species was confidently
flat"); Methods 6.20; Introduction ("confirmed every call in an independent atlas").

**Problem.** In `pseudospace/zonation_classes.py`, `confirm()` touches only up/down calls. **Flat
calls are never checked against the external atlas.** A mouse-only gene is therefore a confirmed
mouse gradient plus a human interval inside ±0.3 log2 computed from two sections of one donor. Their
spread is section-level, not donor-level, variance. Further problems:
- The flat margin is not scaled to the lower human amplitude. A3 scales only the effect floor.
  Truly conserved genes with mouse gradients of 0.5–0.8 log2 and human ≈ 0.58 × mouse can therefore
  be called mouse-only.
- Confirmation power is asymmetric: 12 male mouse donors against 7 mixed-sex Lake donors, whose S2
  pseudobulks come from only 59–140 nuclei in five of the seven donors (`run_manifest.json`,
  `lake_primary_donors`).
- The Introduction's "confirmed every call" is not accurate.

**Evidence.**
- A3 narrows 247:66 to 243:141. A1 (no confirmation) gives 294:133.
- With the full mouse S3 contrast, where mouse power for flat calls is higher, human-only exceeds
  mouse-only (77 against 111; `class_sizes_by_contrast.csv`). This is not reported.
- The new settings A6 and A7 (female snRNA and female microdissection confirmation) give 222 and 167
  mouse-only genes. A7 retains 149 of 247 (60%) (`class_retention_by_setting.csv`). The draft
  reports neither.

**Fix.**
1. Require species-only genes to be flat, or at least not same-sign significant, in that species'
   external atlas.
2. Add A3b, scaling the human flat margin by the noise-corrected amplitude ratio.
3. Show the mouse-gradient distribution of mouse-only genes. How many have |mouse| ≥ 1.5 with human
   inside ±0.3?
4. Report A6/A7 and the S3 − early class counts.
5. Correct the Introduction sentence.

### M4. Species-level confounders that external human data share are not addressed; "tissue state ruled out" overclaims

**Location.** R3 "Ruled out: log compression and tissue state"; Discussion "Why flattening...";
Limitations.

**Problem.** The relabeling control cannot see donor-level or species-level terms, as the paper
correctly states. External replication removes only those confounders that differ between our human
tissue and Lake/Census. Several do not differ:
- **Procurement:** surgical or biopsy tissue with warm and cold ischaemia, against perfused or rapidly
  harvested mouse kidneys.
- **Age:** older patients against young-adult mice. The Census mouse dataset spans the lifespan, so
  donor ages should be stated.
- **Disease and medication.**
- **Physiological state.** Fasting before surgery against ad-libitum chow, and time of day. The
  mouse-only programs are precisely the classic fasting-, diet- and circadian-responsive ones:
  PPARα-driven peroxisomal β-oxidation, SREBP2 sterol synthesis and glutathione synthesis. Xiong
  et al. 2023 (already cited) report that caloric restriction feminizes the male mouse kidney.
  Mouse–human differences in PPARα target regulation are documented (Rakhshandehroo et al. 2009,
  https://pubmed.ncbi.nlm.nih.gov/19710929/).

The tissue-state evidence is also weaker than stated.

**Evidence.**
- Lake human-axis amplitude falls with injury: the median is 0.69 / 0.50 / 0.35 in healthy, CKD and
  AKI donors, and ρ = −0.41 (`lake_donor_scores.csv`, recomputed).
- Workstream 3's independent check concludes "supports tissue-state flattening"
  (ρ = −0.42, 95% CI −0.66 to −0.10; `pt_literature_deep/zonation_check/summary.json`).
- Mouse AKI cuts S3 − early amplitude to 0.41–0.44 of control on mouse axes, and the S3 − early
  AKI–control correlation to 0.59 (`aki_amplitude_ratios.csv`;
  `posthoc_reference_free_amplitude.csv`). This is close to the human/mouse S3c − early ratio of
  0.28–0.33. The draft quotes only the S2 − S1 value (0.74, r = 0.83).
- The tissue state of our own donor is unknown.

**Fix.**
1. Retitle the subsection, for example "Tissue state does not explain the S2 − S1 pattern".
2. Report the S3 numbers.
3. Score our two sections on the Lake injury axis to place the donor among the Lake donors.
4. Add a Limitations paragraph on procurement, age, fasting/diet/circadian state and
   disease/medication, noting that external human atlases share them.
5. Where possible, check whether the mouse-only programs keep their zonation in public mouse kidney
   data from fasted mice or other times of day.

### M5. Coordinate selection was outcome-driven, and the specificity result is coordinate-dependent

**Location.** R1 "Alternatives"; Methods 6.10; Abstract ("as the split-plot design predicts");
Discussion point 1.

**Problem.** The pre-registered paper-coordinate selection rule (`ws1/plan.md`, stage-2b protocol)
**selected DPT13**. In `nb36/selection_rule.csv`:
- DPT13 qualifies (P5 0.80, registration gap 0.14);
- SCF13 fails (P5 0.72 < 0.80, gap 0.29 > 0.15).

DPT was better on every pre-specified coordinate-quality metric: P1 human 0.29 against 0.25, P1w,
gene-fold agreement 0.59 against 0.49, P5, and gap. The downstream specificity gate was then added
because DPT failed it, no candidate qualified, and the fallback (scFates) was kept.

The draft says only that the gate "was applied to all candidates after DPT's result was seen". It
should say that the pre-registered rule chose DPT.

The main R2 claim, that position-dependent screens are specific, fails on DPT: NCDR 0.40, and 36
mouse-against-mouse positional pathway calls against 16. So specificity is a property of the
coordinate, not something the split-plot design guarantees.

The DPT result is also for the matched-only test. The final joint-test pipeline (notebook 37) was
never run on DPT13 or SCF13-ED (only `pt_pathway_final/scfates/` exists). R4 then cites "an earlier
analysis on ... DPT" as *supporting* robustness. A coordinate cannot be both too non-specific to use
and evidence of robustness.

**Fix.**
1. State plainly that the pre-registered rule selected DPT.
2. Run notebook 37 end to end with `--coordinate` DPT13 and SCF13-ED, and report joint-test NCDR,
   overlap with the 36, and within-species calls.
3. Make the coordinate-free positional screen the primary demonstration of specificity. The joint
   test on S1/S2/S3 steps already gives 16 calls for species and 0 and 0 for the relabelings
   (`table_S1_strategy_specificity.csv`). Present the coordinate as adding power, not specificity.
4. Drop DPT from the R4 robustness list, or explain the logic.

### M6. The specificity claim rests on two null draws, and the split-plot argument does not predict it

**Location.** R2; Abstract; Discussion point 1; Methods 6.14.

**Problem.**
- With four specimens there are exactly two null partitions. NCDR = 0.02 is a ratio estimated from
  two draws. It was also used to *choose*:
  - the strategy (12+ strategies; `table_S1`);
  - the test (the matched-only test was replaced by the joint test after the matched test's
    relabeling behaviour was seen);
  - the coordinate (M5);
  - the success criterion (the 5% condition was added after the first run, and decoy FDP was
    replaced by NCDR).
  Reporting the winner's NCDR on the same two draws is optimistic.
- The split-plot df argument (2 against 2(k − 1)) concerns specimen-level tests. T_spatial is a
  structure-level partial F used for ranking, so no specimen-level df enters the pathway test.
- The proper specimen-level split-plot test:
  - calls **nothing** with its pre-specified pathway test;
  - ranks genes only weakly like T_spatial (ρ = 0.21);
  - recovers only 15 of the 36 robust pathways with the post-hoc matched-only variant.
- The specificity actually observed is empirical. Specimen deviations in this dataset happen to be
  mostly offsets: within-species positional calls are 16 and 0, and decoy ρ is 0.009 for position
  against 0.019 for offset. On DPT this no longer holds.

**Fix.**
1. Describe the split-plot argument as motivation, not prediction.
2. Add an independent null with many draws. Options with existing data:
   - **(a)** Repeatedly subsample 2 + 2 designs from Census mouse (12 male donors) and Lake human
     (7 donors) segment pseudobulks. Run average-level and positional screens on species splits and
     mixed relabelings. Use the full-donor analysis as truth to estimate the real false-discovery
     proportion. This would turn R2 into a general result.
   - **(b)** Use the two AKI mice: {Ctrl1 + AKI1} against {Ctrl2 + AKI2} gives two more null draws
     with biological variance on both sides. AKI against control is also a 2 + 2 positive control
     with known S3-centred positional biology.
   - **(c)** Use spatial half-sections as pseudo-specimens.
3. Report the joint test's p-value calibration under relabeling, for example the fraction of
   pathways with p ≤ 0.01. At present only call counts at q ≤ 0.05 are given.

### M7. The pathway layer mostly restates mouse zonation, and its external replication is not specific

**Location.** R3 "Position-dependent pathways recover both the conserved core and the
mouse-specific programs"; R4 "External agreement"; Abstract (the 60 pathways).

**Problem.** A species × position interaction cannot "recover a conserved core". A conserved gene
enters T_spatial only through an amplitude difference. Under the global component, the interaction
is −0.81 to −0.90 × the mouse gradient (R² 0.60–0.82), so T_spatial is large for any gene zonated in
mouse. The robust list is therefore largely a list of mouse-zonated programs.

**Evidence.**
- The draft itself says 32 of 36 are consistent with "mouse zonation, flat human".
- Only 1–3 pathways replicate once the global component is removed.
- **Robust pathways replicate only marginally better than relabel-called pathways**, which are
  specimen artefacts by construction. Mann–Whitney p is 0.016 (Census, male mice), 0.033 (Lake,
  male), 0.19 and 0.21 with female mice (`external_replication_summary.csv`). The text reports only
  the comparison with uncalled pathways (p = 10⁻¹³).
- The "dominant class" of each robust pathway is a plurality among 3–15% of its members. 60–80% of
  members are indeterminate, and two pathways are labelled "conserved" on a tie or one-gene margin,
  for example Estrogen Response Late at 7.96% against 7.96% (`robust_pathways_by_class.csv`).
- Two of the four "beyond flattening" pathways, Glutathione metabolism and Formation of Cornified
  Envelope, lose significance with equal-probe genes: q_probe_equal 0.24 and 0.22. So do the
  xenobiotic sets (Chemical carcinogenesis, Metabolism of xenobiotics by CYP450, Aspirin ADME,
  Retinol metabolism and Biological Oxidations; `table_S2_pathways.csv`).
- The Mann–Whitney p-values (7 × 10⁻²⁰, 2 × 10⁻¹³) treat overlapping, selected pathways as
  independent. The robust set was selected on T_spatial, which correlates with both comparison
  statistics, so part of each separation is built in.

**Fix.**
1. Reframe the robust list as "programs whose mouse zonation human PT does not share".
2. Report the relabel-called comparison.
3. Drop or qualify the "conserved core" reading of R3's pathway paragraph.
4. Replace the pathway-level Mann–Whitney p-values with program-level or permutation-based
   comparisons.
5. Flag probe-sensitive pathways in the text.

### M8. Several literature citations contradict the claims they are cited for

**Location.** Introduction (Breljak); R3 "This agrees with the conserved localisations of SGLT1/2,
AQP1 and OAT1–3"; R5 OAT paragraph; Discussion "Where human PT places apical transport is well
predicted by mouse".

**Problem.** Breljak et al. 2016 (https://pubmed.ncbi.nlm.nih.gov/27053689/; abstract verified)
report that in human kidney "all three OATs were stained more intensely in S1/S2 segments compared
with S3 segment in medullary rays". Our data disagree on two of three:
- human OAT1 (SLC22A6) is *flat*, giving a mouse-only class;
- human OAT2 (SLC22A7) *rises* toward late PT and is a conserved late marker (human S2 − S1 +1.93,
  conserved in S3c − early);
- OAT3 is indeterminate.

AQP1 (Maunsbach et al. 1997, verified) is indeterminate in our classes, not conserved
(`gene_zonation_classes.csv`). Only SGLT1/2 agree.

So protein-level human data contradict two of our transcript-level transport calls, and "mouse
predicts transport" is overstated:
- 18 zonated transporters are mouse-only, including OAT1, NaDC3 and SMCT2 (Slc5a12);
- conserved genes are 56% of zonated transporters, against 30% of all zonated genes **[reviewer
  check]**.
- Human and mouse also use different phosphate transporters (SLC34A3 contributes about 40% in adult
  human and is negligible in mouse; https://pubmed.ncbi.nlm.nih.gov/42494142/), even though
  Slc34a3 zonation is "conserved".

**Fix.**
1. Cite Breljak as a *discordance* and discuss mRNA–protein and antibody-specificity explanations.
2. Remove AQP1 and OAT1–3 from the "agrees with" sentence.
3. Soften "Mouse PT predicts where human PT places transport" to a quantified statement.
4. Distinguish conserved zonation from conserved transporter usage.

---

## B. Minor issues

**R1 (coordinate)**
1. *Table, "Mouse | Male snRNA | Coordinate 0.75".* The fold-protected value is 0.745 (rounds to
   0.74; `nb35/p1_external_concordance.csv`). 0.75 is the unprotected value.
2. *"without changing which pathways were called".* The circular arm called 145 pathways and the
   independent arm 130; robust retention was 57 against 59 (`g6_count_split.csv`). Say "without
   changing the robust list materially".
3. *"Removing a third of the genes ... moved human segments by up to 0.29".* In that refit (fold 2)
   the human S1 < S2 < S3 order also broke (`species_registration_across_refits.csv`). Report it.
4. *"DPT was more stable".* DPT was better on every coordinate metric (see M5); "more stable"
   understates this.
5. *"It also localises differences continuously".* This conflicts with R1's own finding that
   within-segment order is weak (P1w 0.10 in human). Peak positions (0.16, 0.63) should be read at
   segment resolution only.
6. *Depth.* Human S1 ordering correlates with library size (−0.31, −0.38). Pathway calls on the
   equal-depth refit (SCF13-ED) are not reported. Run them, or add log library size as a covariate
   in the gene models.

**R2 (specificity)**

7. *"No average-level screen reached a decoy-estimated FDP of 10% at any threshold".* False for the
   pathway-score offset model (family "average"; 1,087 calls at FDP ≤ 10%). Say "No competitive
   average-level screen".
8. *R2 table mixes tests.* The T_level row is the matched test (NCDR 3.85) and the T_spatial row is
   the joint test (0.02). Give a 2 × 2 table of {T_level, T_spatial} × {matched, joint}: matched
   3.85 / 0.13; joint 0 calls with 46 and 2 relabeled / 60 calls with 2 and 0.
9. *Probe panels.* The strongest reason average-level species comparisons are uninterpretable is not
   specimen number. Per-gene probe efficiency differs between the human and mouse panels.
   Within-species gradients cancel it; species offsets do not. DESeq2's 9,010 "specific" genes are
   specific only in the relabeling sense. State this.
10. *"Specimen differences are offsets".* The two human sections are not biological replicates, so
    "0 positional calls" between them says little about donor-level positional variation.

**R3 (amplitude and classes)**

11. *Deming and orthogonal ratios.* On the Lake axis they are 0.16–0.17, not ~2
    (`amplitude_ratios.csv`, `ratio_deming`, `ratio_orthogonal`). Figure S4a shows them, but the
    text does not explain why they contradict the bracketing reading. With weak correlation, the
    total-least-squares slope tracks variance ratios rather than covariance.
12. *"Human gradients ... about half as large overall"* (Discussion). This holds only for S2 − S1
    (0.47–0.58). S3c − early is 0.84.
13. *Sex check.* It considers male-biased genes only. Mouse-only genes are 14.2% female-biased,
    against 5.5% of indeterminate and 3.9% of "neither" genes **[reviewer check]**. Female-biased
    genes are androgen-repressed in males and can also produce male-specific zonation. Report both
    directions against the background of indeterminate genes, and report A6/A7.
14. *Class enrichments.* Pyrimidine metabolism, Azathioprine ADME and Mucopolysaccharidoses
    (mouse-only), and the human-only "Early Phase of HIV Life Cycle" and "Mismatch Repair", are
    not robust to expression or rest on one gene (Lrig1) (`class_programs.csv`). Remove them from
    Table 1's "top programs".
15. *Figure 3c.* The mouse-only representative Slc5a12 is an S1 labelling and orientation marker.
    Choose a gene not used to define labels.
16. *Denominators.* "Largely species-specific" describes 520 classified zonated genes (7% of 7,407).
    The abstract's list omits the 1,297 "neither" genes. State both denominators.
17. *Confirmation region mismatch.* Mouse S3c calls are confirmed against references whose S3 is
    mostly outer stripe. Note this.

**R4 (robustness and replication)**

18. *"Restricting to genes detected in both species kept 94% ... removing mouse sex-biased genes
    kept 97%".* These are fractions with effect ratio ≥ 0.7 (34 and 35 of 36, counting non-testable
    as failures). The actual passes (ratio ≥ 0.7 and q ≤ 0.10) are 28/36 and 32/36, which is what
    produces 26 core pathways (`table_S2_pathways.csv`). Report the pass counts.
19. *"Every robust pathway keeps at least 70% ... equal probe counts".* True for effect, but 8 of 36
    lose significance (see M7). Say so.
20. *Conventional-only pathways.* "Replicated 6 times ... never at BH ≤ 0.10" omits that, as a group,
    they replicate better than other pathways (Mann–Whitney p = 7 × 10⁻⁴;
    `external_average_level_conventional_only.csv`).

**R5 (genes)**

21. *"In the human atlas both are nearly flat (ACADM +0.46; ACAA2 +0.58)".* +0.58 exceeds the paper's
    own 0.5 log2 zonation floor. Lake male donors give +0.82 and +0.95 (`story_gene_evidence.csv`).
    ACAA2 may be a reversal (notebook 39 calls it one). Do not call it flat.
22. *Gatm human.* Report our own human value (−0.78) and Census human (−0.10). The evidence table
    flags "human_references_disagree".
23. *"Of 39 claimed genes, 30 met every pre-specified test" and the Discussion's "checked list ... under
    pre-specified rules"* (adjacent to the known reversal-label item). These are inconsistent with
    R3:
    - six notebook-38 "headline reversals" are not reversals under R3;
    - Acadm, Acaa2, Acsm3 and Nudt19 are headline but indeterminate in our data;
    - the sterol genes verified in notebook 38 (Hmgcs1, Hmgcr) are indeterminate or conserved,
      while R3's mouse-only sterol genes (Lss, Ebp) were never verified gene by gene;
    - the glutathione story and AOX1 were "supporting" under the rules;
    - peroxisomal Acox1 and Hsd17b4 were not tested.
    Rewrite the Discussion list to match what was verified.
24. *"All five external references agree"* (Dcxr). Six values are listed.
25. *UGT1A9.* UGT1A isoforms share exons 2–5. Show that both panels' UGT1A9 probes target the unique
    exon 1. Otherwise UGT1A9 calls may measure the UGT1A locus.
26. *RBP4.* Besides perirenal-fat ambient RNA, Visium HD lateral diffusion and megalin-mediated
    protein uptake (relevant to the HPA protein signal) should be discussed.
27. *Coordinator language.* Remove "coordinator decision" and "coordinator promoted it" from the
    Results. Use "a post-hoc decision by the authors". The promotion rationale for drug handling is
    weaker than stated: the Glutathione metabolism pathway fails the equal-probe check.

**Methods and statistics**

28. *Count splitting* (6.10) needs a citation: Neufeld et al. 2023, Biostatistics,
    https://pubmed.ncbi.nlm.nih.gov/36511385/.
29. *Cross-species integration.* Harmony with species-confounded batches needs justification and a
    citation to a benchmark of cross-species integration (Song et al. 2023, Nat Commun,
    https://pubmed.ncbi.nlm.nih.gov/37838716/). θ = 6 should be justified.
30. *Positional-DE comparators.* tradeSeq (https://pubmed.ncbi.nlm.nih.gov/32139671/) and condiments
    (https://pubmed.ncbi.nlm.nih.gov/38280860/) are the standard tools for condition × pseudotime
    tests. Cite them, and explain why they were not used.
31. *Pseudobulk debate.* The "conservative but underpowered" characterisation of pseudobulk
    (Zimmerman 2021) is contested (Murphy & Skene 2022,
    https://pubmed.ncbi.nlm.nih.gov/36550119/). Cite both.
32. *Visium HD aggregation.* Bin-to-structure aggregation and diffusion should cite prior Visium HD
    work, for example Bin2cell (https://pubmed.ncbi.nlm.nih.gov/39250728/). Add a per-species
    contamination index: the mean of non-PT markers (Umod, Slc12a1, Aqp2, Slc12a3) in PT structures,
    relative to their home structures. HUK1_MED1 has only 7,223 of 14,186 polygons with bins, and
    unrecorded H&E–Visium registration could dilute human gradients.
33. *"Dilution-free ratio".* When cross-species correlation is low, this ratio is a projection
    ratio, not an amplitude ratio. Name it accordingly.
34. *Literature check.* The check was done with AI agents (6.22, [VERIFY]). Journals will require
    this to be stated plainly. Describe the search as a structured database search with defined
    queries, and archive the queries.
35. *Kim et al. 2011.* A label-retaining-cell study is a weak primary source for S3 vulnerability,
    although its abstract does state the finding. Prefer a dedicated source; Heyman et al. 2010 is
    already cited.

**Clarity and structure**

36. *The method and the biology are entangled.* R3 and R5 use neither the coordinate nor the
    pathway screen. Since public snRNA alone reproduces the central claim (r = 0.23), the paper must
    say what the spatial data add: intact tissue, region matching by depth (which overturned the OAT
    reading), and replication across structures within a specimen. Restructure so that either:
    - the biology (R3) leads, with R1/R2 as methods; or
    - the paper is a methods paper with the biology as an application.
    In the current order (R1 coordinate → R2 specificity → R3 labels-only biology → R4 pathway
    robustness → R5 genes), the pathway story is split around a section that does not use it.
37. *Jargon density.* NCDR, T_spatial, T_level, joint test, P1/P1w, S3c, "confident", "robust",
    "core", decoy FDP, A1–A7. Add a short glossary box, and drop internal names (SCF13, CAC) from
    the main text.
38. *Abstract.*
    - "as the split-plot design predicts" (see M6) should go.
    - "Mouse PT predicts where human PT places transport, but not metabolism" (see M8) should go.
    - Bullet lists are not accepted in most journal abstracts.
    - The abstract should name the single-donor limitation in the sentence that reports the
      finding, not only in "All results are descriptive".
39. *Title option 4* ("Mouse ... zonates metabolic programs that human ... keeps flat") overclaims.
    Most metabolic leads are indeterminate in our data (Acadm, Acaa2, Acsm3, Crot, Nudt19).

---

## C. Claims against evidence: numbers checked

**Matched source tables.** Abstract and Results values checked and correct:
- 12,866 PT structures; per-specimen and per-segment counts; 12,272 in common support; 26,839 and
  37,559 structures.
- S2→S3 transitions 0.63–0.64; S1→S2 transitions 0.26–0.28 and 0.35–0.38.
- P1 values: 0.53, 0.25, 0.73, 0.51, 0.27 and the unprotected 0.31. P1w 0.17 and 0.10. Gene-fold
  agreement 0.34–0.70. P5 0.36–0.46; library correlations −0.31/−0.38 → −0.04/−0.13. Gap 0.29;
  numerical shift 0.44; n_map SD 0.007 and ρ 0.994.
- Absorption 1.01; 91% of 244; 59/66 and 39/42.
- DPT NCDR 0.40, 36 against 16 within-mouse calls, 96 calls, 56/66 and 38/42; CAC 42/66, 30/42 and
  7,167 structures.
- R2 table rows; minimum FDPs 1.00, 0.75, 0.95 and 0.96; DESeq2 9,010 / 2 / 8; decoy a and ρ;
  1,444, 1,509 and 67–1,317 pathway-score calls; within-species 48/16 and 94/0; joint test 60 / 2 / 0
  and 57 (23–78); 59.5 against 15; split-plot 0.6% / 0.4% / λ 0.87, a = 1.10 against VIF 2.07, 41
  calls, 15/36, z 2.6 against −0.5, p 7 × 10⁻²⁰.
- Amplitude ratios 0.15 (0.11–0.19), 0.15, 0.16, 0.23 and 1.96 (1.41–3.05); 0.14 and 1.66;
  disattenuated r 0.14 (0.10–0.19), 0.23 and 0.83; SD ratios 0.58, 0.47 and 0.84; tertile 0.18;
  probe-balanced 0.16; AKI 0.74; ρ −0.24 and −0.41; injury halves 0.18 and 0.19.
- Classes 155 / 52 / 247 / 66 / 1,297 / 5,590; A3 243 / 141; 219 / 44; 48 / 35; 1,622 and 86;
  sex-bias 5.3%, 15.5% and 199.
- Funnel 60 / 58 / 51 / 36 / 26; 17 programs; peaks 0.16 and 0.63 (10 and 26); shapes 26 / 4 / 6;
  averaging loss 0; 27 of 36 and 17 of 26 in the earlier list.
- External: 17 and 6 of 36; global slope −0.81 to −0.90, R² 0.60–0.82; 1–3 beyond the global
  component; gene table 0.67 / 0.57 (0.62) / 0.64 and 0.73 / 0.60 (0.66) / 0.72; Mann–Whitney
  2 × 10⁻¹³ and 1 × 10⁻⁹; 184 and 6; Lake genome-wide 0.43 / 93.1% / 335 and 0.70 / 96.7% / 276,
  0.67 / 93.7%.
- R5: 39 / 30 / 5 / 4; a = 0.86 (0.64–1.35); Gatm 6.8 / 6.3 / 12.1 / 11.0, 23.5 against 4,270 TPM,
  +0.51; Acadm +2.1/+2.3; Acaa2 −4.0/−4.2; ACADM +0.46; ACAA2 +0.58; UGT1A9 −0.67 and ≤ 0.8 TPM;
  GCLC +0.61 against +3.3/+2.9; Dcxr +1.53 (+1.59/+1.48), −0.46/−0.54, −1.94/−2.03, rat −1.45;
  RBP4 +3.46; AOX1 +1.13; OAT1 970 TPM.

**Mismatches or misleading presentation.** Each is detailed above.

| # | Location | Draft says | Source says |
|---|---|---|---|
| 1 | R1 table | fold-protected mouse snRNA 0.75 | 0.745 |
| 2 | R1 count split | "without changing which pathways were called" | 145 against 130 calls |
| 3 | R1 alternatives | "DPT was more stable" | pre-registered rule selected DPT (`selection_rule.csv`) |
| 4 | R2 | no average-level screen reached FDP 10% | pathway-score offset: 1,087 |
| 5 | R4 | 94% and 97% kept | pass counts 28/36 and 32/36 |
| 6 | R4 | ≥ 70% effect under equal probes | true, but 8/36 lose q ≤ 0.10 |
| 7 | R3 | tissue state "ruled out" | Lake human-axis ρ −0.41; WS3 "supports tissue-state flattening"; mouse AKI S3 ratio 0.41–0.44 |
| 8 | R3 | correlation 0.14 (S2 − S1) only | S3c − early 0.29 not reported |
| 9 | R3 | "every alternative rule" | true for A1–A7; S3 − early contrast is 77 against 111 |
| 10 | R5 | ACAA2 +0.58 "nearly flat" | above the paper's own 0.5 floor; Lake male +0.95 |
| 11 | R5 / Discussion | "checked ... under pre-specified rules" (sterol, glutathione, peroxisomal, AOX1) | see minor 23 |
| 12 | Discussion | human "about half as large overall" | S3c − early 0.84 |
| 13 | R3 / Intro | "agrees with ... AQP1 and OAT1–3" | AQP1 indeterminate; OAT1 mouse-only; OAT2 conflicts with Breljak |
| 14 | R1 table | rat omitted | rat concordance 0.10–0.13, 55% direction agreement |
| 15 | R4 external | robust ≫ uncalled | robust against relabel-called: p 0.016–0.21 |

---

## D. Missing analyses, ranked by importance (feasibility with existing data)

1. **Anatomical validation of the human segment labels** (M1). Glomerulus-attached PT polygons;
   distance to glomerulus by label; species-separate clustering; R3 without integration genes;
   S3c − early as co-primary. *Feasible now:* polygons and glomerulus labels exist (notebook 33 has
   the depth code).
2. **Same-species cross-dataset ceiling for the cross-species correlation, and a decision on the rat
   data** (M2). *Feasible now* from saved external contrasts and notebook 39's machinery.
3. **Symmetric classification** (M3). External confirmation of flat calls; a scaled flat margin;
   amplitude-matched classes; A6/A7 reported. *Feasible now:* notebook 40 settings.
4. **Notebook 37 end to end on DPT13, SCF13-ED and segment steps**, with the pre-registered outcome
   stated (M5). *Feasible now:* `--coordinate` exists.
5. **An independent, many-draw null for R2** (M6). Subsampled 2 + 2 designs from Census mouse and
   Lake human, and mixed control/AKI mouse relabelings. *Feasible:* donor-level pseudobulks exist;
   the AKI mice exist in the mouse-only pipeline.
6. **A positive control for the positional screen.** AKI against control mice (2 + 2), where S3
   injury programs should be called. *Feasible.*
7. **Tissue state and measurement quality of the human sections** (M4, minor 32). Injury score
   placed on the Lake axis; a spillover/contamination index per species; registration QC for the
   human polygons. *Feasible.*
8. **Programs most exposed to physiological state.** Check whether mouse-only sterol, PPARα and
   glutathione zonation holds in public kidney data from fasted mice or other times of day.
   *Partly feasible* (public data); otherwise discuss.
9. **More human donors with spatial data.** Lake et al. 2023 / KPMP Visium and Slide-seq human
   kidneys at coarser resolution: segment-level contrasts from spatially mapped PT spots in
   additional donors. *Feasible with public data, with lower resolution.*
10. **Protein or in-situ validation of position-defining genes** (GATM, GCLM, DCXR, UGT3A1, RBP4,
    UGT1A9, SLC22A6/7). HPA images give presence only. Rat segment proteomics (Limbutara et al.)
    could be checked for the rodent side. *New experiments for human in situ.*

---

## E. Is the central claim supported?

**Partly.** The pieces that hold:
- the cross-species correlation of segment gradients is low in our data and in independent atlases;
- the canonical apical transport markers form the clearest conserved set;
- mouse-only calls outnumber human-only calls under every variant of the class rules.

The pieces that have not been established:
- that the low correlation reflects species rather than human label definition (M1) or dataset
  heterogeneity (M2; the mouse–rat value is as low);
- that "flat in human" is a population property and not a property of one donor's two sections or of
  a global amplitude reduction (M3);
- that the mouse-only programs are not state-dependent rodent programs shared by all lab-mouse
  references (M4).

The strongest alternative explanations a referee will raise, in order:
1. The human S1/S2 boundary is an integration artefact.
2. Cross-species and cross-dataset correlations are generally this low (rat).
3. Asymmetric flat calling, together with global attenuation, inflates the mouse-only class.
4. Procurement, age and fasting/diet/circadian state differ between humans and lab mice.

The draft addresses region, probe count, log compression, injury (for S2 − S1 only) and male-biased
genes. It does not address alternatives 1, 2 or 4, and addresses 3 only through the effect floor.

A defensible version of the claim today would be:

> In segment-labelled data from one human donor and in public atlases, the S1→S2 gene gradients of
> human and mouse PT are weakly correlated, and a set of canonical apical transporters is the main
> shared component. Whether this reflects species or label definition remains to be tested.
