# Glucose handling along the PT (notebook 92)

Frozen on 2026-10-15, before any statistic of notebook 92 was computed. The user reported from a meeting with kidney
collaborators: "they mentioned slg2 inhibitors drugs, how we currently don't understand how thye work as a miracle
drug and could be improved through this method, how glucose changes through nephon. this is extremly hot topic
currently".

**What the data can and cannot address.** No specimen here was treated with an SGLT2 inhibitor or was diabetic. The
notebook can therefore do four things:
1. Map where along the healthy PT glucose is taken up, exported and made: two control mice, one human donor.
2. Read what ischemic injury does to that map from notebook 91's screen: two IR mice against two controls.
3. Run a glucose balance driven by the measured profiles. It is an illustration, and it asks whether a continuous
   profile changes what a stepped one predicts.
4. Size the drug experiment from these mice.

It cannot say how SGLT2 inhibitors act. That needs treated tissue.

## 0 · What was seen before freezing

- **Panel membership and PT-wide detection** of 33 candidate genes per specimen (a scratch feasibility script). These
  are whole-PT shares; no positional profile of any pre-named gene was viewed in this session.
  - SGLT2 (Slc5a2) is detected in 0.69–0.70 of mouse PT units and 0.72–0.77 of human units.
  - SGLT1 (Slc5a1) is detected in 0.41 of mouse units and 0.03 of human units.
  - MAP17 (Pdzk1ip1) is detected in 0.96–1.00 of units in both species.
  - GLUT2 (Slc2a2) is detected in 0.59–0.73 of units.
  - In IR against controls: G6pc 0.42–0.43 vs 0.70–0.78; Pcx 0.70–0.76 vs 0.91–0.95; Bdh1 0.20–0.21 vs 0.48–0.51;
    Hmgcs2 0.17–0.19 vs 0.08–0.09; Slc9a3 0.47–0.50 vs 0.62–0.66.
  - G6pc pairs with no human gene in the shared export (the human panel names it G6PC1). LDHA is not on the human
    panel.
- **Earlier results that touch these genes.**
  - Notebook 90 reported Slc5a2 as continuous (ΔBIC 42) and among the genes conserved between species.
  - The expert reading before notebook 91 named Hmgcs2 among fasting genes up in IR.
- **Which genes built the coordinate.** These are the integration HVGs of either export plus 90's orientation genes.
  - Slc5a2, Slc5a1 and Slc2a1 helped build both coordinates (Slc5a2 also oriented the labels).
  - Pdzk1ip1 (MAP17) and Slc2a2 (GLUT2) built neither, so they are fold-protected.
  - Following notebook 90's rule, any statement that rests on where a gene changes along the axis rests on
    fold-protected genes. SGLT2, SGLT1 and GLUT1 are reported, but their positions are partly circular.

## 1 · Literature (full-text reads, 2026-10-15)

Three reads were run, each from full text where it could be reached. Abstract-only sources are marked AO. Every PMID
was checked on PubMed.

**Position along the PT.**
- Rodent SGLT2 protein sits in the brush border of S1 and S2 (S1 > S2) and is absent from S3. Sources: Sabolic 2012,
  PMID 22262063, Figs 1–2 and 6–7; Vallon 2011, PMID 20616166, Fig 3.
- The two antibody papers disagree on how far SGLT2 reaches into S2.
- Rat SGLT1 runs S1 < S2 < S3 (Balen 2008, PMID 18524944, AO). Mouse SGLT1 is reported higher in S2 than S3, while
  human SGLT1 is strongest in S3 (Vallon 2020 review, PMID 32144488).
- Human SGLT2 is in S1/S2 and SGLT1 in S3 (Vrhovac 2015, PMID 25304002, AO).
- No quantitative axial profile of either transporter exists. Gorboulev 2012 (PMID 22124465) implies an overlap window
  in late S1 to S2.
- MAP17 activates SGLT2 about 150-fold without changing its surface amount (Coady 2017, PMID 27288013). Its position
  along the PT has never been mapped.
- GLUT2 is strongest in S1/S2; GLUT1 runs along the whole PT, highest distally (Vallon 2020).
- In rabbit, PEPCK is equal in early and late convoluted tubule and 50% in the pars recta, and hexokinase rises along
  the nephron (Vandewalle 1981, PMID 7246739, AO).
- mRNA is not protein: in mouse, female Sglt2 mRNA is about 49% higher while SGLT2 protein is about 2.4-fold lower
  (Sabolic 2012, Fig 7).

**Capacity.**
- In euglycaemic mice, SGLT2 reabsorbs 97% of filtered glucose and SGLT1 3% (Rieg 2014, PMID 24226519).
- Micropuncture in wild type: 78 ± 6% is reabsorbed by the early PT and 93 ± 1% by the late PT. In Sglt2−/− mice:
  0% by the early PT, 21 ± 6% by the late PT, and 36 ± 8% for the whole kidney (Vallon 2011, Figs 5–6).
- At high empagliflozin, mice still reabsorb 44 ± 3% (Rieg 2014).
- Human SGLT2 and SGLT1 have similar glucose affinity (K0.5 4.9 vs 1.8 mM; Hummel 2011, PMID 20980548). So the
  division of labour between them is set by abundance and position, not affinity.

**Models.**
- The Layton, Vallon and Edwards nephron models represent the PT as S1–S2 with SGLT2 + GLUT2 and S3 with SGLT1 + GLUT1
  (Layton 2015, PMID 25855513; Layton 2016, PMID 26764207; Hu 2021, PMID 34072329).
- Their densities are calibrated to whole-organ data rather than measured. No model uses a continuous axial profile.
- Under SGLT2 inhibition the models move reabsorption and O2 use to S3. Layton 2016 (diabetes, chronic inhibition)
  predicts S3 active QO2 +26%, rising to +34% if SGLT1 in S3 is not downregulated.
- The SGLT1 density along S3 is the uncertain parameter they name.

**AKI.**
- Clinically, SGLT2 inhibitors lower AKI:
  - serious AKI adverse events OR 0.64 (Menne 2019, PMID 31815931);
  - EMPA-KIDNEY HR 0.78 (PMID 36331190).
- In mice, dapagliflozin before IR protects via HIF1 (Chang 2016, PMID 27391020).
- But Sglt2−/− mice are not protected from IR (Nespoux 2020, PMID 32116018), and Sglt1−/− mice recover better
  (Nespoux 2019, PMID 30995111).
- After IR, whole-kidney Sglt2 mRNA falls about 57–60% at 2–3 weeks; Sglt1, Glut1 and Glut2 return to normal.
- No study reports where along the PT any of these changes sit.
- **Gluconeogenesis and glycolysis after injury:**
  - Pck1 and Fbp1 fall and Hk1 and Pfkp rise after IR. Patients with AKI switch the kidney from lactate uptake to
    release (Legouis 2020, PMID 32694833; read as the preprint full text).
  - Gluconeogenic enzymes fall with CKD (Verissimo 2022, PMID 35273087).
  - Glycolytic enzymes rise in cortex and outer stripe after IR (Lan 2016, PMID 27000065).
- **Ketones after injury:**
  - HMGCS2 is proximal-tubular, and both IR and sham surgery induce it to fasting levels.
  - Kidney Hmgcs2 loss worsens IR injury (Feola 2026, PMID 41568909).
  - Empagliflozin's renoprotection in diabetic mice needs Hmgcs2 (Tomita 2020, PMID 32726607, AO).

**SGLT2-inhibitor omics and novelty.**
- **Positional hypotheses tested only in bins or models.** Lower cortical and higher S3/outer-medullary workload and O2
  use (Layton 2016; Laursen 2021, PMID 34386735: cortical R2* falls, medulla unchanged).
- **The drug's footprint in healthy mice sits downstream.**
  - Eight weeks of empagliflozin in healthy mice enlarges cells only in SGLT2-negative S3 and collecting duct (Sinha
    2023, PMID 37033639).
  - Microdissected SGLT2-positive S1 from non-diabetic mice changed only 2 genes with dapagliflozin (Kim/Vallon 2025,
    PMID 39589189; GSE279174). Microdissected S3 from the same mice changed 838 genes: proliferation and stress up,
    inflammatory/NLRP3 down (Kim/Vallon 2026, PMID 42213650; GSE325466).
- **Single-cell studies resolve the PT by subcluster or S1 against S2.**
  - Schaub 2023 (PMID 36637914; human youth T2D): glycolysis, gluconeogenesis and TCA suppressed, mTORC1 lower.
  - Wu 2022 (PMID 35709763): S1 against S2 as the drug-specificity control.
  - Wu 2022 (PMID 34678510): S3 share restored.
- **Novelty.** The read found no study that maps SGLT2-inhibitor effects, or glucose-handling genes, on a continuous
  PT coordinate. No public SGLT2-inhibitor-treated spatial dataset exists in any species. The paired S1 and S3
  microdissection datasets are the closest external positions.

## 2 · Genes and sets (named before any profile was viewed)

| Role | Mouse genes (protein) |
|---|---|
| Uptake | Slc5a2 (SGLT2), Pdzk1ip1 (MAP17), Slc5a1 (SGLT1) |
| Exit | Slc2a2 (GLUT2), Slc2a1 (GLUT1) |
| Other sugar transporters | Slc5a10 (SGLT5), Slc5a9 (SGLT4) |
| Sodium coupling | Slc9a3 (NHE3), Atp1a1 |
| Production (gluconeogenesis) | Pck1, G6pc, Fbp1, Pcx, Aldob |
| Glycolysis | Hk1, Hk2, Pfkp, Pfkfb3, Pkm, Ldha |
| Lactate | Slc5a12 (SMCT1), Slc16a1 (MCT1), Slc16a3 (MCT4) |
| Glucose sensing | Txnip, Arrdc4, Mlxipl (ChREBP) |
| Ketones | Hmgcs2, Bdh1, Oxct1 |
| Insulin | Insr |
| Regulators | Hnf4a, Ppara, Ppargc1a |

**Gene rules.**
- Genes below detection drop out of the analyses that need them; the notebook reports them.
  - Mouse profiles need 90's usable rule (≥ 5% in both controls).
  - Injury readings need 91's shared rule (≥ 2% in all four mice).
  - Human profiles need ≥ 5% in both sections.
- Human symbols come from the ortholog map, with G6pc → G6PC1 by hand.

**Sets (library of notebooks 71 and 91).**
- Cellular Hexose Transport; Gluconeogenesis; Glycolysis/Gluconeogenesis; Glycolysis (Reactome and Hallmark); Glucose
  Metabolism.
- mTORC1 Signaling; Insulin signaling; AMPK signaling; Energy-dependent regulation of mTOR by LKB1-AMPK.
- Pyruvate Metabolism; Butanoate metabolism (ketone bodies).
- Hallmark Hypoxia; HIF-1 signaling; PPAR signaling.

## 3 · Analyses

**Labels.**
- Notebook: `analysis/notebooks/92_pt_glucose_handling.ipynb`; logic `92.glucose_handling.1`; seed 92.
- Results: `results/pt_glucose_handling/`.
- Module: `pseudospace/glucose_handling.py`, with `tests/test_glucose_handling.py`.
- Figures carry plots only.

**Inputs.**
- Notebook 90's human-vs-mouse construction, with its guards: 71's universe, 72's transitions, stage-01 native
  panels, and the usable genes equal to 90's ladder universe.
- Notebook 91's AKI construction, with its guards: cuts 0.339 and 0.657; the shared genes and L reproduce 91's.
- The sha256 of every table read from 90 and 91 is recorded here.

**H1 · Healthy mouse shapes.** The pre-named genes' rows of 90's ladder table are read, not refitted:
- class and continuous-vs-discrete reading;
- switch position and width, bump position;
- segment steps;
- fold-protected flag.

The centre of mass of each gene's mean control deviation profile is added.

**H2 · Where uptake hands over.**
- *Genes:* SGLT2, MAP17 and GLUT2 (falling) and SGLT1 and GLUT1 (rising).
- *Fit:* per control mouse, 90's switch model on the gene's deviation profile, giving a logistic midpoint and a 10–90%
  width.
- *Uncertainty:* 500 resamplings of the mouse's units, with 95% intervals. Neighbouring units can share a tubule, so the
  intervals are optimistic; the two mice are the replicates.
- *Exports:* the fit is repeated on both exports (human-vs-mouse and AKI) for the same two control mice.
- *Place:* a midpoint within 0.15 units of 1 or 2 is "at" that label boundary; otherwise it is inside a segment.
- *Gaps:* SGLT1's midpoint minus each falling gene's. Negative is an overlap of the two halves, positive a gap.

**H3 · Set positions.** 90's set positions (centre of mass against expression-matched random sets) for the pre-named
sets, all members and fold-protected members.

**S · Human.**
- 90's end contrasts E (late − early, level-free) for the pre-named genes in mouse and human.
- The human deviation profiles in human segment units (the human sections' own label transitions).
- 91's species gradient cosine C for the pre-named sets.
- Human SGLT1 is pre-registered as **not testable**: 3% detection, and the sections lack the outer medulla, where most
  of S3 runs.

**I · Injury.** Read, not refitted, from 91:
- level L, location M and the replication r_P, with their relabelings;
- call, peak position, injury shape class and location;
- gradient survival A and its class;
- induced or lost membership.

Sets come from 91's library: level and location z, q and calls, the set's survival A, and the set's rank among all
91's sets. Δ(x) is recomputed for the figure, with a guard that it reproduces 91's L.

**D · Experiment size.**
- *Per-mouse statistics (AKI export):* each shared gene's level ℓ (mean over valid bins of its log2 rate, centred per
  bin by the median reference gene) and location profile (the centred rate minus ℓ).
- *Noise:* the within-condition SD from the two controls and the two IR mice (2 df per gene). It is pooled in
  expression deciles by its median and corrected by χ²₂'s median. The location noise is the RMS of within-pair
  profile differences / √2.
- *Smallest detectable effect:* a two-sided two-sample t-test with 80% power at α 0.05 and 0.001, for 2–8 mice per
  group.
- *Smallest permutation p:* the number of other splits per design and the smallest p they allow.

**X · The drug's own footprint (GSE279174, GSE325466).**
- *Data:* TPM tables of microdissected SGLT2-positive early PT ("S1") and LTL-positive SGLT2-negative late PT ("S3",
  deep cortex and outer medulla) from the same male DBA mice. The groups are Sglt1 wild type or knockout, non-diabetic
  or Akita, and vehicle, dapagliflozin or semaglutide for 2 weeks. Only Sglt1 wild-type mice are used.
  - The S3 table is log10 TPM with blanks for zero. It is converted to linear and renormalised to 10^6 per sample over
    its own genes.
- *X0 (axis validation):* in the wild-type vehicle mice present in both tables:
  - on common genes, renormalised together, with mean TPM ≥ 5, take the mean over mice of log2(S3 + 1) − log2(S1 + 1);
  - correlate it (Spearman) with 90's mouse end contrast E, for all genes and for fold-protected genes, with a
    1,000-draw gene bootstrap;
  - 90's label step (D S3−S1) is reported beside it.
- *Programs:* within each table:
  - contrasts: dapagliflozin against vehicle in non-diabetic mice; Akita against wild type on vehicle; dapagliflozin
    against vehicle in Akita;
  - test: Welch t on log2(TPM + 1) for genes with mean TPM ≥ 5, BH within the contrast;
  - membership: q ≤ 0.05 and |log2FC| ≥ 0.5, split into up and down, giving 12 programs;
  - a program with fewer than 10 genes in a universe is not testable there.
- *X1 (healthy position):* each program's centre of mass on the healthy axis against 2,000 expression-decile-matched
  random sets (90's B rule). The p is the larger of the two mice's when their z agree in sign, and BH runs within the
  programs, for all members and fold-protected members.
- *X2 (injury):* each program's level L and location M in 91's decomposition against 2,000 matched random sets. It
  uses 91's level and location set rules with BH within the programs.
- *X3 (gene agreement):* for each table and contrast, Spearman ρ of the drug log2FC with the IR Δ at the late PT (bins
  ≥ 1.95) and with L, for the condition and both relabelings.

**M · Glucose balance (illustration).**
- *Densities:* relative SGLT2 and SGLT1 densities per bin are 2^(mean control deviation), scaled to a maximum of 1.
  MAP17 stands in for SGLT2 in the fold-protected version. The step version replaces each density by its mean within
  S1, S2 and S3.
- *Balance:* luminal flux falls by Michaelis–Menten uptake (K0.5 4.9 and 1.8 mM). Flow falls linearly to one third
  along the PT; every bin is the same length.
- *Calibration (each version separately):*
  - with SGLT2 absent, 40% of the load is reabsorbed at 8 mM (mice: 36% in Sglt2−/−, 44% at high empagliflozin);
  - 1% escapes at 16 mM, the threshold Layton's model is calibrated to.
- *Conditions:* plasma 8 and 20 mM, with SGLT2 inhibited 0%, 90% or 100%.
- *Variants:* no water reabsorption; K0.5 halved; densities below 5% set to 0; threshold 12 mM.
- *Readouts:* SGLT1's share of the load, the share taken up in S3, and the luminal concentration at S2/S3.

## 4 · Predictions and decision rules

Each prediction is scored as one of:
- **holds:** the call and its sign or place match;
- **direction only:** the sign matches without a call;
- **fails:** an opposite call;
- **not testable:** the gene is missing or below detection.

**Healthy mouse (H1–H3).** Direction comes from 90's end contrast E in both mice (≤ −0.25 falls, ≥ 0.25 rises).
- HP1. SGLT2, MAP17 and GLUT2 fall along the PT.
- HP2. SGLT1 and GLUT1 rise along the PT.
- HP3. The falling midpoints sit before the S2/S3 boundary (< 1.85), in both mice and both exports.
- HP4. MAP17 tracks SGLT2: midpoints within 0.3 units, in both mice.
- HP5. Gluconeogenic genes are higher early than in S3 (Vandewalle, AO). Hexokinase rises along the PT.
- Open questions, reported without a prediction:
  - the width of the SGLT2 and MAP17 fall (90's sharp < 0.5, graded 0.5–1, gradual > 1);
  - whether the halves of the uptake handoff overlap (literature: implied overlap, never measured);
  - whether mouse SGLT1 peaks inside S2 (90's bump class with its peak in S2; Vallon 2020) or keeps rising into S3
    (rat).

**Human (S).**
- SP1. SGLT2, MAP17 and GLUT2 fall along the human PT as in mouse: human E ≤ −0.25 in the mean of the two sections.
- SP2. A later human decline (Vallon 2020's contrast) is reported as the human against the mouse centre of mass. Each is
  in its own species' label units, so it is descriptive.

**Injury (I).** Predictions from the AKI read. Level holds with a 91 level call; location holds with a 91 location
call and the peak in the named stretch.
- IP1. Slc5a2 falls (level).
- IP2. Slc5a1 falls with its location peak in S3 (peak ≥ 1.95).
- IP3. MAP17 and GLUT2 fall with SGLT2: the same sign of L, and their Δ(x) correlates with SGLT2's at r ≥ 0.5.
- IP4. Slc2a1 rises, peaking in S3.
- IP5. Gluconeogenic genes fall (level), and their zonated members' gradients are attenuated (A < 0.6).
- IP6. Glycolytic genes rise, peaking at S2/S3 or in S3 (peak ≥ 1.75).
- IP7. Slc5a12 falls; Slc16a1 and Slc16a3 rise in S3.
- IP8. Hmgcs2 rises, more in S1/S2 than S3 (location peak < 1.95, or Δ late below L). Sham surgery induces Hmgcs2 too
  (Feola 2026), so a rise against unoperated controls cannot be put down to injury.
- IP9. Sets: Gluconeogenesis and Cellular Hexose Transport fall (level). Glycolysis and HIF-1/Hypoxia rise, with a
  location peak in S3.

**Model (M).**
- MP1. M-check. With SGLT2 absent, Sglt2−/− micropuncture places 21% of 36% (0.58) of reabsorption before the late
  proximal tubule. For each version, the model's share before S2/S3 is reported against 0.58.
  - "Late proximal" is the last surface loop of the convoluted tubule, so the comparison is approximate. No version is
    expected to win.
- MP2. Continuous against step. The method-relevant reading: does the continuous profile move SGLT1's share of the load
  under 90% inhibition, or the luminal concentration at S2/S3, by ≥ 10% relative to the step version in the base case?
  This is reported for every variant.

**Drug footprint (X).**
- XP0. Our axis agrees with microdissection. Spearman ρ of the S3 − S1 contrast with the end contrast E, on
  fold-protected genes, is > 0 with its 95% interval above 0.
- XP1. Program positions on the healthy axis are reported without a prediction. The drug's S3 programs are drug
  responses, not S3 markers.
- XP2. The non-diabetic dapagliflozin S3-up program is a level-up set in IR with a location peak in S3 (peak ≥ 1.95).
  Kim/Vallon 2026 find proliferation and stress up in S3 under the drug. Notebook 91 found IR proliferation sets
  S3-centred, and that result is known while this prediction is made.
- XP3. The non-diabetic dapagliflozin S3 log2FC correlates positively with the IR Δ at the late PT, and more than with
  either relabeling's.

**What would justify taking this to the collaborators.**
- The healthy map (HP1–HP4 on fold-protected genes) holds in both mice and both exports.
- The handoff is a measurable width that a step model omits.

Injury predictions that fail are reported as such. They concern IR, not the drug.

## 5 · Figures

- **92a, healthy glucose map:**
  - handoff genes' profiles per mouse with switch midpoints;
  - midpoints with intervals on both exports;
  - heatmaps of the pre-named genes in mouse and human;
  - end contrasts, mouse against human;
  - set positions;
  - centre of mass by role.
- **92b, injury:**
  - Δ(x) heatmap;
  - L and M with relabelings;
  - 91's sets with the pre-named ones marked;
  - survival A.
- **92c, experiment size:**
  - detectable change against mice per group;
  - smallest permutation p;
  - gene-level detectable change against the injury level change.
- **92e, the drug's footprint:**
  - microdissected S3 − S1 against the axis end contrast;
  - program positions;
  - programs in the injury decomposition;
  - the dapagliflozin S3 programs' IR location profiles.
- **92d, glucose balance:**
  - densities, continuous and stepped;
  - luminal glucose at 8 and 20 mM with and without inhibition;
  - SGLT1's share across variants.

## 6 · Frozen inputs (sha256)

| File | sha256 |
|---|---|
| 90 `a_gene_classes.csv.gz` | 48d3394334fd18c7fb7c6b7728bfec49a75303b74be49831b47ffffe6bdf156b |
| 90 `d_end_contrasts.csv.gz` | b4a5646f71f998fa241ef8fb05f30f9a2151a331ebfeb5ed58ac0ea80d9effe1 |
| 90 `b_set_positions.csv` | 907eed01e3b35653ee030bb0ce63593d5c525e8514988e79a574580e0f122f92 |
| 90 `b_set_positions_fold_protected.csv` | fd79f411f2ffd50f624f9ca3a90075e597f100cf55c6adaba9c4682c50d352dc |
| 91 `i_genes.csv.gz` | d6ae303b5139ea82c32aefeb6a61276f6d383ea0a562478c8e5159d70540e3b7 |
| 91 `i5_library_sets.csv` | 6502c1984ac2345452d718e779bac16f3966726a97cd1f21af4d9eee7ca4d181 |
| 91 `s2_species_sets.csv` | 414bd7d06aa662b8f58243778d622ddde768bbd56fe41aac6d25c05fa30b8636 |
| 91 `i_induced_genes.csv` | b3ae25a0393a6cd132ce57ab16ac1b7a03ed26fa587daca3076e7bd499baf6c3 |
| 91 `i_lost_genes.csv` | c4468348b0bad9b87f1e1bf45670560522e7cf6b0f89c8eafd991ba5017ae079 |
| 71 `gene_sets_aki_vs_control.json` | 034112512c80d4236497c5f7677829e99c67eafcd33ed645429b16cd93a0cdc4 |
| GSE279174 TPM table | c9c6cafb1f1367eb7f2e8486206b15790089431ab4151ea7718dcd362b87be51 |
| GSE279174 series matrix | eaacd8b08cecbb8a4947c118766aaca280f14e13ee7b956d92a2b7c749d1dc34 |
| GSE325466 log10 TPM table | 74f057f80e0b773d1c49842616734ba31e35c4c55785dd98890ee6ba448ae036 |
| GSE325466 series matrix | eb8ee419b4e9aebe233e54f694b9afba9edb373d88e460f987e4c62c9d84c72c |

The GEO files sit under `data/external/sglt2i_lcm_vallon/<accession>/`, outside Git. Before freezing, only their
headers, sample groups and column sums were read: the S3 table's linear sums run 0.98–2.25 × 10^6, hence the
renormalisation.

## 7 · Caveats

- Two control and two IR mice, one human donor, all male. Everything is descriptive.
- mRNA is not transporter protein or flux, and the sex inversion in mouse SGLT2 shows the two can diverge.
- The IR time point, harvest times and whether controls were sham-operated are not recorded here (lab metadata
  pending). Feola 2026 shows sham surgery alone induces Hmgcs2.
- The human sections are cortex only.
- The balance model's length scale is the label axis, not micrometres.
- The microdissected mice are DBA males treated for 2 weeks; ours are of unrecorded strain. The microdissected "S1" is
  SGLT2-positive early PT (S1 and part of S2), not a label-matched S1.
