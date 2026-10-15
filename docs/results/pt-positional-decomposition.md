# A discovery screen along the PT: level, location and gradient survival (notebook 91)

Frozen on 2026-10-15, before any statistic of notebook 91 was computed. After an expert reading of notebook 90, the
user asked: "for the purposes of a paper though, we need some type of filter or analysis that shows the interesting
biology that we can then explore more throughly. that is that thing, either pathway analysis or gene anylsis or some
for of continous satistic or anything that helped you find these things".

**What the expert reading found, and how.** It ran on 90's tables and the AKI profiles in scratch scripts that are not
in the repository. Recovering these findings is not independent evidence; they are named in advance (section 6) as
recovery checks. Anything else the screen finds is new.
- **S3 identity erased in IR.**
  - Cyp7b1 is detected in 94–99% of control S3 units and 1–3% of IR S3 units; Slc22a7 falls from 97–99% to 6–12%.
  - S2 genes (Slc22a6, Slc22a8, Cyp2e1) stay on at the end of the IR curve.
  - Across genes that did not build the AKI coordinate, the IR mice keep 46–60% of the control gradient
    (S2→S3: 24–35%).
  - Found by printing the label genes' profiles.
- **Injury reaches nearly every PT unit, but less in late S1.**
  - Kim-1 is in 60–100% of IR PT units, with a lower share in late S1 than in early S1 or S2 (both IR mice).
  - Found by tabulating hand-picked injury markers per window.
- **Whole-PT shifts that may not be injury.**
  - Clock genes are down (Dbp, Tef, Hlf, Nr1d1, Per1–3, Ciart).
  - Fasting and PPARα targets are up (Hmgcs2, Cyp4a14, Angptl4, Pdk4).
  - Male-biased genes are down (Slco1a1, Kap).
  - Found with hand-picked panels.
- **Species.**
  - Three urate transporters are among the 74 reversed genes, as are several drug targets (Mme, Hpgd) and two
    hypoxia targets (Egln3, Ddit4).
  - Found by reading 90's reversed list.
- **A null result.** Once the mouse gradient is accounted for, the mouse injury change in the late PT does not predict
  the human gradient. Found by a regression.

A paper needs the filter, not the reading: one analysis that scores every gene and pathway, with a stated null, and
whose ranked outputs are the candidates for follow-up. This protocol defines it.

## 1 · The decomposition

A condition effect along the axis, Δ(x), splits into three parts. Each part answers a biological question.
1. **Level L:** the mean of Δ over the axis. Does the whole PT change?
2. **Location P(x) = Δ(x) − L, with size M = RMS of P.** Where along the PT does the change sit beyond its average?
   A shape ladder (90's) classifies P as stepped at the labels, gradient, switch, bump or trough.
3. **Gradient survival A:** the share of a gene's reference gradient still present in the other condition. A ≈ 1
   kept, 0 erased, < 0 reversed. This is the only part a species comparison can use, because probe panels make levels
   incomparable.

Two filters then say why a level change happens:
4. **Kidney-wide:** does the gene change the same way in non-PT tubule units (DCT, thin limb) of the same mice? A
   change shared by other tubules points to a whole-animal cause.
5. **External state signatures:** does the change look like time of day, fasting or sex in public mouse kidney data?

Precedent: Lamian (PMID 37949861, read in full; `docs/research/continuous-analysis-literature.md`) splits a mean
shift from a trend difference. Gradient survival, the shape ladder on the effect, and filters 4 and 5 are added here.

## 2 · Labels

- Notebook: `analysis/notebooks/91_pt_positional_decomposition.ipynb`; logic version `91.positional_decomposition.1`.
- Results: `results/pt_positional_decomposition/` (`tables/`, `figures/`, `stage_cache/`, `run_record.json`).
- New module: `pseudospace/positional_decomposition.py`, with `tests/test_positional_decomposition.py`.
- Seed 91. Figures carry plots only.

## 3 · Part I · Injury (AKI export, four mice)

**Inputs (90's section F, unchanged).**
- The AKI export (pass-1 fingerprint `bb4a80c1b484`) and its support positions.
- The control mice's logistic transitions give segment units; guard: cuts equal 90's 0.339 and 0.657 to 3 decimals.
- U_PT units (`primary_identity` Proximal_Tubule), 0.1-unit bins on [0, 2.5), library = non-mt counts.
- Binned log2 rates (`pp.bin_profiles`, pseudocount 0.5, ≥ 20 units per bin). Valid bins are those valid in all four
  mice.

**Genes.** Non-mt genes detected in ≥ 2% of PT units in both mice of at least one condition.
- *Shared:* ≥ 2% in all four mice.
- *Induced:* < 2% in both controls, ≥ 2% in both IR mice.
- *Lost:* < 2% in both IR mice, ≥ 2% in both controls.
- *Reference genes* (for centring): detection ≥ 5% in all four mice (90's usable).

**Contrasts.**
- The condition: IR2A2 + IR2A4 against Ctrl1A2 + Ctrl1A4.
- The four IR × control pairs.
- Relabeling R1: IR2A2 + Ctrl1A2 against IR2A4 + Ctrl1A4.
- Relabeling R2: IR2A2 + Ctrl1A4 against IR2A4 + Ctrl1A2.

Δ_K(x) = mean case profile − mean reference profile, centred per bin by the median reference gene.

### I1 · Level and location (shared genes)

- **Per gene and contrast:** L, P(x) and M as in section 1.
- **Location replication:** r_P is the Pearson correlation between P from IR2A2 and P from IR2A4, each against the mean
  of the controls.
- **Calls:**
  - *Level change:* |L| ≥ 0.5; the four pairs share L's sign; both relabelings have |L_R| ≤ 0.5 |L|.
  - *Location change:* M ≥ 0.4; r_P ≥ 0.5; both relabelings have M_R ≤ 0.5 M.
  - A gene can carry both calls.
- **Injury shape.**
  - The model is 90's ladder (`fit_ladder`, `ladder_class`).
  - Replicates are the two IR mice's Δ against the mean of the controls.
  - Weights per bin are the harmonic mean of that IR mouse's units and the controls' summed units.
  - Label fractions come from the controls pooled.
  - Readout: classes and locations among location-change genes. The same ladder on each relabeling, with that
    relabeling's two case specimens as replicates, gives the null class counts.
  - "Discrete" here means an injury effect constant within segments and stepped at the labels: the segment model of
    injury.

### I2 · Gradient survival (shared genes)

- **Deviation profiles per mouse (90's):** profile minus its mean over valid bins, minus the median reference gene.
- **Zonated in controls:** the two control deviation profiles correlate at r ≥ 0.5, and the range of their mean is
  ≥ 1 log2.
- **A per gene, cross-fitted.**
  - Template control t; held-out control h; IR mouse i.
  - A_{i,t} = Σ_x d_i·d_t / Σ_x d_h·d_t.
  - The gene's A is the median over the four (i, t).
  - Using the held-out control keeps the controls' noise out of both numerator and denominator alike. A gradient
    estimated from the same controls would make injury look like flattening by regression to the mean.
  - Relabelings: the same, with each relabeling's case and reference specimens.
- **Gene classes:**

  | Class | A |
  |---|---|
  | reversed | < −0.25 |
  | erased | −0.25 to 0.25 |
  | attenuated | 0.25 to 0.6 |
  | kept | 0.6 to 1.5 |
  | amplified | > 1.5 |
- **Genome-wide survival:** the ratio of sums over zonated genes (Σ d_i·d_t over Σ d_h·d_t).
- **Per stretch:** the same ratio for the S1→S2 contrast (bins < 0.8 against 1.2–1.8) and the S2→S3 contrast
  (1.2–1.8 against > 2.1), for genes zonated in that stretch (control contrast ≥ 0.5 in both controls, same sign).

### I3 · Kidney-wide filter

- **Units:** DCT and Thin_Limb units of the same four mice in the AKI export. These are the two non-PT identities with
  ≥ 150 units in every mouse.
- **Per identity:** pseudobulk log2 CPM per mouse (non-mt library). IR mean − control mean, centred by the median of
  that identity's reference genes (≥ 5% detection in all four).
- **Gene call:**
  - *kidney-wide:* both identities share L's sign with |FC| ≥ 0.5 |L|;
  - *PT-only:* both have |FC| < 0.25 |L|;
  - otherwise *mixed*.
- **Readouts:** the call for level-change genes and the mean FC of level-change sets.
- **Caveat:** non-PT units sit inside the PT support and can carry PT spill-over. A PT-only change can therefore leak
  into them, so this filter is biased toward "kidney-wide".

### I4 · Induced and lost programs along the axis (unit level)

- **Programs:**
  - *Induced:* induced genes detected in ≥ 10% of units in both IR mice.
  - *Lost:* lost genes detected in ≥ 10% in both controls.
- **Per unit:** the share of the program's genes detected.
- **Per bin (≥ 20 units) and mouse:** the mean share.
- **Interior dip rule (induced program, per IR mouse):**
  - take the minimum over bins with centres 0.35–1.95;
  - it is a dip when the maxima on both sides exceed it by ≥ 0.1 × the mouse's range;
  - *replicated* when both IR mice have a dip and their minima lie ≤ 0.3 units apart.
- **Sensitivity:** drop units with `curve_end_extension` > 0.
- **Lost program:** its control profile shows where those genes live, and its IR profile shows what remains.
- **Descriptive:** detection profiles of the 12 induced genes with the highest IR detection.

### I5 · Sets: the discovery output

- **Library sets:** 71's AKI sets (`gene_sets_aki_vs_control.json`).
  - L and P use shared members, at least 10 per set.
  - A uses zonated members, at least 5 per set.
- **External programs** (own BH):
  - **Published injury:** 80's 24 programs (`programs_members.csv`).
  - **Time of day:** GSE277302 males. The ZT18-high and ZT6-high programs are genes with |log2 FC| ≥ 0.5 and Welch
    t-test p ≤ 0.01 (6 against 6), read as in notebook 44.
  - **Fasting:** GSE267280. The fasting-up and fasting-down programs are genes with |log2 FC| ≥ 1 at a mean CPM ≥ 10.
    There is one sample per condition, so these are point estimates.
  - **Sex:** Xiong 2023. Male-high and female-high programs are genes significant in any of the PT-S1/S2/S3 female
    vs male tables, by its `Direction`.
- **Set statistics:**
  - the mean L over members;
  - the location profile, the mean P(x) over members, with size M_set = its RMS and its peak at the argmax of |mean P|;
  - set A, the ratio of sums over zonated members.
- **Null.**
  - 2,000 random sets matched member-for-member by expression decile (`pp.strata_permutations`). The deciles are
    the controls' mean log2 CP10k among shared genes, or among zonated genes for A.
  - z = (observed − null mean) / null SD, and p is normal from z.
  - Tails: two-sided for L and A, upper for M_set.
  - BH correction runs within each statistic.
- **Calls:**
  - *Level set:* q ≤ 0.05; the four pairs' set L share its sign; both relabelings have |L_R| ≤ 0.5 |L|.
  - *Location set:* q ≤ 0.05; the two IR mice's set profiles correlate at ≥ 0.5; both relabelings have M_R ≤ 0.5 M.
  - *Gradient set:* q ≤ 0.05; both relabelings have set A ≥ 0.67. It is *erased* below the genome-wide survival and
    *kept* above it.
- **Redundancy:** sets with Jaccard ≥ 0.5 are grouped (`pp.redundancy_groups`); lists report the strongest set per
  group.
- **Signature agreement (descriptive):**
  - the Spearman correlation of gene L with the time-of-day, fasting and sex log2 FCs, for the condition and both
    relabelings;
  - the same correlation using the DCT and thin-limb FCs.

## 4 · Part II · Species (hvm export)

**Input.** 90's end contrasts (`d_end_contrasts.csv.gz`, sha256 `b4a5646f…effe1`). E is computed per specimen as the
mean of the last 8 bins minus the mean of the first 8, on that specimen's deviation profile. Expression strata use
90's `a_gene_classes.csv.gz` (sha256 `48d33943…f156b`) mouse log2 CP10k; genes without it form their own stratum.

### S1 · Genes

- **Graded** in a species: |E| ≥ 0.5 in both its specimens, with the same sign (90's rule).
- **Survival:**
  - for mouse-graded genes, A_h = mean E_human / mean E_mouse;
  - for human-graded genes, A_m = mean E_mouse / mean E_human.
  - Classes as in I2.

### S2 · Sets

- **Sets:** 71's hvm sets (`gene_sets_human_vs_mouse.json`), with ≥ 10 members graded in either species.
- **Statistics:**
  - the gradient cosine C_set = Σ E_m E_h / √(Σ E_m² Σ E_h²) over those members, on species means. It runs from −1
    (reversed) to +1 (conserved);
  - the amplitude ratio Q_set = log2(RMS E_h / RMS E_m).
- **Specimen checks:**
  - C for each of the four mouse × human specimen pairs;
  - the within-species cosines (A2 vs A4, COR1 vs MED1) over the same members.
- **Null:** 2,000 random sets of graded-in-either genes, matched by expression stratum; z, normal p two-sided, BH
  within statistic.
- **Calls:**
  - *Reversed program:* C_set < 0; q ≤ 0.05; all four pair cosines < 0; both within-species cosines ≥ 0.5.
  - *Conserved program:* C_set above its null; q ≤ 0.05; all four pair cosines > 0; both within-species cosines
    ≥ 0.5.
  - *Species-dominant program:* Q_set q ≤ 0.05, with the dominant species' within-species cosine ≥ 0.5. Q < 0 is
    mouse-dominant, Q > 0 human-dominant.

### S3 · Does the human gradient carry the mouse injury pattern? (repeats a scratch check)

- **Gene level:** E_human regressed on E_mouse and Δ_late, where Δ_late is the condition Δ averaged over bins ≥ 1.8.
  The same model is fitted with each relabeling's Δ_late in its place.
- **Set level:** the Spearman correlation across sets between C_set (species) and set A (injury).
- Descriptive.

## 5 · Reading rules and caveats

- Two mice per condition and one donor in two sections. Every q is a ranking device; nothing here is confirmatory.
- **Shared confounders.** IR and control may differ in harvest time, surgery and feeding. The relabelings cannot
  remove a difference shared by both IR mice; I3 and the external signatures are there to flag it.
- **Coordinate.** The package drops much of the injured S3 (`pt-package-issues.md`, item 7), so the late IR bins hold
  survivors. The coordinate's label genes respond to injury, so IR positions are not injury-invariant. I2 and I4 are
  measured on the coordinate. The anatomical check (glomerulus distance, depth) belongs to the next protocol.
- **Resolution.** Position differences below about 0.3 units are below the curve's resolution (85, C14).

## 6 · Recovery checks (named now; a miss is reported, not tuned)

1. **I2:** Slc7a13 is erased or reversed, and the genome-wide S2→S3 survival is below the S1→S2 survival.
2. **I4:** the induced program has a replicated interior dip in S1.
3. **I5 or I3:** the circadian sets (Reactome Circadian Clock, KEGG Circadian rhythm) or a time-of-day program are
   level sets, and the clock genes are kidney-wide.
4. **S2:** Hallmark hypoxia is reversed or mouse-dominant.

## 7 · Figures

- **91a · Injury decomposition:**
  - level against location size per gene, by call;
  - gradient survival, condition against the relabelings;
  - top library sets per component;
  - external programs (L, M and A z-scores);
  - PT level against the kidney-wide FC.
- **91b · Where injury acts:**
  - injury shape classes by location, with the relabeling nulls;
  - the induced and lost program profiles per mouse;
  - the location profiles of the top location sets.
- **91c · Species gradient programs:**
  - E mouse against E human;
  - set cosines against the null, with calls;
  - the top reversed, conserved and species-dominant sets;
  - the injury-link check.

## Notebook 91 results

`analysis/notebooks/91_pt_positional_decomposition.ipynb` (logic `91.positional_decomposition.1`) writes to
`results/pt_positional_decomposition/`. It followed the protocol frozen in 8b739ae. Guards G1–G3 passed. It runs in
about 0.5 min cold. The module is `pseudospace/positional_decomposition.py`, tested by
`tests/test_positional_decomposition.py` (8 tests).

**How the protocol was read (also in `run_record.json`).** Three post hoc additions, labelled:
1. **I1:** the level and location rules applied to each relabeling as the contrast, as a specificity count.
2. **I5:** the set z of L and M recomputed under each relabeling with the same nulls, to calibrate the random-set
   null.
3. **I4:** each induced gene's early-S1 against S3 detection, after figure 91b showed two patterns. Panel e now orders
   all 28 induced genes by this ratio, instead of the 12 most detected.

### I1 · Level and location (10,361 shared genes)

- **Level change:** 3,159 genes (2,897 level only, 262 with location too). Location change: 362 genes (100 location
  only).
- **Specificity:**
  - With R1 or R2 as the contrast, the same rules (without the pair and replication parts) give 40 and 3 level
    changes, against 3,159.
  - They give 0 and 1 location changes, against 362.
- **Injury shape** (ladder on the two IR mice's Δ). Non-flat genes number 6,737 for the condition, against 2,948 and
  2,768 under the relabelings.
- **Location-change genes are continuous:** 243 favour a continuous shape, 30 favour steps at the labels, and 89 are
  undecided.
  - Switches sit inside S2 (50), at S2/S3 (44) and inside S3 (31).
  - Bumps sit inside S3 (15) and inside S1 (12).

### I2 · Gradient survival (2,010 genes zonated in the controls)

| Pooled survival | IR vs control | R1 | R2 |
|---|---|---|---|
| whole PT | **0.52** | 0.92 | 1.15 |
| S1→S2 (1,739 genes) | 0.66 | 0.91 | 1.15 |
| S2→S3 (1,493 genes) | **0.31** | 1.14 | 1.23 |

- **Gene classes:**

  | Class | IR vs control | Relabelings |
  |---|---|---|
  | reversed | 64 | 6, 12 |
  | erased | 543 | 22, 18 |
  | attenuated | 729 | – |
  | kept | 658 | – |
  | amplified | 16 | – |
- **The S3 gradient is what injury removes.** Within S2→S3, under a third of the healthy contrast survives.
- **Male-biased and S3 genes fall by level, not by shape.** Slc7a13 (L −4.56), Slco1a1 (−3.64) and Kap (−2.03) keep
  their shape (A 1.0–1.1). Their loss is a drop along the whole PT, so they are "kept" here. Cyp7b1 and Slc22a7 fall
  below 2% detection in IR and are in the lost program (I4).

### I3 · Kidney-wide filter: uninformative

- DCT and thin-limb units in the PT export follow the PT changes almost one to one:
  - Spearman with L is 0.88 and 0.82, against 0.04–0.21 under the relabelings;
  - the slope is 1.07 and 1.27.
- 2,702 of 3,159 level-change genes are therefore "kidney-wide", and 10 are PT-only.
- Non-PT units inside the PT support carry PT spill-over (the protocol's caveat), so the filter cannot separate a
  whole-animal cause from PT injury here.
- **Next step:** distal tubules from the whole section (the package's stage-01 units and pass-1 labels), with a
  spill-over gauge (Havcr1 detection in distal units).

### I4 · Induced and lost programs (unit level)

- **Induced program:** 28 genes, off in both controls and on in ≥ 10% of both IR mice's PT units. The most detected
  are Havcr1, Fgb, Krt20, Lcn2, Uhrf1, Mcm5, Smpdl3b, Timp1, Clcf1 and Sox9.
- **The program is U-shaped along the PT:**
  - early S1 0.27–0.35;
  - minimum at 0.85 units, 0.08–0.12;
  - end of S3 0.46–0.52.
- **The dip replicates** in both IR mice at the same bin (excess 0.50–0.57 of the range), and holds without units
  beyond the curve's start (0.40–0.43).
- **Post hoc zones:**
  - 22 genes are S3-high: the proliferation genes (Mki67, Top2a, Mcm5, Uhrf1, Ccne2, Rad51) plus Sox9, Sprr1a,
    Lrp8 and Cbr3.
  - 5 genes are along the whole PT with the dip: Havcr1, Lcn2, Krt20, Timp1, Crlf1.
  - 1 gene is early-S1-high: Fgb.
- **Lost program:** 12 genes (Cyp7b1, Ppp1r1b, Cyp2d9, Abca13, Serpina1d, Slc22a19, …). In controls they live in late
  S2 and S3 (detected share rising from 0.05 to 0.47); in IR they are near 0.

### I5 · Sets

**Library sets** (1,449 sets; 261 non-redundant level groups, 121 location groups):
- **Level up:** 507 sets.
  - Top: E2F targets, Myc targets, EMT, DNA replication, mitotic anaphase, S phase, TNF-α/NF-κB, nuclear envelope
    reassembly, programmed cell death.
- **Level down:** 31 sets.
  - Top: KRAS signalling (down), branched-chain amino-acid catabolism, valine/leucine/isoleucine degradation,
    tryptophan catabolism, butanoate, bile acid metabolism, pantothenate and CoA, peroxisome, SLC transport.
- **Location:** 236 sets. 114 of 121 non-redundant groups peak in S3.
  - Gains peaking in S3: cell cycle, TNF-α, p53, apoptosis, inflammatory response, hypoxia.
  - Losses peaking in S3: PPAR signalling, peroxisome, peroxisomal lipid metabolism, cholesterol metabolism and
    biosynthesis, plasma lipoprotein assembly and remodelling, fatty-acid metabolism.
  - Injury takes lipid handling and peroxisomal oxidation away specifically in S3, beyond its PT-wide level.
- **Gradient survival:** 6 erased (mitotic and E2F sets, whose control gradients are small and whose relabeled ratios
  are unstable), 0 kept. Set survival tracks the genome-wide 0.52 (sets: median 0.45, SD 0.18). The erasure is
  genome-wide rather than pathway-specific.
- **Calibration (post hoc):**

  | | IR vs control | R1 | R2 |
  |---|---|---|---|
  | Median set z L | 1.37 | 0.44 | 0.08 |
  | Sets with q L ≤ 0.05 | 554 | 205 | 175 |
  | Sets with q M ≤ 0.05 | 466 | 33 | 51 |
  | Full level rule | 538 | 15 | 1 |
  | Full location rule | 236 | 0 | 0 |

  The random-set null alone is anticonservative for L, because it ignores correlation between genes in a set. The
  per-set relabeling clause carries the specificity, so every q in I5 is a ranking device.

**External programs (own BH):**
- **Sex:** the male-biased PT program (Xiong 2023, 496 genes) falls by level (L −0.44, z −12.5), the strongest
  external signal. The female-biased program rises slightly (L +0.14) and is not a level set.
- **Fasting:** fasting-up genes rise (L +0.37, z 5.8; also a location set peaking at 2.05) and fasting-down genes fall
  (L −0.82, z −7.0). The IR mice carry a fasting signature.
- **Time of day:**
  - The ZT18-high (122 genes) and ZT6-high (85) programs do not change level (L 0.00, −0.03), and gene L correlates
    with the time-of-day log2 FC at ρ 0.07 (relabelings 0.01, −0.06).
  - The core clock genes do fall: Ciart −2.12, Nr1d1 −1.12, Dbp −1.10, Tef −1.03, Hlf −0.98; Cry1 rises +0.59.
  - This argues against a large harvest-time difference between IR and control mice. It cannot exclude a small phase
    offset, which would show first in the high-amplitude clock genes, or injury damping the clock.
- **Published injury programs:** 27 of 30 are level sets in the published direction; 16 are also location sets,
  peaking at 2.25–2.35.
  - Kirita's injured S1/2 and S3, the core IRI set and Gerhardt's failed-repair and VCAM1 programs all rise most in
    S3.
  - Gerhardt 2023's S1-up list is the one miss (L −0.10).

### Part II · Species

- **Gene survival:**

  | Class | Mouse-graded genes, in human (1,342) | Human-graded genes, in mouse (627) |
  |---|---|---|
  | reversed | 267 | 160 |
  | erased | 640 | 143 |
  | attenuated | 261 | 102 |
  | kept | 154 | 141 |
  | amplified | 20 | 81 |

  "Reversed" here is A < −0.25 on species means, looser than 90's rule (both specimens graded).
- **Sets: no calls.** Of 381 sets, none is reversed, conserved or species-dominant at q ≤ 0.05.
  - The set cosines spread around the genome-wide 0.23 (z mean −0.20, SD 0.93).
  - The most negative z are oxidative phosphorylation, xenobiotic metabolism, fatty acids, PPAR and tryptophan
    (z −2.2 to −2.9, q > 0.4).
  - Species gradient differences are gene-level, not pathway-coherent beyond random genes.
- **Injury link (S3):** none. β for Δ_late is 0.012 (relabelings 0.069, −0.057), and the residual ρ is 0.055
  (0.038, −0.035), as the scratch check found. Across sets, species cosine and injury survival correlate at ρ 0.20
  (380 sets, no null): sets with robust gradients score higher on both.

### Recovery checks (section 6)

1. **Partial.** S2→S3 survival 0.31 is below S1→S2 0.66. Slc7a13 is "kept" (A 1.02): its loss is a PT-wide level drop.
2. **Recovered.** The induced program has a replicated interior dip at 0.85 in both IR mice.
3. **Not recovered.**
   - The circadian sets are not level sets (Reactome L −0.05, KEGG −0.16), because their members move in opposite
     directions.
   - The time-of-day programs do not change.
   - The clock genes are "kidney-wide", but I3 is uninformative.
4. **Not recovered.** Hallmark hypoxia's species cosine 0.21 equals the genome-wide value.

### What the screen gives the paper (descriptive)

1. **A filter with a stated null.** Every gene and pathway gets a level, a location (with a continuous or stepped
   shape) and a gradient survival. Relabeled specimens bound each call: 40 and 3 level genes and 0 and 1 location
   genes under relabeling, against 3,159 and 362.
2. **Injury erases the S3 gradient.** The IR PT keeps 0.52 of the healthy gradient overall and 0.31 from S2 to S3,
   against 0.9–1.2 under relabeling.
3. **Where injury acts is continuous and S3-centred.**
   - Location-change genes favour continuous shapes 243 to 30.
   - 114 of 121 location groups peak in S3: proliferation, TNF-α and p53 up; PPAR, peroxisome and lipoprotein
     handling down.
4. **Late S1 is the least injured stretch.** Injury switches on genes from zero along the whole PT, with a replicated
   minimum at 0.85 units. Proliferation is concentrated in S3.
5. **The external filters name non-injury contributors.** The male-biased program and a fasting signature shift in IR;
   the time-of-day transcriptome does not. These go beside every AKI claim.
6. **Null results:** no pathway-level species program beyond random genes; no mouse-injury pattern in the human
   gradient; the kidney-wide filter needs whole-section distal units.
