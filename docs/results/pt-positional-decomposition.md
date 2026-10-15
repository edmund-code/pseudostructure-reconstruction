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
