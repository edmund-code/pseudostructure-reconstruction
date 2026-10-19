# Is the loss of S3 identity after IRI the collapse of the androgen program? (notebook 111)

**The question.** After IRI, our mice lose the late-PT gradient: S2→S3 survives at 0.31 against 0.66 for S1→S2
(notebook 91), and this replicates in Gharaie's spatial data (notebook 105). S3 marker genes such as Cyp7b1 and
Slc22a7 fall below 2% detection, and the male-biased PT program falls (notebooks 91 and 94). Several of the lost S3
genes are male-specific (Cyp7b1 is about 1,000-fold male-biased in S3; Xiong 2023). All our mice are male.

Is the "S3 identity loss" partly the collapse of the androgen receptor (AR)-dependent male program, rather than only
an injury response at S3?

**What the literature says.** The scout read these in full unless marked; the full review is in
`docs/research/androgen-s3-literature.md`.
- **AR sets the male S2/S3 program directly** (Xiong 2023).
  - After castration, 85.5% of core male-biased genes fall.
  - Nephron-specific Ar knockout lowers 82.5% of them.
  - AR knockout perturbs nothing in S1.
  - In the authors' per-segment knockout tables, AR-maintained genes number S1 18, S2 200, S3 306 (the scout's
    count).
- **Surgery alone crashes testosterone** (Soljancic 2013, rat).
  - Testosterone falls 90% by 3 h after sham surgery and after I/R alike.
  - It recovers from 12 h in sham but not after I/R.
  - Castration makes IRI worse in that study.
- **Females lose S3 identity genes on the same 4 h–2 d schedule** (Dixon 2022, female mice). That set is not the AR
  program (1 of 30 genes AR-maintained). So a sex-independent identity loss exists.
- **AR is required for part of the injury response** (Allison 2025, preprint): PT-specific Ar deletion abolishes
  Sox9 induction and reduces injury.
- **The same lost genes have been described without a sex reading.**
  - Kirita 2020's SCENIC regulon activity has Ar as an S3 identity regulon that goes negative in every injured state
    (Dataset S3, not discussed in the paper).
  - Gerhardt 2021 found the female S3 gene Slc7a12 ectopically on in male IRI (RNAscope) and did not mention
    hormones.

**What is known before freezing** (from the scout's work and earlier notebooks; this is not new evidence here):
- notebook 94's male-biased program (K7) in Kirita is already −4.9 at 4 h and −5.9 to −6.2 at 12 h–14 d;
- the scout reports Cyp7b1 at −4.9 and Slc7a13 at −6.5 within Kirita PTS3 nuclei at 4 h;
- notebook 91's check H11d found Slc7a12 not gained in our IR mice (L −0.73).

No AR knockout set has been joined to our IR effects.

Frozen on 2026-10-18.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/111_pt_androgen_s3.ipynb`; logic `111.androgen_s3.1`; seed 111.
- Results `results/pt_androgen_s3/`.
- **AR instrument:** `data/external/xiong2023_ar_ko/`, male nephron-lineage Ar knockout against wild type, per PT
  segment.
  - **AR-maintained in segment s:** log2FC ≤ −1 with q < 0.05 (lost in the knockout).
  - **AR-repressed:** log2FC ≥ +1 with q < 0.05 (gained in the knockout).
  - The S3 sets are primary; S2 is secondary.
- **Sex instrument (descriptive):** `data/external/xiong2023_sex_bias/` (female against male wild type).
- **Our IR effects:** notebook 91's `i_genes` table:
  - L, the level change: the mean of the four IR–control specimen pairs;
  - its per-pair columns, and the relabeling levels `L R1` and `L R2`;
  - `Δ late`, `log2 CP10k (controls)`, `level change`.
- **Control gradient:** E_ctrl per gene from the AKI export controls, the same pooled estimator as notebook 104:
  log2 CP10k of units with x ≥ 0.7 minus x ≤ 0.3, the mean of Ctrl1A2 and Ctrl1A4.
- **Transport membership:** notebook 76's T_PT.
- **External, Kirita 2020:** notebook 94's `pseudobulk_counts.csv.gz` (per mouse × PT type × time).
- **External, GSE144622** (`data/external/vinas2020_gse144622/`): sorted PT cells at 24 h; IR against sham in males
  and in females.

## 2 · Tests (fixed now)

**Universe:** notebook 91's genes that the knockout table tested.

**Matched null:** random genes matched on the joint decile of E_ctrl and the decile of control expression, 2,000
draws (`matched_means`).

**P1 · Feminisation: do AR-repressed genes rise while AR-maintained genes fall?**
- Statistic: the mean L of the S3 AR-repressed set and of the S3 AR-maintained set, each against its matched null.
- **De-repression holds** when the AR-repressed mean L is above the null's 95th percentile, and above the same set's
  mean `L R1` and `L R2` (the relabelings).
- **AR loss holds** when the AR-maintained mean L is below the null's 5th percentile and below both relabelings.
- **Feminisation holds** when both hold. If both sets fall together, the change is not feminisation.
- `Δ late` of each set is reported (S3-weighting).

**P2 · Does AR dependence explain the IR effect beyond the gradient and transport?**
- Model per gene: L ~ E_ctrl + control log2 CP10k + T_PT membership + the S3 knockout log2FC (continuous), by least
  squares.
- Null: 2,000 permutations of the knockout log2FC within the joint deciles.
- **Holds** when the knockout coefficient is > 0 (AR-maintained genes fall more) with one-sided p ≤ 0.05, and exceeds
  the same coefficient fitted to `L R1` and to `L R2`.

**P3 · Timing in Kirita (within identity-retaining PTS3 nuclei).**
- **Mice:** those with ≥ 30 PTS3 nuclei.
- **Per time point:** the mean log2 CPM difference (IR mouse against the mean of control mice) over the S3
  AR-maintained set, and over a **matched identity set**.
  - The identity set is genes S3-high in Kirita controls (PTS3 against PTS1 log2 CPM ≥ 1).
  - It excludes genes that are AR-dependent (|log2FC| ≥ 0.5 in any knockout segment) or sex-biased (q < 0.05 in any
    segment).
  - It is matched to the AR set on the control PTS3 expression decile.
- **AR gap** = AR-maintained change − identity change.
- **Androgen withdrawal predicts** a gap near 0 at 4 h, more negative at 12 h–2 d, and persisting at 14 d.
- **Holds** when the gap at 12 h or 2 d is more negative than at 4 h by ≥ 0.5 log2, and the 14 d gap is ≤ −0.5.
- Expected outcome (from the scout): the AR set is already down at 4 h, so this counts against withdrawal.

**P4 · Sex specificity in GSE144622 (difference in differences).**
- For the S3 AR-maintained set and a control set matched on sham PT expression (baseMean decile), restricted to the
  genes in both tables:
  - DiD_sex = [mean lfc(AR) − mean lfc(control)] in males − the same in females.
- Null: 2,000 matched control draws.
- **Male-specific collapse holds** when the male AR-minus-control difference is below the null's 5th percentile and
  the female one is not.

**P5 · Steroid-receptor confound (descriptive).** The mean L of a glucocorticoid-response set (Fkbp5, Tsc22d3, Per1,
Sgk1, Zbtb16, Dusp1, Klf9, Errfi1). AR, GR and MR share the response element, and surgery raises glucocorticoids.
This list is canonical and was chosen without a literature check; it is descriptive only.

## 3 · Outputs and figures

- **Tables:** `sets.csv`, `p1.csv`, `p2.csv`, `p3_kirita.csv`, `p4_gse144622.csv`, `p5_gr.csv`,
  `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 111a: mean L of AR-maintained and AR-repressed sets against their nulls and the relabelings;
  - 111b: per-gene L against the S3 knockout log2FC, coloured by E_ctrl;
  - 111c: Kirita time course of the AR gap;
  - 111d: GSE144622 male against female changes.

## 4 · Caveats

- The knockout is developmental (Six2-Cre), not acute withdrawal. AR-dependent genes may differ from genes that fall
  within 2 days of androgen loss.
- The knockout depletes Ar more in S3 than in S1 (log2FC −8.6 against −5.2), which could inflate S3 gene counts.
- GSE144622 female PT responds much less at 24 h (497 against 4,101 genes), which is why P4 uses a difference in
  differences against matched controls.
- Our IR effect is a level change over the whole PT. Position is reported, not modelled.
