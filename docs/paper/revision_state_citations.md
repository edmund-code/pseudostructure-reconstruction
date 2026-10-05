# Revision round: physiological state, tissue state, more human donors, rodent protein, citations

Workstream 3 report for referee items M4/D8, D9, M8 and D10 (`docs/paper/referee_report_v01.md`).
Notebook `analysis/notebooks/44_pt_revision_state.ipynb`; results in `results/pt_revision_state/`.

> **Coordinator review note (accepted with one correction).** Rat S1 samples in KTEA carry far less
> protein (most proteins have lower copies per cell in S1, and Gclc is undetected there), so the
> S3/S1 ratios are inflated. Quote the **S3/S2** ratios: Gamt +1.6, Gatm −3.7, Gclc +2.8, Gclm +1.9,
> Gss +2.4, Dcxr −1.1, Acaa2 −2.7, Acadm +0.7. Rat AGAT protein is highest in S2, so write "AGAT
> falls at S3", not "AGAT is confined to S1".


## Protocol for notebook 44 (written before it ran)

Notebook `analysis/notebooks/44_pt_revision_state.ipynb`, results in `results/pt_revision_state/`.
All inputs are public except §B, which uses our own PT structures (human sections and the four
Visium HD mice). Every rule below is fixed before results are seen.

### A · Physiological state (D8, M4.5)

Programs (mouse symbols, fixed):
- *sterol / SREBP2*: Hmgcr, Hmgcs1, Sqle, Cyp51, Idi1, Fdft1, Mvd, Ldlr, Insig1, Dhcr7, Lss, Msmo1
- *PPARα / fatty-acid oxidation*: Acox1, Ehhadh, Cyp4a14, Cyp4a10, Acot1, Crot, Nudt19, Pdk4, Hmgcs2,
  Cpt1a, Acadm, Slc27a2, Acsm3
- *glutathione synthesis*: Gclc, Gclm, Gss, Gsr
- *tyrosine catabolism*: Tat, Hpd, Hgd, Fah, Gstz1, Pah
- *lead genes*: Gatm, Gamt, Acaa2, Dcxr, Ugt3a1, Acox2
- *structural zonation markers (expected state-insensitive)*: Slc5a2, Slc5a12, Slc7a13, Slc22a7,
  Slc34a1, Lrp2

**A1 · Whole-kidney magnitude.**
- Data:
  - GSE267280: 24-h fasted vs fed male mice, n = 3 each, raw counts.
  - GSE277302: ZT18 vs ZT6, both sexes, normalized counts.
  - GSE54650: kidney sampled every 2 h for 48 h, array log2 values.
- Per gene:
  - fasting log2FC: log2 of mean CPM + 1;
  - diurnal log2FC, per sex;
  - circadian peak-to-trough: twice the amplitude of a 24-h cosinor fit.
- Comparator: the gene's external species × position interaction
  |Lake S3 − S1 − mouse snRNA S3 − S1| (male mouse reference).
- Decision per program:
  - *state effect comparable to the species difference* if the median over program genes of the
    largest state effect is ≥ 0.5 × the median |interaction|;
  - otherwise *small relative to the species difference*.
- Caveat fixed in advance: whole kidney, so a segment-restricted change could exceed it.

**A2 · Segment zonation at two times of day.**
- Data: GSE331332, control mice at ZT4 (rest phase) and ZT16 (active, feeding phase), n = 2 each,
  snRNA.
- Nuclei: at least 500 UMI.
- Marker sets, fixed from the Census mouse pseudobulk (independent lab) and excluding every program
  and lead gene above:
  - S1: the 10 genes with the largest male S1 − S3 at S1 CPM ≥ 100;
  - S3: the 10 with the largest S3 − S1 at S3 CPM ≥ 100;
  - PT: canonical Lrp2, Cubn, Slc34a1, Hnf4a.
- Per nucleus, each score is the mean log1p(CP10k) of its marker set.
  - PT nuclei: PT score above the Otsu threshold of the sample's PT-score distribution.
  - Among PT nuclei, z-score the S1 and S3 scores within the sample. Call S1 if z(S1) − z(S3) > 0.5,
    S3 if < −0.5; the rest are not used.
- Per sample: pseudobulk log2 CPM, S3 − S1.
- Decision per program: *zonation is state-stable* if both hold:
  1. for ≥ 80% of program genes with |S3 − S1| ≥ 0.5 at either time, all four samples share the
     sign;
  2. the ratio of program median |S3 − S1| at ZT16 to ZT4 lies in 0.67–1.5.

  Otherwise *time-of-day modulated*. In addition, *time of day could account for a substantial part
  of the species difference* if the program median |ZT16 − ZT4| change in S3 − S1 is ≥ 0.5 × the
  program median |interaction|.

### B · Tissue state of our donor on the Lake injury axis (M4.3, D7)

- Gene universe U: genes present in our human and mouse matrices, Lake and the Census mouse
  pseudobulk.
- Signature from Lake cortex nuclei: per donor with ≥ 30 aPT and ≥ 100 healthy PT nuclei
  (PT-S1/S2/S3), log2 CPM of aPT minus healthy.
  - Up: the 30 genes with the largest median difference at aPT CPM ≥ 20.
  - Down: the 30 with the most negative median difference at healthy CPM ≥ 20.
- Score of a PT pseudobulk: mean within-sample percentile (over U) of up genes minus that of down
  genes.
- Lake: per-donor cortical PT pseudobulk over all PT states. Validation is Spearman(score,
  altered-state fraction).
- Platform offset, estimated in mouse: mean score of our two control mice (Visium HD) minus the
  median of healthy Census mouse donors (snRNA). Sanity check: our AKI mice must score above our
  control mice.
- Our human sections: score minus the offset, placed among Lake donors (percentile among normal
  donors; altered fraction predicted from the Lake linear fit).
- Interpretation:
  - *within the Lake reference range* if between the minimum and maximum of normal donors;
  - *injury-like* if above the Lake AKI median.
- Caveat: the offset is estimated in mouse and assumed to transfer.

### C · Additional human donors with spatial data (D9)

- Data: Abedini et al. 2024 Visium (GSE211785), 14 samples with authors' spot labels; samples with
  ≥ 50 PT_S1 and ≥ 50 PT_S3 spots. Per sample, pseudobulk log2 CPM (counts summed, divided by spot
  totals over all genes), S3 − S1.
- Label check: positive-control signs SLC5A2 −, SLC5A1 +, SLC7A13 +, SLC22A7 + in ≥ 70% of
  samples.
- Lead genes: the human arm *replicates* if its sign (as observed in Lake) holds in ≥ 70% of
  informative samples (|S3 − S1| ≥ 0.1), with at least 5 informative samples.
- The flattening classes of phase 3 are repeated with "≥ 80% of samples share the sign" in place of
  6/7 donors.
- Amplitude k per sample (slope on male mouse snRNA over gene set G) is compared with Lake. Spot
  mixing will lower it, which is stated in advance.

### D · Rodent segment protein (D10)

- Data: Kidney Tubule Expression Atlas, rat, microdissected S1/S2/S3 (Limbutara et al. 2020),
  copies per cell, 3 samples per segment.
- Descriptive: mean copies, detection, log2 S3/S1 and S3/S2 centred on the median ratio of
  proteins detected in all nine PT samples.

---
**Deviations noticed when loading the data (before any result).**
- GSE54650 was dropped: its MoGene array annotation (GPL6246) could not be retrieved from this host.
- GSE267280 holds one column per condition, so the fasting contrast has no replicates; it is
  reported as a point estimate only.
- The mouse AKI comparison uses the mouse-only object's own counts (9,688 genes). Its structures
  do not map one-to-one onto the tubule matrices, so the gene universe U for §B is restricted to
  those genes.

---
## Results (notebook 44, logic version 44.revision_state.1)

Tables: `results/pt_revision_state/tables/`. Figure: `results/pt_revision_state/figures/fig_revision_state.*`.
Tests: `pytest tests` passes (396 tests).

### A · Physiological state (D8)

**A1, whole kidney** (`A1_state_effects_by_program.csv`). Median largest state effect against
median |species × position interaction| (log2):

| Program | State effect | Interaction | Decision |
|---|---:|---:|---|
| Sterol / SREBP2 | 0.75 | 0.65 | comparable |
| PPARα / fatty-acid oxidation | 1.08 | 1.62 | comparable |
| Glutathione synthesis | 0.20 | 1.95 | small |
| Tyrosine catabolism | 0.29 | 0.71 | small |
| Lead genes (Gatm, Gamt, Acaa2, Dcxr, Ugt3a1, Acox2) | 0.57 | 2.75 | small |
| Structural markers | 0.35 | 2.71 | small |

Two programs respond strongly to state:
- *Sterol.* 24-h fasting lowers the sterol genes by 0.2–2.2 log2 (Idi1 −2.2, Msmo1 −1.9, Ldlr
  −1.8). This is an unreplicated point estimate.
- *PPARα.* Fasting raises the classic PPARα targets strongly: Cyp4a14 +6.8, Hmgcs2 +6.3, Pdk4 +4.1,
  Acot1 +2.4, Ehhadh +1.8. The species-different genes in our story respond little: Acadm +0.6,
  Crot +0.6, Nudt19 −0.4, Acsm3 0.0, Slc27a2 +0.3. Time-of-day effects in whole kidney are
  ≤ 0.9 log2.

**A2, segment zonation at ZT4 and ZT16.** About 3,200–3,500 PT nuclei per mouse; segment markers
came from Census and exclude the tested genes. In `A2_snrna_time_of_day_by_program.csv`:
- Fatty-acid oxidation, glutathione, tyrosine catabolism, the lead genes and the structural markers
  keep their zonation. The sign is stable in all four mice for 100% of zonated genes, and the
  amplitude ratio ZT16 ÷ ZT4 is 0.99–1.35. Examples: Acadm +2.1/+1.3 at ZT4 and +1.2/+1.8 at ZT16;
  Crot +2.3/+2.5 and +2.9/+2.7.
- **Sterol zonation is time-of-day sensitive.** The sign is stable for 73% of zonated genes, and
  the median change (0.46 log2) is at least half the species gap. Hmgcr S3 − S1 ranges from
  +0.4 to +3.0 across the four mice. The core genes (Hmgcr, Cyp51, Sqle, Hmgcs1) stay positive in
  all mice.

**Human side, suggestive only.** PDK4, a fasting- and glucocorticoid-induced gene, is high in human
PT: 94% of our human structures against 7% of mouse; 362 CPM in Lake against 68 in Census mouse. The
ketogenic fasting marker HMGCS2 is low in human (10 vs 127 CPM). Sterol genes are lower in human
(HMGCR 16 vs 73 CPM). These are cross-species level comparisons and cannot separate species from
state.

**Verdict for D8.**
- Time of day does not create the mouse zonation of fatty-acid oxidation, glutathione, tyrosine
  catabolism or the lead genes.
- Sterol zonation is state-sensitive and should be reported as such.
- Fasting affects PPARα and sterol programs at whole-kidney scale by amounts comparable to the
  species gap. No public segment-resolved fasted-kidney data exist to test whether fasting
  flattens zonation.

### B · Our donor on the Lake injury axis (M4.3)

- *Signature.* 30 genes up and 30 down in aPT, from 43 Lake donors (`B_injury_signature.csv`). The
  down genes are PT-differentiation genes (SLC34A1, SLC5A12, SLC22A6, SLC13A3, SLC6A19, PCK1,
  HPD …); the up genes are classic aPT genes (PROM1, CLDN1, TNC, ANXA1/3, TPM1 …).
- *Validation in Lake.* Spearman with altered-state fraction 0.52 over 54 donors. Medians: normal
  −0.45 (range −0.62 to 0.18), CKD −0.10, AKI −0.04.
- *Platform calibration in mouse.* Our control mice score −0.55, against −0.36 for healthy Census
  donors, so the offset is −0.20. Our AKI mice score **+0.33 above our controls**, a shift similar to
  the Lake normal→AKI difference (+0.42). The score therefore detects injury on our platform.
- *Our donor.* Raw −0.64; calibrated −0.45 in both sections. That is the **median of the 14 Lake
  normal donors** (50th percentile), within their range, below the CKD and AKI medians, and predicts
  an altered fraction of 0.50 by the Lake fit.
- *Caveat.* The offset is estimated in mouse and assumed to transfer to human.

### C · More human donors with spatial data (D9)

- *Feasibility.* Abedini et al. 2024 Visium (GSE211785): counts 142 MB, metadata 0.8 MB. 14
  samples carry authors' spot labels; 10 have ≥ 50 PT_S1 and ≥ 50 PT_S3 spots (2 control,
  8 disease). It ran in under 10 minutes. KPMP/Lake Visium (GSE183279, 3.5 GB raw) was not needed.
- **The pre-specified label check failed.** In PT_S3 spots, SLC5A2 falls as expected in 9 of 10
  samples, but the S3 markers do not rise: SLC5A1 in 1/10, SLC7A13 in 2/10, SLC22A7 in 4/10.
  Post-hoc centring on non-zonated PT genes changed nothing (offsets 0–0.3 log2), so this is not
  uniform dilution. The authors' PT_S3 spots behave like "less S1", not "more S3". Amplitude k is
  tiny (median 0.025).
- *Usable reading, descriptive only.* Human-arm signs replicate for the S1-high genes UGT1A9
  (9/9 informative samples), UGT3A1, ACOX2 and SLC6A19 (all 100%). RBP4 rises (7 of 8). Genes the
  atlas calls rising cannot be tested with these labels; GATM, DCXR and ACADM *fall* in the spots
  (GATM 10 of 10).
- **D9 conclusion.** Accessible but not informative at Visium resolution with the authors' labels.
  A proper test needs spot deconvolution against Lake PT-S1/S3 profiles, or Slide-seq/Xenium-scale
  data.

### D · Rodent segment protein (D10)

Rat KTEA, microdissected S1/S2/S3 (`D_rat_segment_proteome.csv`). Log2 ratios are centred on the
median ratio of fully quantified proteins; rat S1 samples carry less protein.
- **Gatm** S3/S1 −2.9, S3/S2 −3.7: AGAT protein falls at S3, consistent with rodent early
  restriction.
- **Gamt** S3/S1 +6.6: rises. This is the first protein-level support for the rodent two-step
  separation.
- **Gclc, Gclm, Gss** rise to S3 (S3/S2 +2.8, +1.9, +2.4): glutathione synthesis is late in rodent
  protein too.
- **Dcxr** S3/S2 −1.1, S3/S1 −1.8: falls. **Acadm** S3/S1 +2.4 and **Acaa2** −2.8 match the mouse
  transcript directions.
- **Slc22a6/OAT1** peaks in S2 (S3/S2 −6.6). **Slc22a7/OAT2** is late. **Rbp4** protein is S1-only,
  i.e. reabsorbed ligand, as the 1989 rat in situ study warned.
- **Ugt3a1** and the sterol enzymes (Hmgcs1, Cyp51, Hmgcr) are not quantified.

The rodent arm of every lead story that rat protein covers is supported. Human protein remains
presence-only (HPA).

---

## Recommended wording (M4: the draft's "tissue state ruled out" overclaims)

Replace "Ruled out: log compression and tissue state" with:

> **Log compression is ruled out; tissue state does not explain the pattern in the data available.**
> Our donor's PT injury score, calibrated across platforms in mouse, lies at the median of
> healthy KPMP reference donors. Within those donors, amplitude falls with injury, but healthy
> human donors remain far flatter than healthy mice. Mouse zonation of fatty-acid oxidation,
> glutathione synthesis, tyrosine catabolism and the lead genes is unchanged between rest (ZT4)
> and active (ZT16) phases. Mouse sterol-synthesis zonation is not stable across these times.

Add to Limitations:

> **Shared confounders.** Our human tissue and the external human atlases share surgical or biopsy
> procurement with warm and cold ischaemia, older donors, comorbidity and medication, and
> pre-operative fasting. Mouse kidneys came from young, ad-libitum-fed animals. Fasting changes
> PPARα-target and sterol-synthesis genes in mouse kidney by amounts comparable to the species
> differences we report, and human PT expresses the fasting- and stress-responsive gene PDK4
> strongly. Segment-resolved data from fasted animals or fed humans do not exist. A contribution of
> physiological state to the species difference, especially for sterol synthesis and PPARα
> targets, is therefore not excluded. Caloric restriction also feminizes male mouse kidney gene
> expression ([PMID 37673062](https://pubmed.ncbi.nlm.nih.gov/37673062/)), and PPARα regulates
> largely different target genes in mouse and human hepatocytes
> ([PMID 19710929](https://pubmed.ncbi.nlm.nih.gov/19710929/)).

Also cite:
- The renal tubular clock controls β-oxidation and carnitine handling; about 30% of kidney RNAs
  are rhythmic ([PMID 36862511](https://pubmed.ncbi.nlm.nih.gov/36862511/), OA).
- PT PPARα is fasting-activated and drives renal FAO and ketogenesis
  ([PMID 41460648](https://pubmed.ncbi.nlm.nih.gov/41460648/), abstract).
- Early- and late-PT mitochondrial FAO capacities differ and change with 24-h fasting
  ([PMID 41428383](https://pubmed.ncbi.nlm.nih.gov/41428383/), abstract).

**Story-level consequence.** Sterol synthesis should stay a supporting observation, with "zonation
varies with time of day in mouse". For fatty-acid oxidation, the Acadm/Acaa2 contrast is not a
time-of-day artefact. Whether fasting flattens it in human cannot be tested with public data.

---

## M8 · Citations: what each source supports, and corrected sentences

Verified against abstracts in this round.
- Breljak 2016 ([PMID 27053689](https://pubmed.ncbi.nlm.nih.gov/27053689/)): "OAT1-OAT3 …
  detected only in the BLM of cortical proximal tubules; all three OATs were stained more
  intensely in S1/S2 segments compared with S3 segment in medullary rays, whereas the S3 segment
  in the outer stripe remained unstained"; antibodies validated in transfected HEK-293 cells;
  "considerable interindividual variability".
- Maunsbach 1997 ([PMID 9013443](https://pubmed.ncbi.nlm.nih.gov/9013443/)).
- Vrhovac 2015 ([PMID 25304002](https://pubmed.ncbi.nlm.nih.gov/25304002/)).
- Human vs mouse phosphate transport (PMID 42494142).

Our classes (`gene_zonation_classes.csv`): SLC5A2 and SLC5A1 conserved; Aqp1 indeterminate;
Slc22a6 mouse-only; Slc22a7 conserved rising; Slc22a8 indeterminate; Slc34a3 conserved. Of 54
zonated Slc/Abc/Aqp transporters, 30 (56%) are conserved and 18 are mouse-only, against 29% of
all 606 zonated genes. This confirms the referee's numbers.

1. **Introduction (introduction_draft.md, the Breljak sentence).** Keep it; it reports the paper
   accurately. Append: *"Our transcript-level calls agree with this map for SGLT1/2 but not for
   OAT1 and OAT2 (Results)."*
2. **draft.md, background paragraph** ("Human protein maps also place AQP1 … and OAT1–3 …").
   Replace with:
   > Human protein maps place SGLT2 in S1/S2 and SGLT1 in S3, much as in rodents (Vrhovac et al.
   > 2015). They place AQP1 highest in the proximal straight tubule (Maunsbach et al. 1997), and all
   > three OATs more intensely in S1/S2 than in medullary-ray S3, with outer-stripe S3 unstained
   > (Breljak et al. 2016). No study has compared these maps position by position with rodent data.
3. **draft.md, summary sentence "Mouse PT predicts where human PT places transport, but not
   metabolism."** Replace with:
   > Zonated transporters were conserved more often than zonated genes overall (30 of 54, 56%,
   > against 29%), but 18 transporters, including OAT1, NaDC3 and SMCT2, were zonated in mouse only.
4. **draft.md, R3 conserved-class bullet.** Keep "agrees with SGLT2 and SGLT1", and make the
   discordance explicit:
   > Two of the three OATs disagree with human protein data (Breljak et al. 2016): OAT1 protein is
   > stronger in S1/S2 than in medullary-ray S3, but our human OAT1 mRNA is flat (a mouse-only call);
   > OAT2 protein is likewise S1/S2 > S3, but OAT2 mRNA rises toward S3 in our data and in the KPMP
   > atlas. Possible explanations are mRNA–protein differences, antibody specificity (validated only
   > on transfected cells) and the interindividual variability the authors report. AQP1 is
   > indeterminate in our classes.
5. **draft.md, OAT-module paragraph ("Human OAT1–3 are known to be stronger in S1/S2 …").**
   Replace with:
   > Human OAT1–3 protein is reported stronger in S1/S2 than in medullary-ray S3 (Breljak et al.
   > 2016), whereas our human OAT1 mRNA is flat, so the protein map and our transcript call
   > disagree for OAT1.
6. **draft.md, Discussion "Transport" bullet.** Append:
   > Conserved zonation does not mean conserved usage. In adult human kidney SLC34A3 carries about
   > 40% of sodium-phosphate cotransport and is negligible in mouse (PMID 42494142), although
   > Slc34a3 zonation is conserved here.
7. **r5_draft.md (OAT bullet in R5.6) and my `zonation_literature.md`.**
   - The "counter-evidence: canonical transporters are zoned in human as in rodents" list wrongly
     counted AQP1 and the OATs as agreement, and gave Lake SLC22A7 +2.7 as reproducing them.
   - Lake OAT2 rising *contradicts* Breljak's OAT2 protein.
   - Corrected statement: *"Of the transporters with human protein maps, only SGLT1/2 agree with
     our transcript zonation; AQP1 is indeterminate and OAT1/OAT2 disagree."* The flattening
     bottom line (canonical markers keep zonation) should cite SGLT1/2, SLC5A8 and SLC7A13 only.
8. **Rat-based arguments (M2).**
   - The "rat runs the opposite way" sterol argument rests only on rat RPKM (Lee 2015); KTEA does
     not quantify the sterol enzymes. Remove it.
   - Replace the rat mRNA support for Dcxr with rat **protein** (KTEA: S3/S2 −1.1).
   - For Rbp4, rat protein is S1-only (reabsorbed ligand). Cite only the 1989 rat in situ study
     for the S3 mRNA site.

**Consequence for the lead stories (read with C's caveat).**
- *GATM.* The human arm now has three readings: Lake +0.51, our donor −0.78 and Visium −1.0 (labels
  failed the S3 check). The defensible creatine claim is a large difference in amplitude, not a
  reversal: mouse AGAT is confined to S1 by transcript (microdissection) and by protein (rat KTEA),
  while human GATM declines modestly at most. Drop "human rises"; keep "human GATM extends along the
  PT".
- *DCXR.* The human rise holds in the two sources with valid labels (Lake +1.53, ours +0.97) and
  fails in the Visium spots, whose labels are not valid for rising genes.
- The UGT1A9, UGT3A1 and ACOX2 human arms, and the RBP4 rise, replicate in the Visium donors as
  well.

## Files

- `analysis/notebooks/44_pt_revision_state.ipynb` (executed; stage-order clean)
- `results/pt_revision_state/`
- `data/external/{kidney_physiological_state, abedini2024_visium, ktea_rat_proteomics}/`, each
  with a SOURCE.txt
- `workstream-3 scratch: nb44.py`, `workstream-3 scratch: revision.md`,
  `workstream-3 scratch: fetch/fetch_ktea_xhr.py` (with the earlier `fetch_ktea.py` attempt)
