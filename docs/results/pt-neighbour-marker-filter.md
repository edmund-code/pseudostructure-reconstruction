# Which gene-level findings are neighbour-structure markers? A cell-type specificity filter (notebook 109)

**The question.** Notebook 104 rebuilt every PT unit from bins at least 6 µm inside its polygon. Gradients and AKI calls
survived. Notebook 108 then found that distal-disease genes graded inside PT units also survive the core rebuild:
- SLC12A1 (thick limb) is late in both species;
- SLC12A3 (DCT) is early in mouse;
- AQP2, KCNJ1 and SCNN1A are late in human.

These genes follow the PT's neighbours. So core survival is necessary but not sufficient evidence of PT origin. This
notebook adds the missing filter: a gene is a **neighbour marker** when a single-nucleus atlas places its expression
mainly outside the PT. It then screens every gene list the paper relies on.

**Precedent** (methods scout, full text): smiDE (Vasconcelos et al. 2026) excludes genes whose expression in
neighbouring other-type cells exceeds that in the analysed type (overlap ratio > 1). Its MERFISH case removed 31.5% of
naive neighbour-contrast hits that way. Our rule is the atlas version of that ratio.

Frozen on 2026-10-18, before any atlas profile was computed.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/109_pt_neighbour_marker_filter.ipynb`; logic `109.neighbour_filter.1`; seed 109.
- Results `results/pt_neighbour_marker_filter/`.
- New helper `h5_group_sums` (notebook 98's streaming sum, moved to `pseudospace/external_replication.py` with a test).
- **Mouse atlas:** Kirita 2020 control nuclei summed per broad type, read from notebook 105's cached reference (its
  control columns). PT = PTS1 + PTS2 + PTS3.
- **Human atlas:** KPMP Lake 2023 raw counts streamed from `raw/X`.
  - Nuclei from Healthy_reference_tissue, Healthy_stone_donor and Tumor_nephrectomy, all regions.
  - Summed per `subclass.l1`: PT, TAL, DTL, ATL, DCT, CNT, PC, IC, POD, PEC, EC, FIB, IMM, VSM/P. PapE and NEU are
    dropped (< 400 nuclei).

## 2 · The rule (fixed now)

- CP10k is taken per type on the atlas genes.
- A gene is a **neighbour marker** in a species when the maximum CP10k over non-PT types is ≥ 1 **and** ≥ 4 × its PT
  CP10k. In mouse, NewPT is excluded from both sides.
- Human genes use the human atlas, mouse genes the mouse atlas. A gene absent from the atlas is "not assessable".

## 3 · Lists screened and readings

| List | Source | Species |
|---|---|---|
| Graded early and late genes | notebook 90 (E ≤ −0.5 or ≥ 0.5 in both specimens) | each |
| Conserved, reversed, mouse-only and human-only genes | notebook 90 | both atlases |
| AKI level and location calls | notebook 91 | mouse |
| Conserved handoff genes (11 pairs) | notebook 101 | both |
| Named literature positions (H4) | notebook 108 | both |
| Programs INJ, TPT, CYC | notebook 96 | mouse |
| Driver genes of each method's calls | notebook 93 (positional F-test, clustering, row N; AKI and human vs mouse) | per contrast |

**R1 · How much of each list is neighbour marker?**
- Reported per list against an expression-matched null (2,000 draws from the universe, by expression decile).
- A list is **clean** when ≤ 5% of its assessable genes are neighbour markers.

**R2 · Do the two contamination checks agree?**
- The AUC of notebook 104's edge enrichment ε (mean over a species' specimens, human vs mouse run) for neighbour
  markers against other genes, per species.
- **Agree** when AUC ≥ 0.6.

**R3 · Method exposure (specificity version of yardstick 8).** For each method's driver genes:
- the share that are neighbour markers;
- the null;
- the ratio.

The ranking is compared with notebook 104's edge exposure.

**R4 · Flag list.** Every flagged gene in the lists above is written out with its top atlas type. The paper must not
rest a claim on a flagged gene without saying so.

## 4 · Outputs and figures

- **Tables:**
  - `atlas_specificity.csv.gz` (gene, species, PT CP10k, top non-PT type and CP10k, ratio, neighbour marker);
  - `lists.csv`, `flags.csv`, `r2_agreement.csv`, `r3_exposure.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 109a: share of neighbour markers per list against its null;
  - 109b: edge enrichment of neighbour markers against other genes.

## 5 · Caveats

- Nuclear RNA: a gene whose transcripts leave the nucleus fast can look less PT-specific than it is.
- A gene expressed in both PT and a neighbour at similar levels is not flagged. Its gradient can still carry some
  neighbour RNA.
- The 4× rule is fixed, not tuned. Some genuine PT genes with higher expression elsewhere (for example in thick limb)
  will be flagged, and that is the conservative direction.

## Notebook 109 results

`analysis/notebooks/109_pt_neighbour_marker_filter.ipynb` (logic `109.neighbour_filter.1`) writes to
`results/pt_neighbour_marker_filter/`. It followed the protocol frozen in c55224c and ran in 8 s. Streaming the KPMP
atlas took 15 s.

**Atlases and guards.**
- G1: the mouse atlas has 13 broad Kirita types.
- G2: the human atlas holds 104,490 KPMP reference nuclei over 14 subclasses.
- By the frozen rule, 2,848 of 27,133 mouse atlas genes and 2,186 of 35,432 human ones are neighbour markers.
  SLC12A1, SLC12A3 and UMOD are flagged in both species; LRP2 and SLC34A1 are not.

**How the protocol was read.**
- Mouse lists use the Kirita atlas and human lists use KPMP.
- Expression deciles come from notebook 62's inventories.
- R2 uses notebook 104's human vs mouse run.
- A first scratch run cached the atlas type labels as an object array that could not be reloaded. The cache was
  rebuilt with string labels; the numbers are unchanged.
- **Post hoc (section 5.1):**
  - R1 flagged 7–28% of most lists, close to the expression-matched null.
  - So a stricter "suspect" class is added: a neighbour marker that is also edge-enriched in notebook 104 (ε ≥ 1 in
    both specimens of the species, or in all four AKI mice).

**R1 · Neighbour markers by the frozen atlas rule** (share flagged, null mean in brackets; clean when ≤ 5%):

| List | Mouse atlas | Human atlas |
|---|---|---|
| Graded early | 0.28 (0.19) | 0.16 (0.10) |
| Graded late | 0.14 (0.19) | 0.19 (0.11) |
| Conserved | 0.26 (0.19) | 0.18 (0.10) |
| Reversed | 0.28 (0.19) | 0.16 (0.10) |
| Mouse only / human only | 0.16 / 0.26 | 0.07 / 0.14 |
| Conserved handoff genes | **0.00: clean** | **0.00: clean** |
| Named literature positions (H4) | 0.07 | **0.00: clean** |
| AKI level / location calls | 0.20 / 0.21 (0.20 / 0.19) | – |
| Programs INJ / TPT / CYC | 0.30 / **0.02: clean** / 0.11 | – |

The atlas rule alone is broad. A gene expressed four times more in any other cell type is flagged even when it has a
real PT gradient, so most lists sit near the null. Graded genes, especially in human, carry 1.5–1.8 times more
neighbour markers than expression-matched genes.

**R2 · The two contamination checks agree.** Notebook 104's edge enrichment is higher for neighbour markers:
- AUC 0.65 in mouse (median ε 0.25 against 0.04);
- AUC 0.78 in human (0.52 against 0.02).

**R3 · Exposure of each method's driver genes to neighbour markers:** 1.10–1.23 times the null. The positional
F-test is highest (1.22 and 1.23), against 1.10–1.17 for clustering and 1.14–1.16 for position-split ORA. This ranking
is the reverse of notebook 104's edge exposure, so the two yardsticks disagree on method order. Both put the
differences between methods within about 1.1–1.9 times.

**Post hoc · Suspect genes** (neighbour marker and edge-enriched):

| List | Suspect | Genes |
|---|---|---|
| Graded early (mouse) | 8 / 729 (1.1%) | Egfl7, Ifi44, Svep1, Exoc3l2, Dlc1, Itga9, Igfbp3, Synpo |
| Graded late (mouse) | 1 / 583 (0.2%) | Rerg |
| Graded early (human) | 13 / 376 (3.5%) | FN1, CALB1, FLT1, MGP, RHCG, HSD11B2, SLC8A1, ACTA2 and others |
| **Graded late (human)** | **20 / 250 (8.0%)** | **SLC12A1, UMOD, EGF, KCNJ1, AQP2, WFDC2, STK39** and others |
| Conserved (human atlas) | 11 / 189 (5.8%) | SLC12A1, EGF, UMOD, CALB1, RHCG, HSD11B2, MGP and others |
| Reversed | 0 / 72 (mouse); 2 / 74 (human: NAV2, SFRP1) | |
| Conserved handoffs, H4 positions, INJ, TPT | 0 | |
| AKI level calls | 37 / 3,078 (1.2%) | collagens, Fbn1, Sparc, Kdr, H2-Aa, H2-Eb1 (interstitium, endothelium, MHC II) |
| AKI location calls | 1 / 353 | Igfbp3 |

**What this answers.**
1. **Neighbour-structure RNA reaches a small, identifiable set of gradient genes.**
   - It is most visible at the human S3 end: thick-limb and collecting-duct genes (SLC12A1, UMOD, EGF, KCNJ1, AQP2)
     are 8% of human late genes.
   - It also reaches a few early genes: DCT/CNT genes (CALB1, SLC8A1) and vascular genes (FLT1, MGP, ACTA2).
   - The "conserved" list holds 11 such genes, because the anatomy that places the thick limb next to S3 is conserved
     too. These are conserved neighbours, not conserved PT gradients.
2. **The findings the paper relies on are clean:**
   - the 11 conserved transporter handoffs;
   - the literature positions;
   - the transport program;
   - the injury program;
   - the reversed gradients in mouse (0 of 72);
   - AKI location calls (1 of 353).
   - The AKI level calls carry 37 interstitial and immune genes (1.2%), the expected infiltrate and matrix signal.
3. **As a method:**
   - Atlas specificity and bin-level edge enrichment agree (AUC 0.65 and 0.78), and each alone is too broad or too
     narrow.
   - The combination flags 0–8% of a list. It should be applied before any gene list is published.
   - The flagged genes are in `flags.csv` and `posthoc_suspect.csv`.
