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
