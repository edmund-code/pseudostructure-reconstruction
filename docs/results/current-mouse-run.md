# Results — mouse nephron pseudospace, healthy vs AKI

Summary of the figures and tables produced by `6_mouse_only_pseudospace_final.ipynb`,
run 2026-08-13. Written as a handoff for interpretation.

Everything below is **description of what the outputs show**, with the numbers behind each
claim. Where I flag a concern it is marked as such and separated from the observation. Nothing
here is a re-analysis — it is a reading of the files in `results/mouse_only_v5/`.

---

## 0. What produced these results

| | |
|---|---|
| Samples | `Ctrl1A2`, `Ctrl1A4` (Control) · `IR2A2`, `IR2A4` (AKI / ischemia-reperfusion) |
| Unit of observation | segmented tubule cross-sections from Visium HD, **not single cells** |
| Biological replicate | specimen — **n = 2 vs 2** |
| Tubules after pass-1 QC | 43,848 |
| Tubules in the DPT trajectory | 39,574 |
| PT cohort used for healthy-vs-AKI | 6,628 healthy · 10,812 AKI |
| Genes tested | 9,657 |
| Clustering | Leiden, **resolution 0.7 → 13 clusters** |
| Integration | Harmony, θ = 6, on `X_harmony` |
| Pseudospace | Scanpy DPT, PT-anchored root, oriented by an early→late marker axis |

### Configuration caveats that affect interpretation

These are real and worth resolving before drawing firm conclusions:

1. **The 13 cluster labels are the pipeline's auto-generated argmax suggestion, taken verbatim.**
   The notebook prints a suggested `COARSE_LABELS` dict from per-cluster marker-panel scores,
   with the instruction to *verify against the dotplot*. The dict in the notebook is
   character-identical to that suggestion. Two assignments look wrong (§1.2).
2. **The historical notebook had its `COARSE_LABELS_FINGERPRINT` guard commented out.** The
   active workflow now enforces the reference fingerprint and rejects stale labels before DPT.
   Any new clustering must be reviewed and relabelled rather than reusing these numeric IDs.
3. **Permutation p-values have a hard floor of 2/6 ≈ 0.333.** With n = 2 vs 2 there are only
   C(4,2) = 6 sample-level relabelings, and they form three mirror-image (healthy/AKI swap)
   pairs whose statistic is identical. An exact inclusive-tail p therefore stops at 2/6, and
   **no result here can be significant at conventional thresholds.** (The earlier 1/6 floor came
   from counting strict-greater outcomes only, which under-counts ties.) Every ranking below is
   by *effect size*. This is a design limit, not a defect — the work is descriptive.

---

## 1. Cell typing and trajectory structure

### 1.1 Integration — `celltyping/integration_qc_silhouette.csv`

| embedding | batch ASW (lower better) | bio ASW (higher better) |
|---|---|---|
| `X_pca` (uncorrected) | 0.0106 | 0.1508 |
| `X_harmony` (θ=6) | **−0.0212** | **0.1823** |

Harmony moved batch ASW below zero (batches well mixed) *while raising* bio ASW. That is the
right direction on both axes simultaneously — the usual worry with an aggressive θ is that batch
mixing improves by destroying biology, and that did not happen here.

### 1.2 Cluster identity — `celltyping/coarse_marker_dotplot.png`, `coarse_cluster_top_markers.csv`

Most clusters have clean, unambiguous marker support:

| cluster | top markers | label | reading |
|---|---|---|---|
| 0 | `Alpl, Slc5a2, Gatm, Slc4a4, Slc34a1` | PT-S1 | ✅ `Slc5a2` is S1-specific |
| 1 | `Cyp4b1, Slc22a6, Acy3, Slc13a3, Slc27a2` | PT-S2 | ✅ |
| 2 | `Slc6a18, Kap, Napsa, Mep1a, Slc27a2` | PT-S3 | ✅ |
| 3 | `Tmem52b, Klk1, Tmem213, Clcnkb, Slc8a1` | DCT2 | ⚠️ mixes DCT (`Slc12a3`,`Trpv5`) with intercalated-cell (`Tmem213`) |
| 4 | `Hsd11b2, Aqp2, Ly6e, Tmsb4x, Aqp3` | CCD | ✅ principal cell |
| 5 | `Klk1, Umod, Car15, Clcnkb, Slc12a1` | cTAL | ✅ |
| 6 | `Sema3g, Podxl, Synpo, Nphs2, Plat` | Podocyte | ✅ (also carries endothelial markers — glomerular tuft) |
| **7** | **`Krt19, Tmsb4x, Sprr1a, Ahnak, Krt18`** | **ATL** | ❌ see below |
| 8 | `Tmsb4x, Plet1, Aqp2, Akr1b3, Bcam` | IMCD | ✅ |
| 9 | `Cryab, Vim, Tnc, Thbs1, S100a6` | SmoothMuscle | ⚠️ this is stroma/fibroblast, not smooth muscle |
| **10** | **`mt-Co2, mt-Atp6, Umod, Mt1, Nccrp1`** | **mTAL** | ⚠️ mitochondria-dominated |
| 11 | `Acta2, Tpm2, Tagln, Myh11, Cald1` | SmoothMuscle | ✅ genuine |
| 12 | `Umod, Slc5a3, Atp1a1, Defb1, Ppp1r1a` | mTAL | ✅ |

Two of these matter:

- **Cluster 7 (250 tubules) labeled `ATL`.** Its markers are keratins and `Sprr1a` — a
  urothelial or injured / de-differentiated epithelial signature, not thin limb. Because `ATL`
  rolls up into the `AL` family, these tubules are **kept in the trajectory** rather than
  excluded. Three independent outputs say this label is wrong (§1.3, §1.4).
- **Cluster 10 (2,426 tubules) labeled `mTAL`.** Dominated by mitochondrial transcripts and
  `Mt1` — the classic stressed / low-quality signature. Also kept in the trajectory.

Cluster 9's mislabel is harmless in effect: stroma and smooth muscle are both removed as
non-tubule either way.

**No cluster was assigned to the descending thin limb (DTL).** The `DTL1/2/3` marker panels are
empty across all 13 clusters, so the reconstructed nephron has no descending limb.

### 1.3 Trajectory topology — `celltyping/paga_segment_topology.png`, `paga_connectivity_matrix.csv`

The PT chain is strong and anatomically correct:
`PT-S1 → PT-S2` = **0.715**, `PT-S2 → PT-S3` = **0.568**, `mTAL → cTAL` = **0.571**,
`CCD → IMCD` = **0.851**.

The `ATL` node is not where anatomy puts it. Its strongest connections are:

| ATL connects to | connectivity |
|---|---|
| **CCD** | **0.801** |
| **IMCD** | **0.635** |
| mTAL (expected neighbour) | 0.214 |
| PT-S3 (expected neighbour) | 0.173 |

Transcriptionally this cluster sits with the collecting duct, ~4× more strongly than with the
segments a thin limb should border.

### 1.4 Pseudospace ordering — `celltyping/dpt_by_segment.png`, `dpt_by_segment.csv`

Median DPT by segment, in the order the trajectory actually places them:

| segment | n | median DPT | anatomically expected? |
|---|---|---|---|
| PT-S1 | 6,530 | 0.019 | ✅ |
| PT-S2 | 5,753 | 0.203 | ✅ |
| PT-S3 | 5,157 | 0.427 | ✅ |
| cTAL | 2,406 | 0.474 | ✅ |
| mTAL | 11,038 | 0.497 | ✅ |
| CCD | 1,883 | 0.725 | ⚠️ before DCT2 |
| DCT2 | 3,597 | 0.879 | ⚠️ after CCD |
| IMCD | 2,960 | 0.912 | ✅ |
| **ATL** | **250** | **1.000** | ❌ should be ≈0.45, between PT-S3 and TAL |

**The PT arm is textbook** — a clean, tight, monotonic S1 → S2 → S3 progression with narrow
interquartile ranges. That is the strongest single piece of evidence that pseudospace is
recovering real anatomy.

Two ordering problems: `ATL` is pinned at the terminal end of the trajectory (consistent with
§1.2 and §1.3 — it is not a thin limb), and `DCT2` sits *after* `CCD` when anatomy runs
DCT → CNT → CCD.

### 1.5 Family sub-trajectories — `heatmaps/subset_recomputed_dpt_diagnostics.csv`

Spearman correlation between each family's recomputed sub-trajectory and its marker axis:

| family | n | ρ (DPT vs marker axis) |
|---|---|---|
| PT | 17,440 | **0.797** |
| DCT | 3,597 | 0.721 |
| CNT_CD | 4,843 | 0.689 |
| **AL** | 13,694 | **0.225** |

Three of four families reconstruct well. **`AL` is the outlier** — and `AL` is exactly the family
that absorbed both questionable clusters (the keratin-high cluster 7 and the mito-dominated
cluster 10, together 2,676 tubules).

### 1.6 Independent fidelity check — `concordance/three_axis_concordance.csv`

Agreement among three axes that are derived independently: DPT, an early→late marker score, and
physical distance to the nearest glomerulus.

| sample | n | DPT vs dist-to-glom | marker vs dist-to-glom | DPT vs marker |
|---|---|---|---|---|
| Ctrl1A2 | 8,782 | 0.465 | 0.729 | 0.695 |
| Ctrl1A4 | 8,022 | 0.419 | 0.667 | 0.732 |
| IR2A2 | 11,316 | 0.382 | 0.522 | 0.780 |
| IR2A4 | 11,454 | **0.201** | **0.237** | 0.847 |

DPT and the marker axis agree well and consistently (ρ 0.70–0.85). Agreement with *physical*
distance-to-glomerulus is weaker and **degrades monotonically with injury severity** —
controls ≈ 0.42–0.47, IR2A4 ≈ 0.20.

Note the pattern: DPT-vs-marker agreement is *highest* in IR2A4 (0.847) precisely where both
agree least with physical position. Two readings, and the data here cannot separate them:
(a) AKI genuinely distorts the tissue so molecular position decouples from physical position, or
(b) both molecular axes are being driven by a shared injury program in the most injured animal.

---

## 2. Marker cascades — `heatmaps/`

### `mouse_control_total_marker_heatmap.png` vs `mouse_aki_total_marker_heatmap.png`

Both conditions show an ordered marker cascade along pseudospace:
`Slc5a2` (S1, peaks ≈0) → `Slc22a6` (S2, ≈0.15) → `Slc7a13` (S3, ≈0.30) → `Cldn10`/`Cldn16`
(TAL, ≈0.50) → `Aqp2`/`Aqp3`/`Slc14a2` (CD, ≈0.70) → `Pvalb`/`Trpv5` (DCT, ≈0.90) →
`Slc4a1`/`Slc26a4` (intercalated, ≈0.95).

**The cascade survives injury** — the AKI panel preserves the same left-to-right ordering. That
is a substantive result in itself: the segment-identity axis is not destroyed by IR.

Two features visible in both panels reflect §1.4: `Sptssb` (ATL) peaks at the extreme right edge
rather than mid-trajectory, and the DCT markers (`Pvalb`, `Trpv5`) peak *after* the CD markers.

Condition differences visible by eye: the AKI PT block is compressed and weaker, and `Cldn10`
shows an anomalous band at the very start of the AKI axis that has no counterpart in control.

---

## 3. Healthy vs AKI over the PT cohort — `healthy_vs_aki/`

The analysis decomposes each gene's pseudospace trajectory into a **level** shift (vertical
offset) and a **shape** change (altered spatial profile), ranked by `shape_rms`.

Distribution of `shape_rms` across 9,657 genes: median 0.027, 75th pct 0.056, max 1.140 — so the
top genes sit ~20–40× above the typical gene. The signal is concentrated in a small set.

### 3.1 Top shape-changed genes — `top_shape_trajectories_zscore_genes.png`

Top 20 by `shape_rms`: `Inmt, Cyp7b1, Cyp4b1, Cyp2e1, Krt20, Mep1b, Pah, Ugt2b38, Slc22a7,
Slc34a1, Krt8, Odc1, Crot, Hsd11b1, Serpina1f, Akr1c14, Dnase1, S100a6, Havcr1, Pck1`.

Of the top 30, **23 are down in AKI and 7 are up** (`Krt20, Krt8, Odc1, S100a6, Havcr1, Anxa3,
Anxa2`).

The trajectory figure is the clearest result in the whole run, and it shows **two distinct
phenomena, not one**:

1. **Loss of the proximal-tubule differentiation gradient.** `Cyp7b1`, `Mep1b`, `Slc22a7` are
   flat and low across all of pseudospace in healthy tissue until ≈0.3, then rise steeply into
   S2/S3. In AKI that rise **does not occur at all** — the curves stay flat. This is not a
   uniform downshift; the *spatial program itself* is absent. `Inmt`, `Cyp2e1`, `Cyp4b1`,
   `Ugt2b38` show the mirror image: a healthy mid-axis peak that is flattened in AKI.
2. **Gain of a spatially-structured de-differentiation program.** `Krt20` and `Krt8` are flat and
   low across healthy pseudospace and rise steeply in AKI toward late PT (0.3 → 0.5). The injury
   response is not uniform along the tubule — it increases distally.

`Havcr1` (Kim-1), the canonical AKI biomarker, appearing in the top 20 and up in AKI is a useful
external sanity check that the contrast is capturing real injury biology.

Within each condition, the per-specimen dashed curves track the pooled curve tightly — the two
animals within a condition agree closely, and the healthy/AKI gap is much larger than the
within-condition spread.

### 3.2 Effect size vs replicate noise — `specimen_separation_genes.csv`

The honest n=2v2 framing: condition gap relative to between-specimen SD.

| gene | condition gap | between-specimen SD | separation t | clears t₍.₉₇₅,df2₎ = 4.303 |
|---|---|---|---|---|
| Inmt | 1.793 | 0.088 | **48.3** | ✅ |
| Mep1b | 1.197 | 0.079 | 16.7 | ✅ |
| Cyp2e1 | 1.080 | 0.135 | 12.4 | ✅ |
| Cyp7b1 | 1.106 | 0.105 | 12.0 | ✅ |
| Cyp4b1 | 0.960 | 0.267 | 9.6 | ✅ |
| Krt20 | 1.372 | 0.187 | 7.0 | ✅ |

All 12 profiled genes clear the reference. The condition difference is 7–48× the variation
between animals of the same condition.

### 3.3 Robustness — `loso_shape_stability_genes.png`

Leave-one-specimen-out: for every top gene, all four LOO estimates bracket the all-specimen
value within roughly ±0.15 `shape_rms`. **No gene collapses when any single animal is dropped**,
so no result is carried by one specimen. With n=2 per arm this is the strongest robustness
statement available.

### 3.4 Is the coordinate itself condition-dependent? — `pseudospace_structure_diagnostics.csv`

The concern: if AKI shifts where tubules sit along pseudospace, a "shape change" could be an
artifact of resampling a fixed curve.

| metric | value |
|---|---|
| Wasserstein distance, healthy vs AKI DPT | 0.0099 |
| KS statistic | 0.0329 (p = 2.7 × 10⁻⁴) |
| marker-axis monotonicity, healthy | 0.840 |
| marker-axis monotonicity, AKI | 0.755 |

The distributions are **statistically distinguishable but numerically almost identical** — a
Wasserstein distance of 0.01 on a unit interval, and a KS statistic of 0.03, are tiny. The p-value
is small only because n ≈ 17,000. The coordinate is close to condition-invariant, and monotonicity
degrades only modestly under injury (0.84 → 0.76).

### 3.5 The circularity check — `physical_axis_sensitivity.png`, `injury_shape_adjudication.csv`

This is the most important control, and the most nuanced result. It re-runs the shape analysis on
a **physical, injury-independent coordinate** (glomerular-density depth proxy).

- Coordinate agreement, DPT vs physical: Spearman **+0.33** (weak)
- Shape-effect agreement: Spearman **+0.70** (strong)

The metabolic genes — `Slc34a1`, `Pah`, `Cyp2e1`, `Cyp4b1`, `Cyp7b1`, `Inmt`, `Ugt2b38`,
`Mep1b`, `Odc1`, `Slc22a7` — sit high on **both** axes. Their shape change is not an artifact of
the marker-derived coordinate.

**`Krt20` and `Krt8` do not replicate.** They have high `shape_rms` on DPT (0.72–0.91) but near-zero
on the physical axis (0.09–0.15). The keratin/de-differentiation signal is at least partly a
property of the DPT coordinate itself — which makes sense, since a de-differentiation program will
influence a marker-derived axis.

`injury_shape_adjudication.csv` flags `Cyp2e1` as *"injury-associated BUT axis-basis
(de-diff confounds coordinate)"*, while shortlisting `Krt20`, `Krt8`, `Cyp4b1`, `Odc1`, `Pah` as
coordinate-independent by a different test (within-pseudospace-bin correlation, sign consistency
0.93–1.00 across 60 bins). **The two tests disagree about the keratins** and that disagreement is
worth resolving.

### 3.6 Injury burden vs position — `injury_confounding_diagnostics.csv`

| metric | value |
|---|---|
| acute injury score, AKI − healthy (Cohen's d) | **1.88** |
| failed-repair PTC score, AKI − healthy (Cohen's d) | **1.07** |
| acute injury score, ρ with pseudospace | 0.169 |
| failed-repair score, ρ with pseudospace | 0.345 |
| acute injury, top-decile median pseudospace | 0.381 (vs 0.191 overall) |
| failed-repair, top-decile median pseudospace | 0.440 (vs 0.191 overall) |

Injury scores separate the conditions strongly (d ≈ 1.1–1.9) but correlate only weakly with
position (ρ ≈ 0.17–0.35). However, the **most injured decile sits ~2× further along pseudospace**
than the cohort median (0.38–0.44 vs 0.19) — injury is concentrated in later PT (S2/S3), which is
consistent both with S3's known ischemia susceptibility and with the distal-increasing keratin
pattern in §3.1.

### 3.7 Pathways — `pathway_level_shape_results.csv`

1,397 pathways tested. Top by `shape_rms`, and the theme is consistent:

| pathway | library | shape_rms | level |
|---|---|---|---|
| Primary bile acid biosynthesis | KEGG | 0.293 | −0.53 |
| Metabolic Disorders Of Biological Oxidation Enzymes | Reactome | 0.248 | −0.13 |
| Synthesis Of Bile Acids/Salts via 7α-hydroxycholesterol | Reactome | 0.216 | −0.45 |
| Mitochondrial Fatty Acid β-Oxidation (saturated) | Reactome | 0.190 | −0.47 |
| Nitrogen metabolism | KEGG | 0.188 | −0.11 |
| MyD88 Deficiency (TLR2/4) | Reactome | 0.185 | **+0.33** |
| Butanoate metabolism | KEGG | 0.175 | −0.74 |
| IRAK4 Deficiency (TLR2/4) | Reactome | 0.166 | **+0.32** |

Oxidative/catabolic metabolism (bile acid synthesis, fatty-acid β-oxidation, biological
oxidation) is **down** with altered spatial shape; innate-immune/TLR signalling modules are
**up**. This is the pathway-level restatement of the gene-level story.

Caveat carried from the run log: only 48% of `KEGG_2019_Mouse` symbols matched the tested PT gene
background, attributed to PT detection filtering.

---

## 4. What the results support

1. **Pseudospace reconstructs the proximal tubule convincingly.** Monotonic S1→S2→S3 DPT ordering
   with tight IQRs, PAGA connectivity 0.72 / 0.57 along the chain, subset-DPT ρ = 0.80 to the
   marker axis, and an ordered marker cascade in both conditions.
2. **The reconstruction is weaker outside the PT** — `AL` sub-trajectory ρ = 0.23, a DCT/CCD
   ordering inversion, no DTL recovered at all, and one cluster placed at entirely the wrong end.
3. **AKI removes the PT spatial differentiation program rather than just lowering expression.**
   Genes like `Cyp7b1`, `Mep1b`, `Slc22a7` lose their S2/S3 induction entirely. Effect sizes are
   7–48× between-animal variation and survive leave-one-specimen-out.
4. **A spatially-graded de-differentiation program appears in AKI**, increasing toward distal PT
   (`Krt20`, `Krt8`), alongside `Havcr1` induction.
5. **The metabolic finding is not an artifact of the coordinate** — it replicates on an
   independent physical axis (shape-effect ρ = 0.70). **The keratin finding is coordinate-dependent
   and needs adjudication.**
6. **Nothing here reaches conventional significance and nothing can**, given the 2/6 permutation
   floor at n = 2 vs 2. This is an effect-size-ranked, descriptive result.

## 5. Open questions for interpretation

1. **Cluster 7 (`Krt19`/`Sprr1a`/`Krt18`, 250 tubules).** Injured/de-differentiated PT, urothelial
   contamination, or something else? It is currently labeled `ATL`, sits at DPT = 1.0, and
   connects to CCD at 0.80. If it is injured PT it belongs in the PT cohort and would change §3
   materially. If it is urothelium it should be dropped. Either way the `ATL` label looks wrong.
2. **Cluster 10 (mitochondria-dominated, 2,426 tubules).** Real mTAL with high mitochondrial
   content — mTAL is genuinely mitochondria-rich — or a low-quality/stressed artifact? It is 18%
   of the `AL` family and a plausible driver of that family's poor ρ = 0.23.
3. **Do the keratins survive?** §3.5's two tests disagree. Resolving this decides whether
   de-differentiation is a genuine spatial finding or a coordinate artifact.
4. **Why does DCT2 order after CCD?** Real biology of the IR cohort, a labeling problem in
   cluster 3 (which mixes DCT and intercalated-cell markers), or a DPT branch-collapse issue?
5. **The IR2A4 concordance drop** (physical-axis ρ = 0.20 vs ~0.45 in controls). Tissue distortion
   from injury, or both molecular axes tracking a shared injury program?
6. **Absent DTL.** Is the descending thin limb genuinely unrecoverable at this resolution, or is
   it hidden inside another cluster? It leaves a real gap in the PT→TAL transition, where PAGA
   connectivity is weakest (PT-S3→ATL = 0.17).

## 6. Figure index

| file | shows |
|---|---|
| `celltyping/coarse_marker_dotplot.png` | 13 clusters × marker panels — the basis for every label |
| `celltyping/paga_segment_topology.png` | trajectory topology; ATL's misplacement |
| `celltyping/dpt_by_segment.png` | pseudospace distribution per segment — clearest anatomy check |
| `celltyping/labelled_umap.png`, `pass2_umap.png`, `total_dpt_umap.png` | embeddings |
| `heatmaps/mouse_{control,aki}_total_marker_heatmap.png` | full marker cascade per condition |
| `heatmaps/mouse_*_pt_*_heatmap.png` | PT-specific cascades, global and subset DPT |
| `healthy_vs_aki/top_shape_trajectories_zscore_genes.png` | **the core result** |
| `healthy_vs_aki/loso_shape_stability_genes.png` | leave-one-specimen-out robustness |
| `healthy_vs_aki/physical_axis_sensitivity.png` | **the key circularity control** |
| `healthy_vs_aki/pseudospace_structure_diagnostics.png` | is the coordinate condition-dependent |
| `healthy_vs_aki/injury_matched_diagnostics.png` | position-matched injury contrast |
| `healthy_vs_aki/volcano_shape_{genes,pathways}.png` | effect size vs permutation p (floor 2/6) |
| `concordance/three_axis_concordance.png` | DPT vs marker vs physical agreement per specimen |
