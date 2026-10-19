# Public spatial datasets of mouse AKI for a position-along-PT replication (scout, 2026-10-18)

**Purpose.** We want an external spatial dataset of mouse ischaemia–reperfusion injury in which PT units can be
placed along S1→S3. With it we can test whether injury effects depend on position (notebook 105). Ideal: early time
points (4 h–3 d), at least 2 animals per condition, processed counts with coordinates.

**How it was read.**
- **Full text read** (methods, design passages, figure legends, data statements):
  - Cheung 2022, Gharaie 2023, Xuanyuan 2025, Dixon 2022 (PMC author manuscript);
  - Melo Ferreira 2021, Polonsky 2024, Wu 2024, Janosevic 2021;
  - Marshall 2022, Wang 2026, Singh 2026, Wijaya 2025.
- **Abstract only:** Rietjens 2026 (bioRxiv refused) and a Kidney Week 2026 abstract (Matsuura/Rabb).
- File sizes come from GEO `filelist.txt` and from the repository APIs.

## Summary

| Dataset (accession) | Technology | Model, sex | Time points | Animals per condition | PT segments | Usable without FASTQ |
|---|---|---|---|---|---|---|
| **Cheung 2022, JCI Insight (GSE200115)** | Visium v1, 55 µm, whole kidney | bilateral IRI 19 min, male B6, 12–16 wk | d0, 12 h, d1, d6, d28 | **3** (1 section each) | S3 zone by Slc7a13/Napsa/Nudt19; no labels deposited | counts (h5); coordinates from the v1 barcode grid |
| **Gharaie 2023, Sci Rep (GSE233512)** | Visium v1, 55 µm, sagittal | bilateral IRI 30 min, male B6, 7–8 wk | control, 24 h | **2** on GEO (4 of 8 areas) | none | full Space Ranger set |
| Xuanyuan 2025, Nat Commun (GSE269719, GSE269622, figshare 28761695) | Xenium 300 genes, cells; Visium CytAssist | bilateral IRI 18 min, male B6 | sham, 4 h, 12 h, 2 d, 14 d, 6 wk | **1 mouse** (left and right kidney) | Healthy/Injured S1, S2, S3, FR-PT | h5ad 2.0 GB; tarballs |
| DynaST-seq (GSE290361), unpublished | barcoded array 50 × 50, 4sU | renal artery clamp 50 min; sex not stated | ~0, 5, 13, 25 h | 2 (unit unknown) | unknown | h5ad |
| Dixon 2022, JASN (GSE182939) | Visium v1 | bilateral IRI 34 min, **female** B6 | sham, 4 h, 12 h, 2 d, 6 wk | **1** | PTs12 / PTs3 | Space Ranger outs |
| Melo Ferreira 2021, JCI Insight (GSE171406) | Visium v1 | bilateral IRI 22 min, 129/SvEv, male | sham, 6 h | **1** | PT(S1/S2), PT(S3-C), PT(S3-OS) | counts and spatial |
| Polonsky 2024, Nat Commun (Zenodo 12709329) | seqFISH 1,300 genes, cells | bilateral IRI, male B6 | control, **d28** | **3** | PTS1/2/3 | CSV |
| Rietjens 2026, JASN (GSE298899) | Stereo-seq | bilateral IRI | sham, **d14** | 4 | PT-S1/S2 | yes |
| Wu 2024, Nat Commun (GSE227044/5) | HybISS 200 genes | bilateral IRI 23 min, male | sham, 4 h, 12 h, 2 d, 6 wk | **1** | per sample | segmentation |
| GSE299516 (GeoMx) | GeoMx ROIs | **unilateral** IRI, male | d1, d7, d14 | 2 | cortex vs medulla | DCC |

**Excluded:**
- Janosevic 2021: a single CLP section.
- Marshall 2022: no AKI model.
- Wang 2026: rat, n = 1.
- Singh 2026: ex vivo cold storage.
- Wijaya 2025: rat cisplatin, restricted data.

## The two used in notebook 105

**Cheung 2022.**
- Fig. 1A legend: "There were 3 biological replicates per time point."
- "Male C57BL/6J mice, 12–16 weeks of age"; the clamps were removed "After 19 minutes".
- One 10 µm section per mouse, whole kidney (cortex to papilla).
- The spot-position file is not deposited. Every barcode maps onto the fixed Visium v1 grid (checked on all
  sections).
- Two day-1 sections are shallow (median 2,113 and 1,133 UMI per spot).

**Gharaie 2023.**
- "We profiled 2 biological replicates (i.e. 2 mice) per condition, and for each biological replicate we performed 2
  technical replicates."
- GEO holds 4 areas: control 1 and 2, IRI 1 and 2. That these are different mice is likely but not stated.
- Sagittal sections through the middle axis.

## Cautions

- **Injury lowers PT segment markers.** Wu 2024 (Results): combined clustering lost PTS1/2/3 "because injury causes
  downregulation of PT segment-defining markers". A marker axis can misplace injured S3. Use injury-stable positional
  genes, or an anatomical coordinate.
- **Spot mixing.** Visium spots are 55 µm mixtures. Outer-stripe PT spots carry thick ascending limb RNA (Melo
  Ferreira 2021: "the PT S3 segment cluster … contained gene expression from the neighboring TAL").
- **Lab independence.** Dixon, Wu and Xuanyuan share the Humphreys lab's protocols with Kirita 2020. Cheung (UAB)
  and Gharaie (Johns Hopkins) are independent.
- **Descriptive precedent.** Xuanyuan 2025 reports injured PT "initially localized within the corticomedullary
  junction and progressively expanded into the cortical region by 12 hours" (one mouse per time point).
