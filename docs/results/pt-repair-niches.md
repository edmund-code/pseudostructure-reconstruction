# What defines the neighbourhoods of proliferating and transport-losing tubules? (notebook 110)

**The question.**
- In the IR mice, the cell-cycle (CYC) and transport (TPT) residuals cluster among neighbouring tubules beyond
  position and depth (notebook 96, Moran's I 0.08–0.11 above the controls).
- Both programs are patchy across tubules, whereas the injury markers are near-universal (notebook 107).
- Is there a repair niche, meaning a neighbourhood that predicts which tubules proliferate or lose transport? If so,
  what is in it?

**What the literature predicts.** A scout read 33 of 38 papers in full; the summary is in
`docs/research/repair-niches-literature.md`.
- **Timing.** Tubular proliferation starts about 24 h after IRI and peaks about 48 h, dominated by the outer stripe.
  - Witzgall 1994: PCNA⁺ 80% outer stripe against 23% cortex.
  - Humphreys 2008 and 2011.
  - Kirita 2020: not cycling at 12 h; a cycling "repairing PT" at 2 d.
- **Who divides.** Dividing cells are the injured, dedifferentiated ones.
  - 72–75% are Kim-1⁺ (Humphreys 2011).
  - 89% of SOX9⁺ cells are HAVCR1⁺, and 40% of SOX9⁺ cells are Ki67⁺ (Kumar 2015).
  - The most severely injured cells arrest (Kirita 2020's severe-injury arrest signature; Witzgall's clusterin⁺ cells
    are never PCNA⁺).
- **Clones do not make clusters among tubules.**
  - Rinkevich 2014: clones "did not cross to adjacent segments within a nephron … nor to neighboring nephrons".
  - Berger 2014: "one would expect EGFP-positive cells to be arranged in clusters … which was not the case".
  - Chang-Panesso 2019: 87% of clones are single cells at day 2.
- **Proliferation clusters near injury within a nephron** (Bordoni 2023, intravital imaging):
  - more cycling cells within 25 µm of necrotic nuclei (2.27% against 1.62%);
  - a downstream wave in the same nephron, driven by casts, not by tissue proximity (laser-injured tubules 0–35 µm
    away were no more proliferative than those 85–135 µm away).
- **The interstitium at 12 h–2 d is fibroblasts and, in the outer stripe, neutrophils; macrophages come later.**
  - Fibroblasts peak at day 2 (Kirita 2020; Wu 2024; Xuanyuan 2025: fibroblasts 7.1% to 16.1% while immune cells stay
    at 3.3%; Polonsky 2024's fibroblast microenvironments score highest at day 2).
  - Neutrophils sit in the outer stripe: 52.2% at 6 h (Melo Ferreira 2021); 76% of neutrophil nuclei come from 12 h
    (Lake 2023 mouse projection).
  - Macrophage dependence of proliferation is a day 3–7 effect (Lee 2011, Lin 2010, Zhang 2012).
  - "Kim-1 positive tubules at 6 weeks … have increased F4/80 positive signal compared with sham and 12-hour time
    points" (Dixon 2022).
- **No study has quantified a neighbourhood effect on tubular proliferation at 12 h–2 d.** The closest is
  Xuanyuan 2025: the injured-PT neighbourhood has the highest spot-level G2M score at day 2.

Frozen on 2026-10-18, before any association between a tubule's programs and its neighbours or interstitium was
computed. Before freezing, the scout tabulated marker detection fractions per specimen and read notebook 104's edge
enrichment. Examples: Mki67 in 15–18% of IR units; Adgre1 in 1–2% of all units; S100a8 and S100a9 only in IR.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/110_pt_repair_niches.ipynb`; logic `110.repair_niches.1`; seed 110.
- Results `results/pt_repair_niches/`.
- New helper `peritubular_ring` in `pseudospace/contamination.py`: for each bin outside every polygon, the nearest
  polygon within a distance. Tested.
- **Units:** the AKI export in the common support (14,048), with notebook 96's scores, residuals (`e INJ`, `e TPT`,
  `e CYC`), depth percentile and centroids.
- **Unit diameter:** the export's `equivalent_diameter_um`.
- **Bins:** as in notebook 104 (the AKI run's stage-01 manifest and bin-to-unit tables).
- **Pixel size:** each specimen's `microns_per_pixel`.

## 2 · Measurements (fixed now)

- **Neighbours:** k nearest PT units of the same specimen by centroid, k ∈ {3, 6, 12, 24}. The neighbour injury of a
  unit is the mean `e INJ` of its k neighbours.
- **Touching pairs:** centroid distance < 1.25 × the pair's mean equivalent diameter.
  - Non-touching neighbours lie at 1.25–3 mean diameters.
  - A unit's non-touching neighbour injury is their mean `e INJ` (missing if none).
- **Peritubular pool of a unit:** its edge bins (depth ≤ 1), plus the outside bins whose nearest polygon is that
  unit within 3 bins (6 µm).
- **Lineage scores:** log1p of 10⁴ × (counts on the lineage's markers) / (total counts) in the peritubular pool,
  standardised on the control units. Markers present in the export:
  - fibroblast/stromal: Col1a1, Col1a2, Col3a1, Dcn, Pdgfrb, Pdgfra, Postn, Timp1, Tnc, Acta2, Crlf1, Npr3, Igfbp5;
  - myeloid: C1qa, C1qb, C1qc, Cd68, Lyz2, Csf1r, Itgam, Adgre1, Cd74, Mrc1, Fcgr1, Ms4a7;
  - neutrophil: S100a8, S100a9, Csf3r, Cxcr2, Mmp8, Retnlg. Lcn2 is excluded because it is tubular.
  - endothelium: Pecam1, Emcn, Plvap, Kdr, Flt1.
- **Hypoxia** (the unit itself, as notebook 96 scores programs): Egln3, Bnip3, Ndrg1, Adm, Slc2a1, Ldha, Pdk1, Ankrd37,
  P4ha1, Vegfa. Residual against the control expectation at position, as `e`.
- **CYC halves:** G2M (Hallmark G2-M) and E2F (Hallmark E2F), each scored and residualised as CYC.

## 3 · Tests and readings

Models are fitted per mouse by least squares on standardised predictors. Every test is calibrated by permuting the
predictor of interest among units within position-decile × depth-tercile strata of that mouse (2,500 permutations).

**D1 · Does neighbours' injury predict a tubule's proliferation, beyond its own injury?**
- Model: e CYC ~ e INJ + e INJ² + depth pct + neighbour e INJ (k).
- **Holds** when, at k = 3 and k = 6, the neighbour coefficient is > 0 with one-sided p ≤ 0.05 in both IR mice, and
  exceeds the coefficient in both control mice.
- **TPT:** the same model on e TPT, with the coefficient expected < 0.
- k = 12 and 24 are reported. An effect that appears only at k ≥ 12 is treated as suspect.
- A negative CYC coefficient is reported as support for the arrest branch, not tuned away.

**D2 · Which interstitial lineage marks the niche?**
- Model: e CYC ~ e INJ + e INJ² + depth pct + fibroblast + myeloid + neutrophil + endothelium + e hypoxia, per IR
  mouse.
- **The fibroblast niche holds** when the fibroblast coefficient is > 0 with p ≤ 0.05 in both IR mice, and is the
  largest positive lineage coefficient in both.
- **Pre-registered expectations:**
  - myeloid near null;
  - neutrophil, if anything, in late S3 and in IR2A2.
- The same model on e TPT is reported (hypoxia expected negative).

**D3 · Is the clustering a niche, or contact within a nephron?**
- D1 for CYC is refitted with non-touching neighbour injury only.
- **"Beyond touching" holds** when the coefficient stays > 0 with p ≤ 0.05 in both IR mice.
- Descriptive:
  - Moran's I of e CYC and e TPT at k = 6, on all pairs and with touching pairs removed;
  - the share of each unit's 6 neighbours that touch.

**Controls (named now).**
- **C1 · Compartment.**
  - CYC is rescored on core aggregates (bins ≥ 6 µm inside the edge, notebook 104).
  - Reported: ρ(core CYC, full CYC) and D1 on the core score.
  - Also reported: notebook 104's edge enrichment of the CYC genes.
- **C2 · Entry or division.** The E2F and G2M residuals along x (descriptive), and D1 on each half.
- **C3 · Timing.** The share of IR units with Mki67 or Top2a detected, against controls. This places our mice in the
  12 h–2 d window (descriptive).

## 4 · Outputs and figures

- **Tables:**
  - `units.csv` (unit, specimen, x, depth, scores, residuals, neighbour terms, lineage scores, touching share);
  - `d1.csv`, `d2.csv`, `d3.csv`, `controls.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 110a: D1 coefficients by k, IR against control mice, with permutation bands;
  - 110b: D2 lineage coefficients per IR mouse;
  - 110c: Moran's I with and without touching pairs;
  - 110d: a section map of e CYC and the fibroblast score, one IR mouse.

## 5 · Caveats

- One section per mouse: neighbours are 2-D. A convoluted tubule cut several times in a section gives several
  similar units, which is why D3 is there.
- The peritubular pool mixes interstitium with the neighbours' edges. Lineage scores are proximity measures, not cell
  counts.
- Two IR mice. The late-S3 interstitial rise is seen in one (IR2A2).
- T and B cells are not measurable (Cd3e and Ms4a1 are about 0), and are not tested.
