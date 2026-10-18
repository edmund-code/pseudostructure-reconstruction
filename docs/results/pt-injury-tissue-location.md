# Nephron position or tissue location: what sets injury between tubules (notebook 96)

**The question.** Our AKI analyses describe injury as a function of one coordinate: where a tubule sits along the PT.
In ischemia–reperfusion, injury "is most apparent in the S3 segment of the proximal tubule in most animal models of
ischemia". Blood flow to the outer medulla "is reduced disproportionately to the reduction in total kidney perfusion",
and "the anatomy of the capillaries in the outer medulla makes them very vulnerable to occlusion". These quotes are
from Bonventre & Yang 2011, JCI 121:4210, PMID 22045571, read in full at jci.org/articles/view/45161.

So S3 vulnerability may be partly a property of the cell (nephron position) and partly of the tissue around it: depth,
and the neighbouring vessels and tubules. Our tubule units have a nephron position *and* a place in the section.
Beyond nephron position, does a tubule's tissue location explain its injury?

The answer matters for the method. If tissue location adds nothing, the continuous nephron coordinate is a sufficient
description of injury. If it adds something, the pseudospace should be read beside the tissue map.

Protocol frozen on 2026-10-18, before any injury score, depth or neighbour statistic was computed on these units.
Notebook 76's code for depth (distance to glomeruli) and its placement of the dropped units were read, not rerun.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/96_pt_injury_tissue_location.ipynb`; logic `96.tissue_location.1`; seed 96.
- Results `results/pt_injury_tissue_location/`.
- Module `pseudospace/tissue_location.py`, tests `tests/test_tissue_location.py`.
- **Units:** the AKI export `results/package_pseudospace/aki/pt_units.h5ad` (pass-1 fingerprint `bb4a80c1b484`),
  through `package_inputs.aki_inputs`, as notebook 93 builds it:
  - the common-support units and their positions (min–max scaled, x ∈ [0, 1]);
  - specimen and package segment labels;
  - counts.
  Centroids (`x_centroid`, `y_centroid`) come from the export's obs.
- **Glomeruli:** stage-01 units with `seg_class == 'glomerulus'`
  (`results/package_pseudospace/runs/aki/01_units/{specimen}_units.h5ad`).
- **Sensitivity units:** notebook 76's `placement_units.csv`, the c15 and c19 units with `on_pt_stretch`. These are
  the injured S3 units the package dropped (`pt-package-issues.md`, item 7). Their placed position is
  `median_value`, scaled with the same support. Their counts come from
  `runs/aki/03_integration/pass1_labelled.h5ad` (layer `counts`).
- **Programs (fixed by notebook 76, never chosen here):**
  - INJ: 76's 18 injury markers (`positive_control_markers.csv`), the primary score;
  - TPT: 76's transport program T_PT (93 genes, `programs_transporters.csv`), secondary;
  - CYC: 76's cell-cycle program C (Hallmark G2-M ∪ E2F in U76, 297 genes), secondary.

## 2 · Definitions

- **Unit expression:** log1p(count / total counts × 10⁴). Total counts are the unit's counts over all export genes.
- **Program score:** the mean over the program's genes of (expression − control mean) / max(control SD, 0.25), with
  the mean and SD over the pooled control units.
- **Injury effect** e = score − m(x), where m is the control expectation at the unit's position: the mean score of
  control units in the 0.05-wide position bin, pooled over both controls.
  - **Cross-fitting for controls:** a control unit's m comes from the other control mouse only, so its own noise is
    not in its expectation.
- **Depth** (notebook 76's rule): the mean distance to the 3 nearest glomerulus polygons of the same section,
  expressed as a percentile among that section's units in the analysis. 0 is next to glomeruli, 1 the deepest.
- **Position bins for nulls:** 20 equal-count bins of x within each section.
- **Neighbours:** the k = 6 nearest units of the same section by centroid, with binary weights, row-standardised.

## 3 · Analyses

**A1 · Nephron position (recovery check, H1).**
- In each IR mouse, Spearman ρ of e with x.
- Holds if ρ > 0 in both IR mice for INJ, matching notebooks 76 and 91 (injury S3-centred).

**A2 · Depth beyond position (H2).**
- In each mouse, ordinary least squares: e ~ natural cubic spline of x (5 df) + depth.
- **Statistics:**
  - β_depth, the change in e from the shallowest to the deepest unit;
  - its partial R².
- **Null:** depth permuted within position bins within the section, 1,000 times. p = (1 + #{|β_null| ≥ |β|}) / 1,001.
- **Controls** use cross-fitted e. They give the anatomical baseline: in healthy tissue, deeper PT at the same position
  may differ for reasons other than injury.
- **H2 holds** if, for INJ:
  - β_depth > 0 with p ≤ 0.05 in both IR mice; and
  - each IR β exceeds both control β by ≥ 0.25 score units.
- **Within S3:** the same model on S3 units only, reported beside the primary.

**A3 · Neighbourhood beyond position and depth (H3).**
- Residual r = e − the A2 fit.
- Moran's I of r on the neighbour graph, per mouse.
- **Null:** r permuted within position bins within the section, 1,000 times.
- **H3 holds** if, for INJ:
  - I exceeds the null's 97.5th percentile in both IR mice; and
  - each IR I exceeds both control I by ≥ 0.05.
- **Correlogram:** I at k = 3, 6, 12, 24 neighbours, reported for scale.

**D4 · Graded or all-or-none (descriptive).**
- Within S3, and within each position tertile, the BIC of a one-component against a two-component Gaussian mixture of
  e, per mouse.
- A two-component fit counts as preferred when its BIC is lower by more than 10. The mixing fraction and component
  means are reported.

**Secondary programs.** A1–A3 for TPT (expected sign negative) and CYC. Reported; H2 and H3 are judged on INJ only.

**Sensitivity.**
- A1–A3 with the c15 and c19 units on the PT stretch added (IR and control alike).
- The verdict is "robust" if the H2 and H3 readings do not change, and "changes" otherwise. The primary stays the
  export's units.

## 4 · Reading rules

- Depth and x are correlated in controls: Spearman 0.43 and 0.62 in control S3 (notebook 76). A2 asks about depth
  *at fixed* x, so its power rests on depth varying at a given x. The depth SD within position bins is reported.
- The package dropped many injured deep-S3 units: deep S3 IR 251 and 32, against control 783 and 700. The primary
  units therefore under-sample the most injured, deepest tubules. The sensitivity arm exists for this reason.
- x/y enter only this downstream analysis, never the reconstruction (the method stays tissue-agnostic).
- Two IR mice: every holding statement needs both.

## 5 · Outputs and figures

- **Tables:**
  - `units.csv` (unit, specimen, x, segment, depth, scores, e);
  - `a1_position.csv`, `a2_depth.csv`, `a3_neighbours.csv`, `a3_correlogram.csv`, `d4_mixtures.csv`;
  - `depth_spread.csv`;
  - `sensitivity.csv`, `verdicts.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 96a: e along x per mouse, coloured by depth;
  - 96b: β_depth and Moran's I per mouse, against their null bands;
  - 96c: the sections' PT units in x/y, coloured by the INJ residual r, one panel per mouse;
  - 96d: correlograms.

## Notebook 96 results

`analysis/notebooks/96_pt_injury_tissue_location.ipynb` (logic `96.tissue_location.1`) writes to
`results/pt_injury_tissue_location/`. It followed the protocol frozen in 2f8167f. Guards G1–G2 passed: 14,048
common-support units, and programs INJ 18, TPT 93 and CYC 297. Runtime is 24 s. Module `pseudospace/tissue_location.py`
with `tests/test_tissue_location.py` (3 tests).

**How the protocol was read** (`how_the_protocol_was_read.csv`):
- the unit library covers all 19,059 genes in both arms;
- the control expectations are built from export control units only;
- the placed units count as S3 beyond the controls' S2/S3 cut, because they have no package label;
- the correlogram nulls use A3's position bins.

**Units.**
- Ctrl1A2 2,680, Ctrl1A4 2,686, IR2A2 3,466, IR2A4 5,216.
- The sensitivity arm adds 31, 26, 274 and 258 placed c15 and c19 units on the PT stretch.

**Depth varies at fixed position.**
- The median within-bin SD of the depth percentile is 0.20 in the controls and 0.27–0.29 in IR.
- Depth tracks position far more in the controls (Spearman 0.73, 0.69) than in IR (0.35, 0.19). The IR sections lost
  most deep S3 units (notebook 76), so their remaining S3 units sit at all depths.

**Verdicts (INJ):**

| Hypothesis | Ctrl1A2 | Ctrl1A4 | IR2A2 | IR2A4 | Holds | With placed units |
|---|---|---|---|---|---|---|
| H1 nephron position: ρ(e, x) | −0.059 | 0.033 | **0.715** | **0.710** | **yes** | yes (robust) |
| H2 depth at fixed position: β (p) | −0.056 (0.003) | −0.022 (0.19) | −0.089 (0.02) | −0.027 (0.34) | no | no (robust) |
| H3 neighbours: Moran's I (null 97.5%) | 0.106 (0.021) | 0.064 (0.021) | 0.141 (0.020) | 0.112 (0.017) | no | no (robust) |

- **H2 fails clearly.** At fixed position, deeper tubules are not more injured. β is negative or null in IR, on a
  score whose IR mean effect is 2.2–2.7.
  - Within S3 alone: β 0.114 (p 0.32) and 0.037 (p 0.61).
  - TPT shows a small positive depth slope, meaning less transport loss deeper (0.044 and 0.048, p 0.001). Ctrl1A2
    shows the same (0.061, p 0.002), so this is anatomy, not injury.
- **H3 fails the frozen margin.**
  - Injury residuals cluster among neighbouring tubules in every section: I 0.06–0.14, against a null of 0.02.
  - IR exceeds the controls by 0.006 to 0.077, short of the 0.05 margin over both controls (IR2A2 − Ctrl1A2 = 0.035).
  - Neighbour similarity in healthy tissue (shared ambient RNA, spill-over, section effects) is as large as injury
    adds.
  - The correlogram falls with k in every section (k 3 → 24: IR2A2 0.166 → 0.096, Ctrl1A2 0.126 → 0.072). The
    clustering is short-range, not territorial.
- **Secondary programs cluster more under injury** (not judged):

  | Moran's I at k 6 | Ctrl1A2 | Ctrl1A4 | IR2A2 | IR2A4 |
  |---|---|---|---|---|
  | CYC | 0.052 | 0.066 | 0.140 | 0.144 |
  | TPT | 0.043 | 0.060 | 0.151 | 0.107 |

  Proliferation and transport loss come in neighbourhoods of tubules beyond position and depth: 0.08–0.11 above the
  controls.
- **D4, graded.** Within S3 and within every position tertile, one Gaussian describes the IR injury effect as well as
  two: no two-component BIC gain greater than 10 in either IR mouse. At the resolution of whole tubule cross-sections,
  injury at a given position is a graded shift, not a mix of spared and injured tubules. The controls prefer two
  components in S3, from a small high-INJ tail: 4% of Ctrl1A2's S3 units.

**What this answers.**
1. **Nephron position, not tissue depth, organises injury between tubules.**
   - At the same position along the PT, depth adds nothing.
   - The S3 vulnerability seen here is carried by position along the nephron. The outer-medullary blood-flow account
     would predict a depth gradient at fixed position, and we do not see one.
   - The deep, most injured S3 units the package dropped are under-sampled. The placed-unit arm, which adds 532 of
     them, gives the same verdicts.
   - The continuous nephron coordinate is a sufficient axis for the injury program.
2. **Injury is graded at the tubule level.** Within position there is one shifted population, not spared and injured
   tubules.
3. **Proliferation and transport loss cluster in small neighbourhoods** beyond position and depth. The injury markers
   themselves cluster only as much as healthy neighbours do. This is a lead for a spatial analysis of repair niches,
   reported, not judged.
4. **Limits:**
   - two IR mice;
   - depth from glomerulus distance is a proxy for the cortico-medullary axis;
   - the IR sections hold few deep S3 units;
   - whole cross-sections average many cells, so cell-level patchiness inside a tubule cannot be seen.
