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
