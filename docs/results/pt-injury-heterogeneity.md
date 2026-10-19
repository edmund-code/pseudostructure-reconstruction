# Does injury make tubules unequal? Position-resolved differential variability (notebook 107)

**The question.** Every AKI analysis so far compares means: how much a gene or program changes at a position along
the PT. Injury could also change how *unequal* tubules are at the same position. Two pictures differ here:
- **Patchy injury** (some tubules hit, others spared) raises tubule-to-tubule variance, most where about half the
  tubules respond.
- **Graded injury** (every tubule shifted) leaves the residual variance unchanged.

Notebook 96 found one Gaussian per position for the injury score (graded). It also found that proliferation and
transport loss cluster in neighbourhoods of tubules, which hints at patchiness for those programs. A variance test
asks the question gene by gene and program by program.

**What the literature offers.** The methods scout read these in full, section C of
`docs/research/methods-scout-contamination-networks-variability-evolution.md`:
- **Variance follows the mean, so differential variability must be residual.**
  - BASiCS (Eling et al. 2018) defines residual over-dispersion as log dispersion minus a smooth trend on the mean.
  - memento (Kim et al. 2024) uses residual variance after a mean–variance fit.
- **Martinez-Jimenez et al. 2017:** in old mice, more variable activation came from a smaller responding fraction,
  not noisier responders. So detection fraction should be reported beside dispersion.
- **DiffVar (Phipson & Oshlack 2014)** treats samples as replicates. No gene-level test is valid at 2 against 2; the
  published methods need about 8–50 samples per group.
- **Inference is therefore at the specimen level**, against balanced relabeling. A cell-level Levene test on 70 cells
  from 11 mice gave 169 DV genes where pseudobulk found none (Jiang 2026).
- **No kidney injury study has tested tubule-to-tubule variance.**

Frozen on 2026-10-18, before any dispersion was estimated.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/107_pt_injury_heterogeneity.ipynb`; logic `107.heterogeneity.1`; seed 107.
- Results `results/pt_injury_heterogeneity/`.
- Module `pseudospace/heterogeneity.py`, tests `tests/test_heterogeneity.py`.
- **Units:** the AKI export, in the common support (14,048 units; `support_positions`), with raw counts and
  `n_spots`.
- **Genes:** notebook 71's AKI universe (`universe_aki_vs_control.csv`).
- **Programs:** notebook 96's INJ (18 markers), TPT (93 transporters) and CYC (297 cell-cycle genes).
- **Position bins:** 8 equal-width bins of x.
- **Relabelings:** the two balanced specimen relabelings, {IR2A2, Ctrl1A4} against {IR2A4, Ctrl1A2} and {IR2A2,
  Ctrl1A2} against {IR2A4, Ctrl1A4}.

## 2 · Estimators (fixed now)

For each specimen s, bin b and gene g, over the units u in that cell:
- **Expected count:** μ_u = N_u · Σy / ΣN, with N_u the unit's total count.
- **Moment dispersion:** α̂ = Σ[(y − μ)² − μ] / Σμ², floored at 10⁻³. Computed when the cell has ≥ 30 units and a
  mean count ≥ 0.5 per unit; otherwise missing.
- **Residual over-dispersion:** ε_gsb = log α̂ − f_sb(log mean count). Here f_sb is a lowess trend (frac 0.3) over
  the genes of that specimen and bin.
  - Subtracting the trend removes shifts shared by all genes: section quality and unit-size distribution.
  - The shift of f_sb itself, IR against control, is reported as global variability.
- **Detection fraction:** the share of units in the cell with y > 0.
- **Program value:** the mean over a program's genes with finite ε.

**Sensitivity (unit size):** the same, using only units whose `n_spots` lie within the pooled interquartile range.

## 3 · Readings

**R1 · Does injury make a program more heterogeneous?** For INJ, CYC and TPT, in each bin:
- Δε = mean ε over IR minus mean ε over controls.
- The bin **holds** when all four IR–control specimen pairs are positive, and Δε exceeds the absolute group
  difference under both balanced relabelings.
- **A program becomes heterogeneous** when it holds in ≥ 3 of the 8 bins.

**R2 · Patchy or graded?** For each program and each IR specimen, the Spearman ρ across bins between ε change and
p(1 − p).
- p is the share of the IR specimen's units whose program score exceeds the controls' 95th percentile at that
  position (scores as in notebook 96).
- Patchiness predicts ρ > 0, with excess variance where about half the tubules respond.
- Reported, not judged (8 bins).

**R3 · Genome-wide (exploratory).**
- For each gene, Δε averaged over bins with finite values.
- Its standing against the two relabelings: |Δε| above both.
- The 50 genes with the largest positive Δε are listed with their programs, and with their detection-fraction change.
- The count of genes beyond both relabelings in each direction is compared between the condition and the
  relabelings.

**Recovery expectation (named now).** Notebook 96 found the INJ score graded (one Gaussian per position). So INJ is
expected **not** to become heterogeneous, while CYC, which clusters in neighbourhoods, may. A heterogeneous INJ
program would contradict notebook 96 and is reported as such.

## 4 · Outputs and figures

- **Tables:**
  - `cells.csv.gz` (specimen, bin, gene, units, mean, α̂, ε, detection);
  - `trend_shift.csv`, `programs.csv`, `r1.csv`, `r2.csv`, `genes.csv.gz`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 107a: ε against log mean for one control and one IR specimen, with the trends;
  - 107b: program Δε along x, the condition against the relabelings;
  - 107c: Δε against detection change for program genes.

## 5 · Caveats

- Moment dispersion is noisy for low counts; the ≥ 0.5 mean-count rule limits this.
- Coordinate error adds variance to steep genes: Var(e) ≈ g′(x)² σ_x². The trend does not remove this when steep
  genes cluster in a program. TPT is steep along x, so read its ε with care.
- The deep S3 units the package dropped (notebook 96) are absent, which biases IR variance down at the S3 end.
- Two IR mice; balanced relabeling gives only two alternative groupings.
