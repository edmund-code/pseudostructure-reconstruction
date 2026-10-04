# Methods: pathway analysis along the PT coordinate (draft for `docs/paper/draft.md`)

Produced by notebooks 12 (gene model), 31 (strategy comparison), 33 (registration and region),
34 (external segment data) and 37 (final analysis). Notebook 37 recomputes every number below from
one PT coordinate, notebook 13's scFates coordinate by default, in a single run.

## Gene universe and gene sets

PT structures, their reviewed segment labels and the shared coordinate come from notebook 13. We
keep structures inside the 1st–99th percentile of every specimen's coordinate range. Raw counts are
rebuilt for the accepted one-to-one ortholog panel and normalized as log1p(count / library × 10⁴).
The library size is the total over orthologs measured in all four specimens.

A gene is eligible when it is measured in both species and detected in at least 2% of structures,
averaged with equal weight over specimens. Gene sets come from Reactome 2022, MSigDB Hallmark 2020
and KEGG 2019 Mouse, mapped through the ortholog table. We test sets with 10–300 eligible members.

## Gene-level model

Each gene is fitted by weighted least squares with a fixed cubic B-spline basis on the coordinate
(6 df, knots at quantiles). Specimens enter as within-species sum-to-zero intercept contrasts, and
weights give every specimen equal total weight. Three nested models are fitted:
- a shared curve;
- the shared curve plus a constant species offset;
- the shared curve plus a species-specific curve.

The partial-F statistics T_level (offset), T_spatial (species × position beyond the offset) and
T_total are used **only to rank genes**. Their residual degrees of freedom count structures, not
specimens, so they are not tests. The species difference Δ(s) = human − mouse is evaluated on a
61-point grid with HC3 standard errors.

## Pathway test (joint test)

**Matched null.** For each gene set we compute the rank-AUC of its members' T_spatial against all
other eligible genes. We compare it with 9,999 random sets of the same size, drawn within strata
defined by tertiles of mean expression, detection and positional coverage. The matched z is the
observed AUC minus the null mean, divided by the null SD.

**Correlation.** The matched null treats genes as independent. Following CAMERA (Wu & Smyth 2012),
we estimate each set's mean pairwise residual correlation ρ from the full model's residuals:
within each specimen, then Fisher-averaged with equal weights and floored at zero. This gives a
variance inflation VIF(m, ρ) for the rank sum.

**Joint test and list.** The joint z is the matched z divided by √VIF, with a one-sided normal p.
BH runs over all 1,513 sets. A set is a candidate when its effect is positive and joint q ≤ 0.05.
The same inflation is applied to every contrast and to every sensitivity run on the same genes.
Each gene-subset check uses the inflation estimated for its own family of reduced sets.

## Specificity control: balanced specimen relabeling and target–decoy calibration

The four specimens (two mice, two cortex sections from one human donor) have exactly three
balanced 2 + 2 splits. These are the species split and two relabelings, each pairing one mouse
with one human section.

**Rerunning every strategy.** Each strategy is rerun unchanged on both relabelings. In a relabeled
model the species curve enters as a nuisance term, so the tested term is a difference between
specimens of the same species (Hou et al. 2023 used the same random-partition control). The two
within-species comparisons are added as further controls.

**Why the relabelings are a valid null.**
- Write a specimen's deviations as an offset a_ij and a positional deviation b_ij(s), independent
  across specimens.
- The species split and both relabelings are the three orthogonal ±½ contrasts of the four
  specimens. Under that independence their noise has equal variance, so each relabeling is one draw
  of the species contrast's section-level null.
- The design is a split plot (Altman & Krzywinski 2015), so the species offset is judged against
  specimens within species (2 df). The species × position term is judged against specimen ×
  position deviations (2(k − 1) df for k positions).
- A term shared by both human sections, such as a donor effect, cancels in the relabelings. It is
  therefore not covered by this control.

**Target–decoy calibration.** For each strategy, the false-discovery proportion at threshold t is
estimated as FDP(t) = mean relabeled count ≥ t ÷ species count ≥ t (Elias & Gygi 2007). It is made
monotone as a q-value. We report the number of species calls at FDP ≤ 10%, with intervals from
resampling pathways (300 draws), and the minimum attainable FDP. Per-pathway inference at specimen
level is not possible with two decoys: a 2-df test produced no calls.

**Strategies compared:**
- whole-PT and per-segment DESeq2 pseudobulk with signed preranked GSEA (gseapy, multilevel);
- pathway-score models tested at structure level;
- gene-GAM offset, total, step-shaped and smooth positional statistics, each with the matched and
  the joint test;
- a size-only null;
- preranked GSEA on T_spatial.

The negative-control discovery ratio (mean relabeled ÷ species calls at BH 0.05) equals the decoy
FDP at a strategy's own threshold.

## Building the reporting set

1. **Candidates.** Joint q ≤ 0.05, positive effect.
2. **Not specimen-sensitive.** Not called by the joint test in either relabeling or either
   within-species comparison.
3. **Stable.** Retained in at least 5 of 7 planned refits: each specimen left out, spline basis 4
   and 8, and a 5% detection rule. Each refit recomputes the matching strata.
4. **Robust.** Passes both of the following.
   - **Registration.** Each human section is warped piecewise-linearly so that its reviewed
     S1→S2 and S2→S3 transitions (logistic fits per specimen) fall on the mouse transitions. All
     genes are refitted. The check passes when the effect ratio is ≥ 0.7 and joint q ≤ 0.10.
   - **Region.** Mouse S3 deeper than the 95th percentile of the same specimen's S1 is removed;
     depth is the mean distance to the three nearest reviewed glomeruli, with centroids verified
     against the segmentation. As a power-matched control, the same number of mouse S3 structures
     is removed at random 20 times. A pathway is region-sensitive if its effect turns non-positive
     or its effect ratio falls below the 5th percentile of the power-matched ratios.
5. **Core.** Robust, and also passes two gene-subset checks: restricting to genes detected in at
   least 5% of structures in each species, and removing genes sex-biased in any mouse PT segment
   (Xiong et al. 2023 public tables; |log2FC| ≥ 1, padj ≤ 0.05). Each check passes when the effect
   ratio is ≥ 0.7 and joint q ≤ 0.10.

**Further annotations, not filters:**
- removing ambient, neighbouring-segment and shared-exon-mapping genes flagged by the literature
  check;
- restricting to genes with equal probe counts in the two Visium HD panels (workstream 3's table);
- 5-bin covariate matching;
- notebook 12's earlier run on a different trajectory method (DPT);
- the specimen-level split-plot test below.

## Programs and positional descriptors

After the list is frozen, robust pathways are grouped into programs.
- **Contributing members** of a pathway are its members in the top 10% of T_spatial; a pathway with
  fewer than three uses all members.
- **Grouping.** Average linkage on 1 − overlap coefficient of contributing members, cut at 0.5.
- **Naming.** Each program is named after its lowest-q pathway and listed with its shared drivers,
  the contributing genes present in at least half its pathways.

Descriptors are display rules, not tests:
- **Local divergence.** D(s) = AUC of the members' |Z(s)| against other genes, minus 0.5; we report
  its peak position and the median member Z at the peak.
- **Shape class** of the median member Δ(s): *reversing*, *localized* or *graded*.
- **Averaging loss.** The median over contributing members of 1 − |mean_s Δ| ÷ mean_s |Δ|, which is
  the share of the positional difference that a whole-PT mean cancels.

## Specimen-level split-plot test

**Pseudobulk.** Raw counts are summed per specimen × coordinate bin: 8 equal-width bins on the
common support, keeping cells with at least 20 structures.

**Gene test.** log2-CPM values get voom-style precision weights (Law et al. 2014). Each gene is
fitted with: specimen + cubic spline(position, 3 df) + species × spline. The 3-df interaction is
tested with a moderated F (empirical-Bayes variance moderation, Smyth 2004). Its residual df come
from specimen × position deviations. limma is not available in the analysis environment, so a
Python implementation (`pseudospace/split_plot.py`) is used and checked against simulations.
Relabelings add species × spline as a nuisance.

**Calibration** is judged under relabeling by the median-F inflation λ and the fraction of genes at
p ≤ 0.01.

**Pathways.** Matched rank-AUC on the moderated F, with the inflation estimated from pseudobulk
residual correlations.

## External replication

**Data.** Public snRNA, summed per donor and segment group:
- **Census human renal cortex** (6 donors, convoluted PT versus S3); pre-specified as the primary
  human reference.
- **Lake et al. 2023 / KPMP healthy cortex** (CELLxGENE v1.5): 7 donors with ≥ 50 nuclei in each
  of PT-S1, PT-S2 and PT-S3; S1 and S2 are pooled as early. Genes are matched by Ensembl id, with
  the symbol as fallback.
- **Census mouse** (12 male and 12 female donors, S1 + S2 versus S3).
Per donor, Δ = log2CPM(S3) − log2CPM(early). Our own contrast, from reviewed segment labels, is
(human S3 − early) − (mouse S3 − early) averaged over specimens.

**Informative genes.** Detected in at least 10% of our structures in each species, with mean
log2CPM ≥ 1 externally in both species.

**Global component.** A global component is measured explicitly: the median contrast across
informative genes, and the slope of the interaction on the mouse gradient (the signature of
attenuated human positions). The primary analysis removes the median per donor and per specimen.

**Test.**
- Per gene, a Welch t compares human and male-mouse donors.
- Signing it by our direction gives an agreement score per gene.
- Each pathway gets a matched rank-AUC of that score against expression-matched random genes. The
  inflation comes from donor-level residual correlations, so donors are the replicates.
- BH runs within the robust list.

**Comparators.**
- uncalled pathways, and pathways called only under relabeling;
- female mice as reference;
- additionally removing the slope component;
- for pathways called only by conventional screens, the analogous average-level test (whole PT,
  signed by our DESeq2 direction).

## Software and reproducibility

Python 3.11, NumPy, SciPy, statsmodels, pandas, gseapy, pydeseq2. Stages are cached by parameters,
inputs and code (`pseudospace.stage_cache`). The notebook logic version is `37.pathway_final.1` and
the seed is 12.

## Caveats (printed beside every result)

- Two male control mice and two cortex sections from one male human donor; specimens are the
  replicates.
- Relabeling controls section-level noise only.
- Human sections are cortex, so mouse outer-stripe S3 has no human counterpart.
- The Visium HD probe panels differ between species.
- All results are descriptive.
