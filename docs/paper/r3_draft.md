# R3 draft: conserved and species-specific PT zonation programs (notebook 40)

Source: `analysis/notebooks/40_pt_zonation_classes.ipynb` (logic version `40.zonation_classes.1`).
Results are in `results/pt_zonation_classes/`; the main figure is `figures/fig_R3_zonation_classes`
and the table is `tables/table_R3_zonation_classes.csv`.

## Results (draft, about 400 words)

**Beyond a conserved transport core, PT zonation programs are largely species-specific.**

Reconstructed positions showed that human and mouse PT diverge mainly in which genes vary along
the tubule (notebook 39). We therefore called segment zonation in each species separately, using
two contrasts. S2 − S1 is cortical in both species. S3 − early uses only cortical-like mouse S3,
so the region matches the human cortex.

Each call was scaled to that species' own noise with a variance-moderated test, and confirmed in
an independent atlas of the same species (12 mouse snRNA donors; 7 healthy human cortex donors
from KPMP). We called a gene "species-only" only when the other species was confidently flat, not
merely non-significant, so the noisier human data cannot inflate the mouse-only class.

Of 7,407 probe-balanced genes measurable in both species, 155 were zonated in the same direction
in both (conserved) and 52 in opposite directions (reversal). Another 247 were zonated only in
mouse and 66 only in human; 1,297 were flat in at least one species and zonated in neither. Most
of the rest (5,590) could not be classified with two specimens per species.

Mouse-only genes outnumbered human-only genes under every alternative rule (Fig. Rb). The gap
narrowed, from 247 against 66 to 243 against 141, when the human effect threshold was scaled to
the lower human amplitude. Most of the excess came from S1 → S2 (219 against 44 genes). For the
region-matched S3 contrast the two classes were similar (48 against 35). Confirming mouse calls in
female rather than male donors kept most calls: female snRNA kept 214 of 247 mouse-only, 151 of 155
conserved and 44 of 52 reversal calls; female microdissection kept 149, 127 and 30.

The **conserved core** is apical solute transport: neutral and cationic amino-acid transport
(Slc3a1, Slc7a8, Slc6a18, Slc7a7), organic anions (Slc22a7, Slc22a12), phosphate (Slc34a3), Slc5a2
and Slc5a10, and the intrarenal renin–angiotensin genes (Agt, Ace, Enpep). These enrichments
survived finer expression matching and removal of low-expression genes.

**Mouse-only zonation** concentrates in metabolism: sterol synthesis (Cyp51, Lss, Ebp), peroxisomal
very-long-chain β-oxidation (Acox1, Acox2, Hsd17b4, Slc27a2), tyrosine catabolism and pyrimidine
catabolism. It also covers several genes highlighted earlier: creatine synthesis (Gatm),
glutathione synthesis (Gclc, Gclm) and the organic-anion uptake genes Slc22a6 and Slc13a3.

Human-only and reversal classes are small. Their pathway enrichments were mostly not robust to
expression, so the human-specific side is best described gene by gene:
- **reversals:** Dcxr, Ugt3a1, Vnn1;
- **zonated only in human:** Sel1l3, Pkhd1;
- **expressed and zonated only in human:** Rbp4, Aox1. These genes are absent from mouse PT.

Of the 36 robust position-dependent pathways (R2), 19 have mostly mouse-only zonated members, 15
mostly conserved members, 1 mostly reversal members and 1 none. Two of them (Steroid hormone
biosynthesis, Nitrogen metabolism) have fewer than 10 classified members. The pathway screen therefore recovers both the shared transport core
and the mouse-specific metabolic programs.

**Caveats** (keep beside the result):
- Two mice and two sections from one human donor.
- Calls require reference confirmation, and the human reference has fewer donors.
- Most genes remain indeterminate.
- Probe-imbalanced genes (1,622, of which 86 fall in a zonated class) and one-species-expressed genes are reported separately.

## Methods paragraph (draft)

**Zonation classes.**
- **Contrasts.** For each specimen, PT segment contrasts (S2 − S1, and S3 − early with mouse S3
  restricted to cortical-like structures) were computed as pseudobulk log2 CPM ratios of summed
  counts. Genes needed ≥ 20 counts in every specimen × segment pseudobulk of both species.
- **Species calls.** Within each species, the specimen-mean contrast was tested with a moderated
  one-sample t. Variances were shrunk within ten expression bins (limma-style squeezeVar), and BH
  was applied across genes.
- **Zonated, flat, indeterminate.** A gene was zonated when q ≤ 0.05 and |contrast| ≥ 0.5 log2,
  and an independent same-species reference agreed in direction (one-sided moderated donor-level
  p ≤ 0.05). The references were CELLxGENE Census mouse snRNA (12 male donors) and Lake et al.
  2023 healthy human cortex (7 donors). A gene was flat when its 95% interval lay within ±0.3
  log2; otherwise it was indeterminate.
- **Classes.** Conserved: both species zonated in the same direction. Reversal: opposite
  directions. Mouse-only or human-only: zonated in one species and flat in the other.
- **Alternatives.** No confirmation; no effect floor; human floor scaled by the noise-corrected
  amplitude ratio; flat margin 0.5; microdissection confirmation.
- **Pathway enrichment.** For each class, the class indicator was tested with a competitive
  rank-AUC against random gene sets matched on expression and detection (9,999 draws). The
  variance was inflated by a CAMERA factor estimated from residual correlations after removing
  specimen × segment means. BH was applied within each class. Robustness was checked with 5-bin
  matching and excluding the lowest expression tertile. Enriched pathways were grouped by gene
  overlap among class members.

## Coordinator check: is the mouse-only excess a male-mouse effect?

All specimens are male, and male-biased mouse PT programs concentrate in S2/S3 (Xiong et al. 2023;
Chen et al. 2023). Using the Xiong et al. per-segment female-versus-male tables (|log2FC| ≥ 1,
padj ≤ 0.05 in any PT segment), among probe-balanced genes:

| Class | Genes | Male-biased in mouse PT | Female-biased |
|---|---:|---:|---:|
| Conserved | 155 | 15.5% | 14.2% |
| Reversal | 52 | 5.8% | 9.6% |
| Mouse-only | 247 | 5.3% | 14.2% |
| Human-only | 66 | 6.1% | 4.5% |

199 of the 247 mouse-only genes are not sex-biased in either direction. The mouse-only excess is
therefore not driven by male-specific mouse programs. Correction: notebook 40 confirms mouse calls
in the 12 **male** mouse snRNA donors only (matching our male mice), not in a pooled-sex atlas.

Female-donor confirmation is now a sensitivity check (notebook 40 v2, A6 and A7). It keeps most
calls:

| Female confirmation | Mouse-only | Conserved | Reversal |
|---|---|---|---|
| Female snRNA | 214 of 247 | 151 of 155 | 44 of 52 |
| Female microdissection | 149 of 247 | 127 of 155 | 30 of 52 |

Mouse-only still exceeds human-only: 222 against 67, and 167 against 69.
