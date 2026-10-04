# Choosing the PT pathway analysis for the paper

Notebook `analysis/notebooks/31_pt_pathway_method_selection.ipynb`. Run 2026-10-04,
logic version `31.method_selection.1`. Inputs: notebook 13's PT scFates coordinate and
notebook 12's gene universe, pathway families and sensitivity runs. Nothing upstream is refit.

Cohort: 2 control mice, 2 cortex sections from 1 human donor, 12,272 PT structures,
11,071 genes, 1,513 pathways. All four specimens are male (Y-linked genes detected in every
one), so sex does not separate the species. Everything below is descriptive.

## The control this adds

Notebook 12 compares screens to each other. It has no measurement of how many pathways a
screen reports when there is **no** species difference to find. With two specimens per
species, the four specimens can be relabeled into two balanced mouse+human groups. Both
relabelings keep the structures, coordinate, weights, gene universe, pathway family and test
unchanged; only the species difference is removed. The species' own smooth shape enters each
relabeled model as a nuisance term, so the tested term is a difference between specimens of
the same species. This is the random-partition control of
[Lamian (Hou et al. 2023)](https://doi.org/10.1038/s41467-023-42841-y).

The **negative-control discovery ratio** (NCDR) is mean relabeled discoveries ÷ species
discoveries. A strategy is called specific when NCDR < 0.25 and the relabeled groupings call
fewer than 5% of tested pathways. Two relabelings screen for gross nonspecificity; this is not
a calibrated false-discovery rate.

## Result: only smooth position-dependent screens are specific

| Strategy | Species | Relabeled (mean) | NCDR | Specific |
|---|---:|---:|---:|:--:|
| 1 Whole-PT pseudobulk + signed GSEA | 158 | 224 | 1.42 | no |
| 2 S1/S2/S3 pseudobulk + signed GSEA | 204 | 299 | 1.47 | no |
| 3a Pathway score GAM, offset | 1444 | 1152 | 0.80 | no |
| 3b Pathway score GAM, shape | 1509 | 281 | 0.19 | no (19% of pathways) |
| 4 Gene GAM offset `T_level` | 27 | 104 | 3.85 | no |
| 5 Gene GAM total `T_total` | 76 | 131 | 1.72 | no |
| 6 Gene GAM × S1/S2/S3 steps | 67 | 42 | 0.62 | no |
| 7 Gene GAM × 6 equal bins | 100 | 42 | 0.42 | no |
| **8 Gene GAM × smooth position `T_spatial`** | **116** | **15** | **0.13** | **yes** |
| 8b `T_spatial`, size-only random sets | 344 | 57 | 0.17 | yes |
| 8c `T_spatial`, preranked GSEA | 311 | 29 | 0.09 | yes |
| Smooth position beyond S1/S2/S3 steps | 176 | 37 | 0.21 | yes |

The gene statistics behave correctly: whole-PT DESeq2 calls 9,010 genes at padj ≤ 0.05 for
species and 2–8 for the relabelings. The nonspecificity is in the pathway layer. A competitive
set test asks only whether member genes **rank** above other genes, and specimen-level
variation moves whole pathways together, which gene-permutation tests assume away
([Goeman & Bühlmann 2007](https://doi.org/10.1093/bioinformatics/btm051);
[Wu & Smyth 2012](https://doi.org/10.1093/nar/gks461)). With two specimens per species a
pathway-level average difference cannot be separated from specimen variation. Structure-level
tests of pathway scores additionally treat 12,272 structures as replicates
([Squair et al. 2021](https://doi.org/10.1038/s41467-021-25960-2)).

Within-species controls agree: comparing the two mice gives 48 `T_level` pathways and 16
`T_spatial`; the two human sections give 94 and 0. Specimen differences are offsets, not
positional remodeling.

**Decision.** Primary analysis = gene-level `T_spatial` from the nested GAM on the scFates
coordinate, tested with notebook 12's covariate-matched rank-AUC. It keeps notebook 12's list
unchanged, is the only strategy passing specificity, answers the positional question, and
reads positions back from D(s)/S(s). Constant-offset results (whole-PT, S1/S2/S3, `T_level`)
are reported as context only: they may be real biology, but this cohort cannot separate them
from specimen variation. Conventional analyses are not wrong; they ask a question this design
cannot answer specifically.

## Final list

129 notebook-12 spatial candidates → 104 not specimen-sensitive → 84 stable (retained in ≥5 of
7 planned sensitivity runs) → **82 correlation-supported**. Versus conventional screens:
54 new with pseudospace (no conventional species call), 26 whose conventional call is also
produced by a relabeling, 1 conventional call that no relabeling produces, 1 not testable
conventionally.

Themes (keyword rules applied after the list was frozen): amino-acid metabolism (12),
xenobiotic & drug metabolism (12), sterol/steroid/bile acid (12), vitamins & cofactors (9),
fatty-acid oxidation & peroxisome (6), glucose & energy metabolism (6), solute transport (6),
other/signaling (19).

70 of 11,071 genes (0.6%) reverse species direction along the PT under the display rule
(|Z| ≥ 3 and |difference| ≥ 0.25 over ≥3 contiguous grid points on each side). A whole-PT
comparison reports one direction per gene and cannot represent these.

## Caveats to keep beside any number here

- Two mice and two sections from one human donor. Specimens, not structures, are the
  replication units; no result is population-level human inference.
- The relabeling control uses two partitions only.
- Themes, drivers and reversal counts are display rules on working statistics, not tests.
- Both human slices are healthy cortex despite the cortex/medulla labels.

Outputs: `results/pt_pathway_method_selection/` (tables, `figures/`, `run_manifest.json`).
