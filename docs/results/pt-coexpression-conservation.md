# Co-expression conservation between human and mouse PT (notebook 97)

**The question.** Every species comparison so far compares levels or shapes. Levels are confounded by the
species-specific probe panels. Shapes gave about a quarter of each gradient shared and 74 reversed genes, but no
pathway-level species program (notebooks 90 D, 91 S).

A third, level-free readout asks: does a gene keep the same co-expression partners in both species? This is a
different question from whether its gradient matches. Co-expression is computed within a specimen, so a per-gene probe
efficiency cancels.

- Is per-gene co-expression conservation reproducible?
- Is it confirmed in independent atlases?
- Does it add to the gradient comparison?
- Does it find pathway-level rewiring that gradients did not?

**The measure** is Crow, Suresh, Lee & Gillis 2022 (NAR 50:4302, PMID 35451481, read in full; with CoCoCoNet, Lee et
al. 2020, PMID 32392296):
- Spearman co-expression;
- for gene g, its top 10 co-expression partners among the shared orthologs in one network;
- the AUROC with which g's co-expression ranks in the other network recover them;
- the average of both directions.

The AUROC is 0.5 by chance and 1 when the neighbourhoods are identical.

Protocol frozen on 2026-10-18, before any co-expression was computed on these data.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/97_pt_coexpression_conservation.ipynb`; logic `97.coexpression.1`; seed 97.
- Results `results/pt_coexpression_conservation/`.
- Module `pseudospace/coexpression_conservation.py`, tests `tests/test_coexpression_conservation.py`.
- **Our units:** the human vs mouse export through `hvm_inputs` (notebook 93's construction, mt genes excluded). It
  gives counts, positions in [0, 1], specimens (Ctrl1A2, Ctrl1A4, HUK1_COR1, HUK1_MED1) and the 71 universe.
  Expression is log1p(count / library × 10⁴), with the library over the measured genes.
- **Atlases:** `data/external/census_pt_cells/{human,mouse}_pt_cells.h5ad`. These are nucleus counts over the accepted
  ortholog panel in mouse symbols, with donors. Human: renal-cortex benchmarking dataset. Mouse: the sex-specific
  kidney dataset.
- **Our earlier readouts, for comparison only:** notebook 90's end contrasts and species categories
  (`d_end_contrasts.csv.gz`) and notebook 71's human vs mouse library sets.

## 2 · Genes and networks

- **Universe G:** genes of the 71 human vs mouse universe that meet both rules:
  - detected (count > 0) in ≥ 10% of units in each of our four specimens;
  - detected in ≥ 5% of nuclei in each atlas species.
- **Total network,** per specimen: the Spearman correlation over that specimen's units, on G.
- **Residual network,** per specimen: the same on residuals. A residual is the unit's expression minus the mean of
  its specimen's units in the same position bin (20 equal-count bins of x). This network asks what co-varies between
  tubules beyond their position.
- **Atlas networks:**
  - donors with ≥ 300 PT nuclei;
  - per donor, the Spearman correlation over its nuclei;
  - per species, the aggregate is the mean of the donors' rank-standardised networks (Crow's aggregation; edges
    ranked over the whole network and divided by the maximum rank);
  - **atlas halves:** each species' donors split in two (alternating, ordered by nucleus count), each half aggregated.

## 3 · Per-gene scores (both layers)

- **Within-species:** W_m = AUROC(Ctrl1A2, Ctrl1A4) and W_h = AUROC(HUK1_COR1, HUK1_MED1). Human is one donor, so W_h
  is across sections.
- **Cross-species:** X = the mean AUROC over the four mouse × human specimen pairs.
- **Rewiring:** D = (W_m + W_h) / 2 − X. Positive means the gene's partners agree less between species than within.
- **Atlas:**
  - X_atlas = AUROC(mouse aggregate, human aggregate);
  - W_atlas = the mean of the two species' half-against-half AUROC;
  - D_atlas = W_atlas − X_atlas.

## 4 · Questions and decision rules (fixed now)

**Q1 · Is there a reproducible neighbourhood?**
- The median W_m and the median W_h, total and residual layers.
- The measure is informative if median W_m ≥ 0.7 on the total layer. 0.7 is "quite good" in CoCoCoNet's reading of
  AUROC.

**Q2 · Is species rewiring confirmed outside our data?**
- Spearman ρ over G between D and D_atlas, and between X and X_atlas.
- **External confirmation** needs ρ(D, D_atlas) ≥ 0.2 with a permutation p ≤ 0.01. The null permutes the atlas gene
  labels 1,000 times.

**Q3 · Does it add to the gradient comparison?**
- Notebook 90's level-free gradient divergence over the shared genes: |E mouse − E human|, the absolute difference of
  the species-mean end contrasts (log2). *Amended before any computation: the frozen text said "1 minus the cosine",
  which is undefined for one end contrast per gene.*
- ρ between D and that gradient divergence.
- D by 90's categories (conserved, reversed, mouse only, human only, none), with a Mann–Whitney test of reversed
  against conserved.
- **Non-redundant** if |ρ| ≤ 0.5 while Q2 holds.

**Q4 · Pathway-level rewiring (91 found none from gradients).**
- cameraPR (inter-gene correlation 0.01) on D over 71's human vs mouse library sets, restricted to G (10–300 genes),
  in our data and in the atlas.
- **A rewired set** has q ≤ 0.05 and D "Up" in both. A **conserved** set has q ≤ 0.05 and "Down" in both.
- The count of replicated sets is compared with the count expected by chance, the product of the two call rates times
  the number of sets.

**Q5 · Tubule-level programs beyond position (residual layer, descriptive).**
- Leiden communities (resolution 1, seed 97) on the mean residual network of the mouse specimens: edges = the top 10
  partners per gene, symmetrised.
- For each module: its size, its mean X (cross-species conservation of its members' neighbourhoods), and its top
  library sets by ORA.

**Method verdict.** Co-expression conservation counts as a useful species method if Q1 and Q2 hold, and either Q3
shows it is non-redundant or Q4 replicates more sets than chance.

## 5 · Outputs and figures

- **Tables:**
  - `genes.csv`: per gene W_m, W_h, X, D, atlas scores and the 90 category, both layers;
  - `networks_summary.csv`, `q1_reproducibility.csv`, `q2_external.csv`, `q3_gradients.csv`, `q4_sets.csv`,
    `q5_modules.csv`;
  - `verdicts.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 97a: distributions of W_m, W_h and X, both layers;
  - 97b: D against D_atlas;
  - 97c: D by gradient category;
  - 97d: set z, ours against atlas.

## 6 · Caveats

- One human donor, so W_h is across two sections. Donor and species are confounded in X.
- Co-expression across tubules is dominated by position (zonation), so the total layer largely compares zonation
  partners. The residual layer removes position.
- Atlas nuclei are sparse, so their Spearman networks are noisy. Aggregation across donors is Crow's remedy.
- The atlases are other people's tissue, sex mix and platform. Agreement is a lower bound.

## Notebook 97 results

`analysis/notebooks/97_pt_coexpression_conservation.ipynb` (logic `97.coexpression.1`) writes to
`results/pt_coexpression_conservation/`. It followed the protocol frozen in 137fd92, amended before any computation in
590b6c5. G1 passed: the export fingerprint is notebook 70's, and `hvm_inputs` rebuilds 71's 11,055-gene universe.
Runtime is 3 min, almost all of it the atlas networks. Module `pseudospace/coexpression_conservation.py` with
`tests/test_coexpression_conservation.py` (3 tests).

**How the protocol was read.**
- Our expression is `hvm_inputs`' Y, and atlas expression is log1p CP10k on the per-nucleus total.
- The atlas halves alternate donors.
- The residual layer has no atlas twin, so Q2 is judged on the total layer.
- Q5 settings are as listed.
- One post hoc section (5.1).
- Figure 97c was redrawn from `genes.csv` after execution to rotate overlapping tick labels. The notebook source
  carries the fix.

**Inputs.**
- Units: Ctrl1A2 2,436, Ctrl1A4 2,557, HUK1_COR1 3,073, HUK1_MED1 3,317.
- Atlas donors with ≥ 300 nuclei: mouse 28 (59,277 nuclei), human 7 (15,059).
- G holds 3,679 genes. 6,019 pass our detection rule and 4,612 the atlases'.

**Per-gene scores (medians).**

| Layer | W_m | W_h | X | D | share D > 0 |
|---|---|---|---|---|---|
| total | 0.670 | 0.675 | 0.488 | 0.202 | 0.98 |
| residual | 0.668 | 0.622 | 0.515 | 0.133 | 0.88 |
| atlas | 0.969 | 0.887 | 0.591 | 0.292 | 0.99 |

W is bimodal in our data: on the total layer, the quartiles of W_m are 0.595 and 0.962.

**Questions.**

| Question | Value | Holds |
|---|---|---|
| Q1 reproducible neighbourhoods (median W_m ≥ 0.7, total) | 0.670 | **no** |
| Q2 rewiring confirmed in the atlases (ρ(D, D_atlas) ≥ 0.2, p ≤ 0.01) | ρ 0.297 (null 99% 0.037), p 0.001 | yes |
| Q3 non-redundant with gradient divergence (\|ρ\| ≤ 0.5) | 0.377 | yes |
| Q4 more replicated sets than chance | 1 rewired (Xenobiotic metabolism) and 1 conserved (mRNA processing), against 0.04 and 0.01 expected; the atlas calls only 1 set each way | yes, formally |

- **Method verdict (frozen): "not shown useful".** Q1 fails narrowly, because the median gene's neighbourhood
  reproduces at 0.67 between two mice.
- **Other results:**
  - Cross-species X is 0.49 in ours against 0.59 in the atlases.
  - Residual co-expression (beyond position) reproduces as well as total within species (0.668), but crosses species
    at near chance (0.515), and ρ(X residual, X atlas) is only 0.066.
- **Rewiring agrees with gradient reversal.** Median D by notebook 90's categories:
  - reversed 0.761 (31 genes);
  - human only 0.397;
  - mouse only 0.304;
  - conserved 0.213;
  - none 0.193.

  Reversed against conserved: Mann–Whitney p < 10⁻¹⁰.
- **Post hoc, genes reproducible in both species** (W_m and W_h ≥ 0.7; 868 genes):
  - ρ(D, D_atlas) is 0.342 against a null of 0.080 (99%);
  - ρ with gradient divergence is 0.278.
  - **Rewired in both**, top decile of D in ours and in the atlases (25 genes): Dst, *Mme*, Gclm, *Sh3gl2*, *Lap3*,
    Mgll, Plekha5, Napsa, *Hpgd*, Mfap3l, Atp6v1b2, Glg1, Aplp2, Ctsh, Slc16a9, Agpat3, Slc41a2, Ccdc88c, Larp4b,
    *Slc17a3*, Acot12, Acox2, Ndfip2, *Ugt3a1*, Celsr1. Italics mark genes on notebook 90's independent list of
    reversed gradients.
  - **Conserved in both** (16 genes): Sptan1, Fam151a, Slc4a4, Itpr2, Apoe, Tpp1, Cotl1, Dpep1, Iyd, Gpx3, Tbx2, Lgmn,
    Slc6a19, Slc7a8, Atp6v0c, Gnptg. These are core PT transport and handling genes; Slc7a8 is also a 90 "conserved"
    gradient.
- **Q5, residual modules.** The fatty-acid and xenobiotic-metabolism module holds 489 genes (ORA fdr ~10⁻¹⁹).
  *(Corrected 2026-10-20: it is the **smallest** of the five modules, not the largest. `q5_modules.csv` sizes are
  1,079 / 839 / 748 / 524 / 489; the largest is module 2, whose top sets are Protein Secretion, Signaling By ROBO
  Receptors and Cellular Response To Hypoxia at fdr 0.11. The 489 size, the annotation and the 0.53 conservation
  are all correct.)*
  Tubules at the same position co-vary in their lipid and drug-metabolism programs, a tubule-level program beyond
  position. Its cross-species conservation is near chance (0.53).

**What this answers.**
1. **As a genome-wide method, co-expression conservation is not shown useful** by its frozen rule: the median gene's
   neighbourhood reproduces between two mice at 0.67, below 0.7.
2. **Where it is measurable, it confirms and extends the gradient comparison.**
   - Species rewiring of a gene's partners replicates in independent atlases (ρ 0.30, and 0.34 for reproducible
     genes).
   - It is only moderately related to gradient divergence (0.28–0.38).
   - It converges on the same genes: six of the 25 doubly rewired genes are reversed-gradient genes from notebook 90,
     found by a different, level-free readout.
3. **Pathway-level rewiring:** only xenobiotic metabolism replicates, consistent with notebooks 62 and 86's
   rodent-expanded drug-metabolism families. This confirms 91's finding that species divergence is gene-level, not
   pathway-level.
4. **Use in the paper:** a supporting, level-free cross-check of the species gene lists, not a discovery method of its
   own.
