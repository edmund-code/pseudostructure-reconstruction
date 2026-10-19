# Does sequence evolution predict which PT gradients are conserved between human and mouse? (notebook 106)

**The question.** Notebook 90 split the 1:1 orthologs graded along the PT into four groups:
- conserved in direction (189 genes);
- reversed (74);
- graded in mouse only (229);
- graded in human only (70).

All of these are level-free end contrasts, so the species' probe panels cannot make them. Is a gradient reversal a
property of genes that evolve fast? Or does it strike constrained genes too, which would point to regulatory change?

**What the literature predicts.** The methods scout read these in full
(`docs/research/methods-scout-contamination-networks-variability-evolution.md`, section D):
- **Strand et al. 2007** (PLoS Genet), across human and mouse brain regions. The correlation of regional log ratios
  is about 0.53 in the lowest dN/dS quartile against about 0.33 in the highest.
- **Liao & Zhang 2006** (Mol Biol Evol) found **no** correlation between dN/dS and a shape-only divergence measure
  (1 − r across tissues). Our end contrast is a shape measure, so a null has precedent.
- **Warnefors & Kaessmann 2013** (Genome Biol Evol): expression divergence correlates with dN (ρ 0.27–0.38). Among
  lowly expressed genes the correlation disappears once expression level is controlled.
- **Huminiecki & Wolfe 2004; Kryuchkova-Mostacci & Robinson-Rechavi 2016:** duplication is the strongest categorical
  predictor of expression divergence.
- **Cardoso-Moreira 2019 and 2020** (Nature; Cell Rep): disease and essential genes are depleted for human–mouse
  expression differences.
- **No study relates the conservation of a zonation or axial pattern to dN/dS, duplication or constraint.**
  Yakubovsky 2026 compared human and mouse liver zonation but did not test sequence evolution.

Frozen on 2026-10-18, before any sequence or constraint statistic was joined to the gradient categories.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/106_pt_gradient_evolution.ipynb`; logic `106.gradient_evolution.1`; seed 106.
- Results `results/pt_gradient_evolution/`.
- **Gradients:** notebook 90's `d_end_contrasts.csv.gz`. It gives E for each of 4 specimens and the category of
  each of 7,971 1:1 genes.
- **Gene identifiers, expression and probes:** notebook 62's `gene_inventory_mouse.csv` and
  `gene_inventory_human.csv` (Ensembl gene id, `n_probes`, `detection_1/2`, `log2_cp10k`).
- **dN, dS:** `data/external/ensembl_compara_human_mouse/` (Ensembl 98, `ortholog_one2one` rows, matched on both
  gene ids). ω = dN/dS. Pairs are excluded when dS is 0 or more than twice the median dS (Strand's cap).
- **Paralogs:** `data/external/ensembl_compara_paralogs/` (Ensembl 98 within-species paralogs). A gene **has a close
  paralog** when it has a within-species paralog with ≥ 50% protein identity (`identity`) in human or in mouse.
- **Constraint:** `data/external/gnomad_v4_1_constraint/` (gnomAD v4.1). LOEUF (`lof.oe_ci.upper`) of the canonical
  transcript, matched on the human symbol.

## 2 · Tests (fixed now)

**Main comparison:** reversed (74) against conserved (189). Both groups are graded in all four specimens, so they are
matched for detectability and reliability.

**Strata:** quintiles of expression, the mean over species of `log2_cp10k`, over the 263 genes.

**P1 · Constraint on the protein (primary).**
- AUC of log ω for reversed against conserved (reversed higher means AUC > 0.5).
- Null: 2,000 permutations of the labels within expression quintiles; the p is one-sided.
- **Holds** when AUC ≥ 0.55 and p ≤ 0.05.

**P1b · Strand-style quartiles.**
- Genes: those graded in at least one species (the four categories).
- Within each ω quartile:
  - r = Pearson(E_mouse, E_human), with E the species mean of its two specimens;
  - reliability-corrected r_c = r / √(rel_m · rel_h), where rel is the Pearson correlation between a species' two
    specimens within that quartile.
- **Holds** when r_c(lowest ω quartile) − r_c(highest) ≥ 0.1, and the difference lies above the 95th percentile of
  2,000 permutations of ω among the genes within expression quintiles.

**P2 · Duplication.**
- Share with a close paralog, reversed against conserved.
- Cochran–Mantel–Haenszel odds ratio across expression quintiles.
- **Holds** when OR ≥ 1.2 with p ≤ 0.05.

**P3 · Loss-of-function constraint.**
- AUC of LOEUF for reversed against conserved (reversed less constrained means AUC > 0.5).
- **Holds** when AUC ≥ 0.55 with a permutation p ≤ 0.05, as in P1.

**Artifact control (after Sepp et al. 2024).**
- A random forest (500 trees, 5-fold stratified CV, seed 106) is trained on these features only:
  - mouse and human probe counts;
  - detection in each specimen;
  - log2 CP10k per species.
- If its cross-validated AUC for reversed against conserved is > 0.6, P1–P3 are **not read biologically**; panel or
  detection artefacts would be live.

**Descriptive.**
- P1 and P3 AUCs for mouse-only and human-only genes against conserved.
- The ω and LOEUF of the named reversed genes (Apoc3, Pdk4, Vnn1, Hpgd and others).

## 3 · Outputs and figures

- **Tables:** `genes.csv` (gene, category, ω, dN, dS, close paralog, LOEUF, covariates), `tests.csv`,
  `quartiles.csv`, `artifact_control.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 106a: ω and LOEUF distributions by category;
  - 106b: r_c per ω quartile, with the permutation band;
  - 106c: E_mouse against E_human coloured by ω quartile.

## 4 · Caveats

- 74 against 189 genes gives modest power (z ≈ 2.5 at AUC 0.6). P1b carries most of it.
- dN/dS comes from one mouse–human comparison: lineage-specific change cannot be placed on either branch.
- "Close paralog" at ≥ 50% identity is a proxy for recent duplication. Ensembl's TSV has no duplication node.

## Notebook 106 results

`analysis/notebooks/106_pt_gradient_evolution.ipynb` (logic `106.gradient_evolution.1`) writes to
`results/pt_gradient_evolution/`. It followed the protocol frozen in d4f04d4 and ran in 15 s.
- G1: notebook 90's categories.
- G2: dN/dS for 6,845 of 7,971 genes (dS cap 1.08); LOEUF for 7,710.

| Test | Result | Holds |
|---|---|---|
| P1 log dN/dS, reversed against conserved | AUC 0.51 (null 95% 0.55), p 0.29; 70 vs 158 genes | **no** |
| P1 sensitivity, 100 − % identity | AUC ≈ 0.5 | no |
| P1b reliability-corrected r(E_mouse, E_human) per dN/dS quartile, low to high | 0.33, 0.38, 0.56, 0.43; lowest − highest −0.09 (null 95% 0.17) | **no**: the most constrained quartile is the *least* conserved |
| P2 close paralog (≥ 50% identity) | 16% of reversed against 19% of conserved; CMH OR 0.78, p 0.50 | **no** |
| P3 LOEUF | AUC 0.49, p 0.41 | **no** |
| Artifact control: random forest on probe counts, detection, expression | cross-validated AUC 0.54 | readable (≤ 0.6) |

- Species-only genes against conserved give the same picture (AUC 0.44–0.48 for dN/dS and LOEUF).
- Median dN/dS is 0.10 for reversed genes and 0.11 for conserved ones.
- Several of the strongest reversals are in highly constrained proteins:
  - Pdk4 (ω 0.06, 93% identity);
  - Ddit4 (0.05);
  - Oat (0.07);
  - Socs2 (0.09);
  - Hpgd (0.13).
- Fast-evolving ones exist too: Apoc3 (0.46), Slc17a3 (0.44).

**What this answers.**
1. **Gradient reversal between human and mouse is not explained by protein evolution, duplication or
   loss-of-function constraint.**
   - Reversed genes are as constrained as conserved ones.
   - Strand's brain-region result does not carry over to a positional shape measure. Liao & Zhang's null for
     shape divergence does.
   - The reversals hit conserved proteins, which points to changes in regulation along the nephron. Examples are
     PPARα/fasting targets such as Pdk4, and stress-response genes such as Ddit4.
2. **The reversed category is not a panel artefact.** Probe counts, detection and expression cannot tell reversed
   from conserved genes (AUC 0.54). This is a new, independent check on notebook 90's species result.
3. **As a method:** sequence features add no prior for which gradients to trust across species. Do not weight
   species calls by conservation of the protein. For the paper, this is a negative result and a guard on the
   reversed-gradient list.
