# Revision: corrected R3 wording and numbers (notebook 43)

Source: `analysis/notebooks/43_pt_revision_classes.ipynb` (logic version `43.revision_classes.1`).
Results are in `results/pt_revision_classes/reviewed/`. The protocol was pre-registered in the
notebook header.

To rerun with workstream 1's labels, use
`--labels results/pt_revision_labels/segment_labels.csv`, with columns `structure` and `segment`.
If that file exists, the notebook picks it up automatically. Results then go to
`results/pt_revision_classes/<label file stem>/`.

## What changed

1. **The "largely different genes" claim survives the same-species ceiling (D2).**
   - Noise-corrected same-species, cross-dataset correlations:

     | Contrast | Within human | Within mouse |
     |---|---|---|
     | S2 − S1 | 0.70 (ours vs Lake) | 0.63–0.88 (0.88 ours vs snRNA male) |
     | S3 − early | 0.70–0.93 | 0.73–0.93 |

   - Cross-species correlations are 0.09–0.20 (S2 − S1) and 0.13–0.24 (S3 − early).
   - **Conservation index** (cross ÷ √(human ceiling × mouse ceiling)):
     - 0.17 (0.12–0.22) for S2 − S1;
     - 0.28 (0.22–0.33) for S3 − early;
     - 0.29 (0.23–0.35) for S3c − early;
     - external-only (Lake against mouse snRNA): 0.27–0.28.

     Every upper bound is below 0.5, so the pre-registered rule keeps "largely different genes".
   - **Rat is excluded.** Its replicate reliability is 0.03 (S2 − S1) and 0.45 (S3 − early), below
     the 0.5 rule; 80% of genes map by symbol. All rat-based arguments are withdrawn.
2. **The mouse-only excess does not survive symmetric classification (D3).**
   - Notebook 40's union classes gave 247 mouse-only against 66 human-only.
   - Matching human thresholds to the lower human amplitude (effect floor and flat margin × 0.58 /
     0.84 / 0.70) reverses this to 59 against 141.
   - The symmetric rule, which adds external confirmation of flat calls, gives 41 against 76.

   | Contrast (symmetric rule) | Mouse-only | Human-only |
   |---|---|---|
   | S2 − S1 | 31 | 58 |
   | S3c − early | 13 | 28 |
   | Full S3 − early | 1 | 110 |

   - Under the symmetric rule, conserved rises from 155 to 196 and reversal from 52 to 75. The
     female confirmations (A6, A7) give 38 against 71 and 31 against 78.
   - Strong mouse-only genes (|mouse| ≥ 1.5, |human| ≤ 0.3): 21 under v1, 1 under the symmetric
     rule.
   - Pre-registered decision: **no consistent asymmetry**. The earlier excess was a threshold effect
     of lower human amplitude.
3. **The pathway layer summarises mouse zonation (M7).**
   - Of the 36 robust pathways:
     - 26 are summaries of mouse zonation (≥ 50% of their zonated contributing members are
       mouse-led);
     - 3 are shared (Nitrogen metabolism, Aspirin ADME, Glutathione metabolism);
     - 2 are mixed (Drug ADME, mitochondrial FAO);
     - 5 are unresolved (fewer than 3 zonated contributing members);
     - none are human-led.
   - 8 robust pathways are probe-sensitive: Glutathione metabolism, Steroid hormone biosynthesis,
     Chemical carcinogenesis, Retinol metabolism, xenobiotics by CYP450, Aspirin ADME, Formation of
     Cornified Envelope, Biological Oxidations.
   - Of notebook 39's four "beyond flattening" pathways, Glutathione metabolism and Formation of
     Cornified Envelope fail the equal-probe check and are dropped. Estrogen Response Early and
     Late remain, but their members are mostly mouse-led.
   - External replication: robust pathways beat relabel-called pathways only marginally
     (Mann–Whitney p 0.016 and 0.033 with male-mouse references; 0.19 and 0.21 with female).
     Report this beside the uncalled comparison.

## Proposed R3 wording

> **Human and mouse PT zonate largely different genes, with a conserved core of strong markers.**
> Same-species comparisons across independent datasets set the ceiling: noise-corrected
> correlations of segment gradients reach 0.70–0.93 within human and 0.63–0.93 within mouse.
> Between species they are only 0.09–0.24. Expressed against these ceilings, the conservation of
> gene-level zonation is 0.17 (95% CI 0.12–0.22) for S1 → S2 and 0.28 (0.22–0.33) for S3 versus
> early PT. The same holds using the external atlases alone (0.27–0.28).
>
> At the level of individual genes, the strongest gradients are mostly shared. Under symmetric,
> amplitude-matched calls, 196 genes are zonated alike in both species and 75 in opposite
> directions. Confidently species-only genes are few, 41 in mouse and 76 in human, and most genes
> cannot be classified with two specimens per species. Species specificity is therefore a
> genome-wide property of the gradients, plus a lower overall human amplitude (noise-corrected SD
> ratio 0.58 for S1 → S2). It is not a large set of mouse-only genes; the earlier mouse-only excess
> was an artefact of applying mouse-scaled thresholds to the flatter human gradients.
>
> The position-dependent pathway screen mainly summarises mouse zonation that human PT does not
> share: 26 of the 36 robust pathways, including fatty-acid oxidation, peroxisome, bile acid,
> branched-chain amino-acid and vitamin metabolism. None is led by human-specific zonation.

## Corrections elsewhere

- **Introduction.** "confirmed every call in an independent atlas" → "confirmed every zonated call
  in an independent atlas of the same species; flat calls were additionally checked in the
  revision".
- **R3.** Report S3c − early as co-primary: conservation index 0.29, SD ratio 0.84.
- **R5 and Discussion.** Remove rat arguments: the sterol "rodent vs mouse" point and the rat
  support for Dcxr and Rbp4.
- **Abstract.** Replace "mouse-only exceeds human-only" (if present) with the conservation index
  and the conserved-core statement above.
