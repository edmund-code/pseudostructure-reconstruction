# Is human PT zonation globally weaker than mouse? (notebook 39, draft results text)

Notebook `analysis/notebooks/39_pt_zonation_amplitude.ipynb`, logic version
`39.zonation_amplitude.1`, results in `results/pt_zonation_amplitude/`. The protocol and decision
rules are pre-registered at the top of the notebook. Every contrast uses segment labels, so the
results do not depend on the reconstructed coordinate.

## Pre-registered results

### Statistics: the dilution-free ratio

We regressed both of our gradients on independent references (S2 − S1, OLS, bootstrap over
donors then genes). The human ÷ mouse amplitude ratio depends on which species the reference
comes from.

| Reference axis | Human ÷ mouse ratio | 95% CI |
|---|---|---|
| Mouse microdissection, male (Chen 2021) | 0.15 | 0.11–0.19 |
| Mouse microdissection, female | 0.23 | — |
| Mouse snRNA, male | 0.15 | — |
| Mouse snRNA, female | 0.16 | — |
| **Lake healthy human cortex** | **1.96** | **1.41–3.05** |

- The naive slope of our human gradient on our mouse gradient is 0.07, so dilution made the
  original "flattening" look stronger than any reference-based estimate.
- Deming and orthogonal regression give 0.15–0.23.
- S3c − early (cortical-like mouse S3): 0.28–0.33 on mouse axes, about 3.1 on human axes
  (Lake, Census human).

### Decisions

| Rule | Outcome |
|---|---|
| D1 | Met as written: two mouse reference groups give ratios below 0.7 with CIs excluding 1. Human-axis references give 2–3. |
| D2 (not log compression) | Met. Top expression tertile 0.18; probe-balanced genes 0.16. |
| D3 (independent human data, on a mouse axis) | Met: 0.14 (0.09–0.20). |
| D4 (S2-specific) | Not met. |
| D5 (tissue state explains it) | **Not met by any of its three routes.** Mouse AKI ÷ control is 0.74 for S2 − S1 (0.42 for S3 − early), well above 0.15. Lake donor injury score against mouse-axis amplitude gives ρ −0.24. The low-injury half of our structures gives the same ratio as the high-injury half (0.18 against 0.19 on microdissection). |

## Post-hoc diagnostics (not pre-registered; they explain the bracketing)

| Diagnostic | Value |
|---|---|
| Symmetric external ratio, Lake ÷ mouse snRNA on a human axis (Census, S3 − early) | 1.66 (1.29–2.21) |
| Same ratio on a mouse axis | 0.14 |
| Reference-free noise-corrected SD of S2 − S1 gradients, human ÷ mouse, ours | 0.58 (0.55–0.61) |
| Same, independent data (Lake healthy ÷ mouse snRNA) | 0.47 (0.43–0.50) |
| Same, S3c − early, ours | 0.84 |
| Same, mouse AKI ÷ control | 0.89 (S2 − S1), 0.71 (S3 − early) |
| Disattenuated cross-species correlation of S2 − S1 gradients, ours | 0.14 (0.10–0.19) |
| Same, external | 0.23 |
| Same, AKI vs control | 0.83 |
| Lake donors: injury score against human-axis amplitude | ρ = −0.41 |
| Lake human-axis amplitude, median: healthy / CKD / AKI | 0.69 / 0.50 / 0.35 |

## Interpretation

The "global flattening" seen in notebook 37 (interaction ≈ −0.8 × mouse gradient) is mostly a
mouse-centric view. Two things are true at once:

1. **Human and mouse PT zonate largely different genes.** The cross-species correlation of
   S1 → S2 gradients is about 0.14–0.23 after correcting for noise. Mouse AKI keeps 0.83, so injury
   does not scramble zonation like this.
2. **There is also a moderate global reduction.** The overall S1 – S2 contrast is about half as
   large in human (0.58 ours, 0.47 independent). The S3 − early contrast is not reduced in the
   same way (0.84).

Injury lowers amplitude in both species (mouse AKI; Lake AKI and CKD donors). It does not produce
a mouse-like loss of correlation, and healthy human donors remain flatter than healthy mice.

## Consequence for R3

- With the mouse-axis ratio r = 0.15, 32 of the 36 robust pathways are consistent with "mouse
  zonation, flat human". Four go beyond it (BH ≤ 0.10): Estrogen Response Early, Estrogen Response
  Late, Glutathione metabolism and Formation of Cornified Envelope.
- Lead genes beyond flattening:
  - **Reversals:** Acaa2, Dcxr, Ugt3a1, Cyp7b1 (S2 − S1) and Nt5e (S3c − early).
  - **Human-stronger or other:** Slc22a7, Slc5a8, Slc22a8, Psat1, Cyp24a1, Hmgcr, Pah, Igfbp4,
    Acox2, Acadl, Ephx1, Hadh, Lpl.
  - **Expressed only in human:** Rbp4, Aox1, Phgdh. These are zonated in human and cannot be
    tested against the mouse.
- Classic lead programs are consistent with mouse-specific zonation: Gatm, Acsm3, Crot, Nudt19,
  Slc27a2, Gclc/Gclm, Cyp51, Slc22a6, Slc13a3.

## Proposed R3 wording

> Human and mouse PT differ in **which** genes vary along the tubule more than in how strongly.
> In both our data and independent single-nucleus atlases, S1-to-S2 gene gradients are only weakly
> correlated between species (noise-corrected r ≈ 0.14–0.23), and the overall S1–S2 contrast is
> about half as large in human cortex (0.47–0.58). Most position-dependent pathways (32 of 36)
> reflect mouse zonation programs, such as fatty-acid oxidation, peroxisomal, creatine and sterol
> genes, that are nearly flat in human PT. A smaller set of genes and pathways carries human
> positional structure of its own: direction reversals (Acaa2, Dcxr, Ugt3a1, Cyp7b1, Nt5e) and
> genes zonated in human but absent or flat in mouse (Rbp4, Phgdh, Aox1, Slc5a8, Slc22a7).
> Injury lowers zonation amplitude in both species but does not explain the human pattern: mouse
> AKI preserves gene-level zonation (r = 0.83), and healthy human donors remain flatter than
> healthy mice.

**Caveats.**
- One human donor in our data.
- Segment contrasts depend on how each dataset annotates S1, S2 and S3.
- The human cortex lacks the outer stripe.
- The mouse AKI comparison uses the mouse-only pipeline's own counts and labels.
- The post-hoc diagnostics were not pre-registered.
