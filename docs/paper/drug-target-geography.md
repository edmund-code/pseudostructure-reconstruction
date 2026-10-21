# Drug-target geography along the PT, and an axis-free way to measure it

2026-10-21. Four drug-target tracks (glucose/SGLT2i, urate/uricosurics, OAT–OCT secretion, endocytic
uptake with phosphate and acid-base), each attacked; an audit and independent re-estimation of the
sharper-in-human set; then three independent figure designs, a judge, and sequential build → inspect
→ refine rounds in which the critic opened the rendered images and looked at them.

## 1 · The methodological unlock: a contrast with no coordinate in it

**SLC5A2/SGLT2 is an axis-construction gene in both species** (`fold-protected` False). Its position
on the reconstructed axis is therefore partly circular, which is why the earlier SGLT2 work fell back
to label-level detection and never ran a human arm at all.

The way round it is to drop the coordinate entirely: `log2((S3 cp10k + 0.1)/(S1 cp10k + 0.1))` from
the committed per-specimen `segment_means`, i.e. **from the H&E segment labels alone** — no
coordinate, no spline, no registration, no purity bins. I verified the five cache copies per
specimen are byte-identical and rebuilt the statistic myself for 10,669 of the 11,055 shared genes.

**Within-arm nulls on the axis-free statistic**, |d| > 0.3 in both members, over all genes:

| | sign disagreements |
|---|---|
| mouse vs mouse | **0 of 2,946** |
| human section vs human section | **14 of 2,664** (0.53%) |
| cross-species | **380 of 915** (41.5%) |

The round reported 0 of 951 and 1 of 485 on a 1,401-gene clean set; over the full set the human null
is 14, not 1. It is still a ratio of roughly 79× against the cross-species rate, and it is the
honest version — the human arm's null is never quite empty, which is consistent with everything else
this cohort has shown about it.

**This matters more than the drug results.** It means the reversal phenomenon does not need the
reconstructed axis at all: the same genes reverse on a statistic built from segment labels and raw
counts. The axis is how we *found* them and how we quantify shape; it is not what makes them visible.

## 2 · The drug map

Axis-free log2(S3/S1), both specimens of each arm agreeing in sign on every row I quote. Verified by
me from the caches:

| | mouse | human | call |
|---|---|---|---|
| SLC5A2 (SGLT2, gliflozins) | −6.12 / −7.00 | −3.14 / −3.70 | **conserved** |
| SLC2A2 (GLUT2) | −3.57 / −3.77 | −1.60 / −1.06 | conserved, exit co-located with entry in both |
| SLC5A12 (SMCT2) | −5.67 / −6.38 | −0.05 / −0.42 | **zonation lost in human** |
| GPT2 | −2.30 / −2.37 | −0.84 / −1.13 | conserved |
| ALDOB | −0.85 / −0.90 | −0.33 / −0.44 | conserved |
| FBP1 | −0.59 / −0.68 | −0.28 / −0.20 | conserved |
| **PCK1** | **−0.73 / −0.64** | **+1.21 / +1.43** | **reversed** |
| G6pc | −3.19 / −3.45 | *no human partner reaches the shared space* | see §4 |
| DCXR | −1.08 / −1.15 | +1.73 / +2.01 | reversed (a validation — see the novelty ledger) |

**The target sits in the same place in both species.** SLC5A2's human contrast sits at percentile 91
of the 45 clean genes with the steepest mouse gradients. Absolute human SGLT2 reaches further down
the tubule, but that is exactly the direction of the cohort-wide flattening, so **the honest claim is
conservation**: a drug targeting the convoluted segment in mouse targets it in human.

**The consequence does not.** Mouse gluconeogenesis is a single proximal block co-located with
SGLT2 — every step S1-weighted (G6pc −3.3, Gpt2 −2.3, Pcx −1.5, Got1 −1.1, Got2 −0.9, Aldob −0.9,
Pck1 −0.7, Fbp1 −0.6). In human the committed step **PCK1 inverts** while **FBP1, ALDOB, GPT2, GOT1
and GOT2 stay proximal**. This is a dissociation *inside one pathway, in one assay*, which is why
neither a global flattening nor a registration error can produce it: flattening shrinks everything
toward zero, which is what the trunk does, while PCK1 crosses zero.

**A second dissociation on the substrate side.** SLC5A12/SMCT2 (Na-coupled lactate, pyruvate and
β-hydroxybutyrate uptake) is mouse S1-restricted at 100–250× and in human is flat at a high level.
So the human PT can take up monocarboxylates along its whole length where mouse confines it to the
convoluted segment. (Mouse runs the reciprocal SMCT1/Slc5a8 distally, but human SLC5A8 detection is
**0.071**, so that half is not measurable — the segment-mean contrast still returns a number for it
and must not be used.)

## 3 · What is not measurable, which is half the result

Three walls, each with the number:

- **SGLT2 → SGLT1 is off the table.** Slc5a1 human detection 0.063. Also gone: GLUT1/SLC2A1
  (0.108/0.098, neighbour ratios 5.67/9.52), so "where does glucose leave the S3 cell" is
  unanswerable; and the entire Na/K-ATPase (Atp1a1, Atp1b1, Atp1b3, Fxyd2 all neighbour-flagged), so
  "is the pump positioned with the transporter" is unanswerable. Glycolysis is effectively unmeasured.
- **The aminoglycoside-uptake machinery has no resolvable gradient in either arm** — LRP2, CUBN, AMN
  and the PDZK1 scaffold all sit on zero. That is an honest negative and it should be published as
  one, because endocytic uptake geography is exactly what a nephrotoxicity reader will ask about.
- **OXCT1's human limb is not measured.** Human detection **0.157** against mouse 0.960, committed
  neighbour-marker flag True (ratio 5.03), raw human cp10k 0.18/0.28/0.95 — a thick-limb gene rising
  at the PT/TAL boundary. The mouse gradient is solid; the human one must be reported as absent, not
  as flat. The same mechanism demotes **SLC9A3/NHE3**: do not carry an NHE3 reversal.

**The distal-admixture asymmetry, quantified** — and this independently confirms the admixture
finding from the sharpening set, by a different route. S3/S1 ratio of thick-limb and distal markers
inside the PT-labelled units, mouse against human:

| | mouse A2 / A4 | human COR1 / MED1 |
|---|---|---|
| Umod | 2.75 / 4.58 | **14.42 / 16.37** |
| Slc12a1 | 2.87 / 3.93 | **11.88 / 13.60** |
| Cldn16 | 1.49 / 1.62 | 4.95 / 4.45 |
| Egf | 1.56 / 2.39 | 4.81 / 4.75 |
| Kcnj1 | 1.08 / 1.12 | 3.30 / 3.51 |
| Aqp2 | 1.23 / 1.04 | 2.66 / 3.55 |

The human distal-labelled compartment carries three to five times the distal-nephron signal that the
mouse one does. Any human gene whose apparent gradient rises distally must clear this, and the round
measured the dependence directly: the human-looks-more-distal shift rises with a gene's atlas
neighbour ratio (median +0.663 and 68.6% shifting distal at ratio ≥ 4, against −0.022 and 48.4% at
1–2). PCK1 (0.19), DCXR (0.78) and SLC5A12 (0.10) are all far below that range, which is why they
survive.

## 4 · A repo bug worth fixing: the terminal gluconeogenic step is invisible

Mouse **G6pc** and human **G6PC1** are the same 1:1 ortholog pair in the same inventory component
(6169), both well detected (mouse 0.555/0.608, human 0.987/0.952) — and **G6pc is absent from the
11,055-gene shared space** under either symbol. So glucose-6-phosphatase, the terminal and
rate-setting step of renal gluconeogenesis and the one that matters most to the SGLT2i question,
cannot appear in any cross-species table this project has produced. Mouse G6pc is also the steepest
gluconeogenic gene we have (−3.2/−3.5). Worth repairing before the glucose story is written.

Same class, not repairable: SLC9A3R1/NHERF1, LDHA, ALDOA and TPI1 are mouse-only on the panel — a
true probe-panel asymmetry.

Also to fix: nb92's `h1_gene_shapes.csv` labels Slc5a12 "SMCT1". SLC5A8 is SMCT1; SLC5A12 is SMCT2.

## 5 · The figures, and what their build log is worth

29 panels under the round's figure directory, built by three independent designs through a judge and
two refine rounds in which the critic opened the rendered images. The build logged **13 places where
the data disagreed with the specification and the builders drew the data**. Several are corrections
to numbers this project has been carrying:

- **SLC2A9 fails detection in the MOUSE arm (0.082); its human detection is 0.802** — the project
  record had this the wrong way round. Likewise **ABCC4** fails in mouse (0.089, human 0.533) and
  **SLC5A8** fails in human (0.071, mouse 0.851).
- **Fold protection:** the unprotected rows among the drug panel are SLC22A12, SLC22A7, SLC22A8 and
  SLC5A2; SLC34A1 and SLC34A3 *are* protected, so four of the eight conserved rows are protected,
  not two.
- **The smearing ladder:** w = 0.25 gives a median ratio of 0.662 and the first width reaching the
  observed 0.453 is w = 0.50, where 0.00% of genes are sharper against 13.09% observed. This matches
  the committed §F2a and corrects a relayed misattribution.
- **The within-species null denominators (85 / 93 / 87) reproduce only when the per-specimen trend ρ
  is taken from the fitted specimen curves**; the stored `rho <specimen>` columns give 81 / 85 / 87
  under the same gate. Quote the gate with the number.
- **CHP1 does not satisfy the "both human marks right of both mouse marks" rule** on the curve
  amplitude the panel plots; its membership comes from the model-free ratio.
- **In the injury design Ctrl1A2 is not clean either** (f_prox 0.181 against f_dist 0.114), not just
  IR2A2.

**A count of mine to update.** The retraction figure classifies the 29 centre-of-mass movers on the
**pure-compartment** contrast per specimen and gets 18 opposite-sign, 6 same-sign more than 1.5×
apart, 4 same-sign within 1.5×, 1 with specimens disagreeing inside an arm. I had committed 16
opposite-sign, which is the same count taken on the **curve-derived** arm contrasts. Both are
correct about their own statistic; they agree on 25 of 29. **Report 18 on the retained statistic**,
since the curve-derived one is the retracted one. The four disagreeing genes are Hao2, Chp1, Msrb1
and Txnip — Hao2 being the gene hypothesis round 2 independently removed as a cell-identity artefact,
and Txnip the known within-arm flip.
