# Novelty ledger: what is already published, and what is ours

2026-10-21. A full-text literature pass over every surviving human-versus-mouse claim, the
drug-target geography, and the reconstruction method itself. Seven reading groups, each followed by
an adversarial prior-art hunter forbidden from searching the way the reader had — by phenomenon
rather than gene name, in other organs, in the pre-2000 literature, and tracing review statements
back to their primary sources.

This replaces the abstract-level novelty classification recorded in
`hvm-continuous-discovery-round.md` §2, which the standing full-text rule had already flagged as a
weakness. The rule paid for itself immediately: a full-text read **deleted the citation behind our
own UGT3A1 framing**.

## 0 · Two reading rules, and one scope violation to declare

**Access labels are preserved exactly as the pass recorded them** — `FULL TEXT READ`,
`ABSTRACT ONLY`, `NOT ACCESSED`. A claim resting on an abstract is marked
*load-bearing-on-abstract* and may not be published without the full text.

**What I verified myself**, as opposed to relaying:

| Source | What I checked, in this session | Verdict |
|---|---|---|
| Mackenzie 2008, PMID 18981171 | Europe PMC record and abstract: *"UGT3A1 is found in the liver and kidney, and to a lesser, in the gastrointestinal tract."* A human cloning, expression and enzymology paper. | **Our citation is wrong.** It reports a tissue distribution, not a species difference. |
| Yakubovsky 2026, PMID 41986723 | Full text via the Europe PMC OA endpoint. Found: *"Lipid metabolism genes such as FASN, APOA5, GPAM, ACSL5 and DPP4 … were strongly pericentral in humans but periportal or weakly zonated in mice"*, with the IF legend *"the non-zonated expression in mice"*. And: *"Correlation between portal biases in mice and humans, defined as log2-transformed expression ratios between periportal and pericentral zones, was unexpectedly low (R = 0.15)"*. | Confirms both the DPP4 scoop **and** that the paper compares **direction only**, never amplitude. |
| Niec 2022, PMID 35728595 | Europe PMC record: not open access, abstract confirms *"apply BayesPrism to deconvolve expression within spatial features and develop SpaceFold to robustly map the niche at high resolution"*. | Deconvolution dependence confirmed from the abstract; the generality quotation is the agent's full-text read, not mine. |
| 8 PMIDs | Metadata and abstracts pulled directly: journal, year, OA status, PMCID. | All as reported. |

**Scope violation, declared rather than buried.** The standing instruction for this project is that
external data is out of scope. Several verdicts below lean on numbers the pass took from external
*datasets* — Cheval 2012's Table S1, the Chen 2021 GEO matrix, the Limbutara 2020 rat segment
proteome — not from publication text. Those numbers are retained here for one purpose only and are
labelled where they appear:

- **Legitimate use — novelty and exposure.** "A reviewer can produce both of our numbers from an
  open supplementary table in one minute" is a fact about priority, and it is the whole point of a
  novelty review. Keeping it is what stops us publishing a discovery claim that is already in print.
- **Forbidden use — arbitration.** No external number may confirm, calibrate or retract any
  measurement of ours. Where the pass used one that way, the verdict below is **restated as a
  question and answered on our own data instead.** That happened once and it mattered: see claim C.

## 1 · The method

| Claim | Call | Decisive prior art | Access |
|---|---|---|---|
| Recover one ordered 1D axis from unordered cross-sections of repeated structures, aggregating counts per segmented structure, using no tissue x/y | **KNOWN AND MATCHES** | Chakravorty 2026, Cell, PMID 42143021: *"The order of tubules … was determined using the de-novo cyclic pipeline in scPrisma"*; *"transcripts aggregated per tubule to generate a tubule by gene matrix"*; 337 tubules, 3 mice, no x/y, and the recovered coordinate shown to be spatially random | FULL TEXT READ (agent) |
| An **open** positional axis inside a repeated epithelial unit from unordered units with no x/y, plus the generality claim (crypts, lobules) | **KNOWN AND MATCHES** | Niec 2022, Cell Stem Cell, PMID 35728595 ("SpaceFold"): a 1D "pseudo-space" axis for intestinal crypts from Visium spots, with the generality claim in print in 2022 | FULL TEXT READ (agent); abstract verified by me |
| The residue: an open axis of position **within** a structure, from partial cross-sections at unknown positions, unknown endpoints, specimen-dependent coverage, in an epithelium whose **cell type does not change** along the axis, calibrated by balanced specimen relabelling at n = 2+2 and transferred across species | **DISCOVERY** | Negative record. No hit on target in any of the pass's or my boolean queries | — |

**Three sentences are now forbidden.** "First to order unordered cross-sections" (Chakravorty).
"Generalises to any repeated structure such as crypts and lobules" (a 2022 SpaceFold quotation).
"Aggregate per structure and order the aggregates" (the hair-follicle literature).

**What this costs and what it buys.** The novelty moves off the curve and onto the *identification
contract* — the five qualifiers above, each load-bearing. The mechanism difference is the defence
and it is clean: **no published open-axis method can run on a single-cell-type epithelium**, because
SpaceFold's axis is built from cell-type composition and needs a deconvolution reference, while
Chakravorty's prior is a closed cycle of known period that does all the identifying work and never
recovers position along a tubule's length. CoPro (PMID 42079111) cannot run without x/y and takes
the nephron ordering as an **ordinal input**.

**Two things to do before submission.** Run scPrisma's linear mode on the tubule aggregates as a
named comparator. And produce the validation the field now judges this class by: inferred position
against each cross-section's own measured diameter or area — **not** cortical depth, which this
project has already ruled out as a selector.

## 2 · The gene-level claims

| Claim | Call | Decisive source | Access |
|---|---|---|---|
| **DCXR** reverses: mouse falls S1→S3, human rises | **KNOWN AND MATCHES** | McEvoy 2022, Nat Commun, PMID 36496458, 19 living donors, Fig. 1e: *"PT2 shows increased expression of S3-abundant genes (DCXR, AGXT, SLC22A7, SLC7A13)"*. Both limbs also sit in Cheval 2012's open Table S1 (human 3.45→7.00, mouse 1.19→0.13) | FULL TEXT READ (agent) |
| **SLC17A3** reverses: mouse rises, human falls | **KNOWN AND MATCHES** | Sakaguchi 2024, PMID 38166558, Table 1: human SLC17A3 PT_S1 6.02, PT_S2 3.82, PT_S3 3.56, with *"secretion could be most active at the beginning of the PT and gradually decline"*; 3 adult male donors. Halperin Kuhns & Woodward 2021 says the same in one review sentence | FULL TEXT READ (agent) |
| **UGT3A1** reverses: mouse rises toward S3, human falls | **DISCOVERY** | The human limb is unmeasured in the only published design that could have seen it. No hit for `"UGT3A1" AND ("S3 segment" OR "pars recta")`, none for `"UGT3A1" AND "zonation"` | FULL TEXT READ (agent) |
| **ABHD14A**: mouse flat-to-falling, human rises | **DISCOVERY** | No human axial measurement exists for this gene in any source, and no renal axial claim at all | FULL TEXT READ (agent) |
| **ABHD17B**: mouse rises, human falls | **DISCOVERY** | Nothing in any organ places any ABHD17 isoform on a tissue axis | FULL TEXT READ (agent) |
| **PCK1** | **KNOWN AND MATCHES** — stays demoted | Conjard 2001 (PMID 11716765) and Drewnowski 2002 (PMID 11939717) give cross-species *agreement*, not a reversal | **ABSTRACT ONLY** — load-bearing-on-abstract |

**DCXR must be rewritten as a validation, in those words.** It is the single best evidence in the
paper that the reconstruction recovers real axial biology in two species without a shared
coordinate: a 2012 microdissection and a 19-donor single-cell study, neither carrying our probe
panels or our axis, get the same answer in both species. Written as a discovery it is the claim that
costs us the paper's credibility, because a reviewer can check it in a minute. Only the
**within-pure-S1 and within-pure-S3** half of our DCXR result is new.

**SLC17A3 likewise** — "our method recovers a published human axial gradient and adds the rodent
contrast". And say out loud for both SLC17A3 and the UGT3A family that the rodent shape is a
**mid-axis maximum**, not a monotone rise: a pure endpoint contrast cannot see that, and RNAscope
will.

**ABHD14A must stop being called a reversal.** The published rodent record has mouse flat, so the
claim is "no mouse gradient, a strong human gradient" — the wording claim A already uses for NT5E.
(Flagged as resting on an external matrix; the wording change is conservative either way, so it is
made now and the mouse-flat premise is carried as literature, not as arbitration of our fit.)

**The ABHD17B / ASCT2 mechanism is dead and should be dropped before any slide is spent on it.** The
route requires apical ASCT2 in the mouse straight segment and the published rodent data have none.
This is the pass's most valuable negative: it killed a mechanism, not a measurement.

## 3 · The framing claims

| Claim | Call | Decisive source | Access |
|---|---|---|---|
| Orthologous genes can **reverse** their axial gradient between human and mouse along a repeated epithelial unit | **KNOWN AND MATCHES** | Massalha 2020, PMID 33332768: *"SLC2A2 … exhibited pericentral zonation, in contrast to its periportal zonation in mice"*; and Yakubovsky 2026, Nature, same phenomenon class with protein validation and four outgroup mammals | FULL TEXT READ (agent); Yakubovsky verified by me |
| Species differ in **where a drug target sits** along the PT, and that is a new kind of finding | **KNOWN AND MATCHES** | Madunić 2017, PMID 28842746: mouse SGLT1 immunofluorescence with a knockout specificity control, *"the BBM of S2 segment was strongly stained"* | FULL TEXT READ (agent) |
| The human PT arm is **flatter** — median amplitude ratio 0.453, 13.1% sharper in human, smearing excluded | **DISCOVERY** | Nothing in any organ reports a species-systematic axial **amplitude** difference. I verified in full text that the paper best placed to have found one compares direction only (R = 0.15) | FULL TEXT READ, verified by me |

**Do not open on the reversal phenomenon.** It is six years old in liver and in Nature this year.
Open on the nephron and on the standard of proof — because that is where we are genuinely ahead:
Yakubovsky never separates a direction reversal from an amplitude difference, runs no within-species
sign-flip null and no amplitude surrogate, and confounds platform, lab and chemistry with species in
its cross-species arm. Our 0-of-85 and 0-of-93 within-species nulls and our amplitude surrogate are
stricter than anything in that literature. Say so, because a reviewer will otherwise assume we did
less.

**Drop "this has never been shown in the kidney."** Bonting 1960 (PMID 13802623, JCI, FULL TEXT READ
by the pass) is quantitative histochemistry along the nephron and a senior referee will remember the
series. Cite it ourselves, note that it contains no mouse, and position the pure-compartment contrast
as the fix for the compartment impurity its own footnotes record.

**Claim B keeps its novelty and must cite the literature that runs against it in the same breath:**
Albadry 2024 (PMID 38818378, FULL TEXT READ) measures cross-species gradient intensity in liver and
finds *"a high similarity between mice and humans"* with rat as the flat outlier — the opposite
direction for our species pair. And Zhai 2003 says mouse PT ultrastructure shows no segmental
variation unlike rat and human (**ABSTRACT ONLY**, so framing only, never evidence).

**One competing mechanism for claim B that we have not addressed.** Klotzer 2025 reports *"markedly
fewer PTS3 cells were detected in human kidneys"*. That is a **composition** explanation, distinct
from the positional smearing already excluded, and it predicts something testable on our own data:
condition the amplitude ratio on per-specimen axial coverage and pure-S3 unit counts and report
whether it moves. Added to open items.

## 4 · Claim C, and the one place the pass overreached

The pass called our 121-gene sharper-in-human set **CONTRADICTED** on membership, on the grounds that
UMOD and PPP1R1A are absent from microdissected rodent PT in two external datasets. The membership
worry is right and the evidence for it was the wrong kind, so the question was re-asked on our own
objects: `analysis/scripts/sharpening_set_admixture_human_vs_mouse.py`, which builds a per-bin
distal-identity score from thick-limb, distal-convoluted and collecting-duct markers in our own
units and asks which genes ride it. **Admixture adds signal where the distal tissue is, so the
hypothesis predicts a positive correlation** — the sign is the test, and it splits the set:

- The distal score is axially structured in both arms and **points the opposite way in each**: ρ with
  bin index −0.618 / −0.388 in mouse, +0.389 / +0.397 in human. It reproduces within arm (ρ 0.741
  mouse, 0.830 human). So the precondition the attack needs is present.
- **Positive riders are enriched 6.3× in the set**: 38 of 121 (31.4%) against 72 of 1,429 (5.0%) of
  the clean set.
- **83 of 121 survive** dropping the positive riders; 64 survive dropping riders of either sign.
- But the phenomenon is **not** created by admixture: across the clean set, sharpening barely tracks
  the distal score (ρ +0.174 in both human sections, +0.056/+0.079 in mouse).

**So the list was contaminated and the result was not.** Out of the headline: Umod (+0.594),
Ppp1r1a (+0.667), Ddr1 (+0.681), Abhd14a (+0.656), Pter (+0.621), Sgms2 (+0.608), Egf (+0.698).
Every gene the pass named as a membership hazard is independently flagged by our own score, which is
the best available evidence that the score works.

**Dpp4 survives, and leads.** It is an *anti*-rider (−0.579 / −0.598): highest where distal
admixture is lowest, which is what a proximal gene does. It is also the gene whose species asymmetry
is independently published in another organ with antibody validation. So the lead sentence of claim C
is a **cross-organ replication**, labelled as such — and the single highest-value experiment from
this entire round follows from it: **stain DPP4 in matched human and mouse kidney cortex.** Antibody
staining carries no probe panel, so a positive result removes the per-gene capture explanation for
the whole set in one slide, and the prediction is already pre-registered in someone else's figure.

Two cautions on the survivors: Ggact (+0.468/+0.540) and Oxct1 (+0.477/+0.704) are borderline by one
section each and should be reported with their scores attached; and Slc7a13 appears in the survivor
list but is on the project's untestable-gene list (human detection 0.12), so the untestable filter
must be applied to this set before publication.

## 5 · Still abstract-only, and therefore not publishable as evidence

Three items would each change a sentence if the full text disagrees. Get them through the library
before the species sections are written.

- **Vrhovac 2015** — the only human measurement under the S1/S2-versus-S3 SGLT2/SGLT1 scheme.
- **Zhai 2003** — mouse PT ultrastructure, the classical claim that runs against claim B.
- **Guder & Ross 1984** — enzyme distribution along the nephron.
- Also abstract-only and load-bearing where cited: **Gaganis 2007** (human renal UGT1A/UGT2B7
  described as distributed throughout the nephron), **Buckley & Klaassen 2007** (whose probe cannot
  distinguish mouse Ugt3a1 from Ugt3a2, so it supports no single-paralog claim), **Conjard 2001** and
  **Drewnowski 2002** (the PCK1 demotion — acceptable precisely because it demotes our own claim).

**One traced finding worth keeping in the discussion:** several "human facts" in the renal transport
field are rodent facts repeated. The S1/S2-versus-S3 SGLT2/SGLT1 scheme rests on one human antibody
study plus pre-1995 rat and rabbit work, and a 2026 paper finds most anti-SGLT2 antibodies
non-specific against knockout controls. That is the context in which our own SGLT1 result was ruled
untestable (human detection 0.061) — and it is an argument for the method rather than an
embarrassment, provided we say which half of the scheme anyone has actually measured in human.
