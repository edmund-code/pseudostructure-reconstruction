# Results: Figs. 5 and 6 (draft, 2026-11-09)

Written in the style of the Fan lab's Results sections (STalign, STdeconvolve). Headings are claims, each paragraph
reports one result, and each subsection closes with a "Thus" sentence. Numbers come from
`fig5-human-vs-mouse-along-pt.md`, `fig6-aki-along-pt.md` and `fig5-6-highlights.md`. Citations are given as
[Author year], with sources in `citation_checks/G_fig56_clinical_hooks.md` and `D_results_literature.md`. Post-hoc
analyses are marked as such in the text.

---

### COAST reveals that drug transporters occupy different stretches of the mouse and human proximal tubule

Much of what is known about how the proximal tubule (PT) divides its work along its length comes from mice, and this
knowledge underpins renal pharmacology. Human kidney atlases resolve the PT into two or three segments [Lake 2023], and
comparisons between species have measured how much of each transporter is present rather than where it sits [human
cortex proteomics, PMC5118637; cross-species proteomics, PMID 31123036]. We therefore asked three questions that need a
continuous position along the tubule. Do transporters occupy the same stretch of the PT in both species? Do paired
transporters hand over at the same place? And does a gene-by-gene comparison miss part of the picture? We placed PT
cross-sections from two mouse kidneys and two sections of one human kidney on a shared pseudostructure axis and
followed each gene from the start to the end of the PT (Methods). Human gradients were consistently shallower than
mouse gradients, for reasons we could not separate in this cohort (Methods), so we compared where genes are
expressed, not how much.

**Transporters sit at different places.** We scaled each gene's curve in each section from its lowest to its highest
value and asked where it is high (Fig. 5a).
- **Glucose.** SGLT2, the target of gliflozins, stayed above half its maximum through 0.66–0.71 of the human PT,
  against 0.38–0.39 in mouse. GLUT2, which carries the same glucose out of the cell, ended at the same place in both
  species (0.55–0.58).
- **Organic anions.** OAT1 and OAT3 peaked at the start of the human PT (0.06–0.16) but in early-to-mid PT in mouse
  (0.24–0.47). The efflux pump MRP2 started at the beginning of the human PT but past its middle in mouse.
- **Organic cations and urate.** OCT2 started later in human (0.51–0.57 against 0.32–0.36). The urate transporter
  URAT1 peaked late in human (0.84) and mid-PT in mouse (0.54).

Each of these differences was at least four times larger than when we swapped one mouse and one human section between
the groups. In mouse, SGLT2, GLUT2, OAT1 and OAT3 were lost from the end of the PT more steeply than the overall
flattening of human gradients could explain. Many drug targets nonetheless kept the direction of their gradient in both
species, including SGLT2, OAT3, URAT1, ACE and PSMA, while a few reversed, including neprilysin (MME, the target of
sacubitril) and 15-PGDH (HPGD) (Fig. 5b).

These differences bear on how mouse results transfer to human. A longer SGLT2 zone would put more of the human PT under
gliflozin action. OAT1, OAT3 and OCT2 bring many drugs into PT cells from the blood [OCT/MATE review, PMC3675737; to
verify], so where they sit determines which cells meet a drug first.

**Switch points move.** Many transporters work in pairs: one member is high early in the PT and is replaced by a
related member later on. We call the point where the two cross the switch point. Of 53 pairs with a switch in all
four sections, 13 switched at a different position in human, and none did when sections were swapped between species.
Eight switched later in human. They include:
- the sugar pair SGLT2 → SGLT5 (+0.21 of the PT length);
- the urate pair SLC22A23 → URAT1 (+0.25);
- SLC22A23 → OCT2 (+0.24);
- the amino-acid pairs y+LAT1 → AGT1 (+0.24) and B0AT1 → B0AT3 (+0.18) (Fig. 5c).

In each, the early human partner stayed on further along the tubule. All seven later-in-human pairs we could check in a
human single-cell atlas [Lake 2023] reproduced. One of them is SGLT2 → SGLT1; its SGLT1 partner is detected too rarely
in our human sections to measure.

**Detail that clustering misses.** We next asked whether these shifts could be read from the S1, S2 and S3 cluster means
alone. Straight lines through the three cluster means placed most human switches close to where the mouse switched.
They recovered only +0.07 of the +0.25 URAT1 shift, +0.06 of the +0.24 OCT2 shift and +0.02 of the +0.21 SGLT5 shift,
and one monocarboxylate pair moved the wrong way (Fig. 5d; post hoc). Several of these switches move within S2, which
three cluster means cannot resolve, so these species differences are visible only on a continuous axis.

**Gene families the one-to-one comparison cannot see.** Species comparisons usually keep only genes with one copy in each
species. Some PT genes instead belong to families that expanded in mouse, such as CYP4A, CYP2D, UGT2B and SLCO1A. We
asked whether these genes carry PT zonation of their own.
- In mouse, 9.9% of genes from families with more mouse than human copies were strongly zonated, and 14.2% from
  families expanded in both species. Only 3.8% of one-to-one genes were (odds ratio 2.2, P < 10⁻⁵) (Fig. 5e).
- Mouse genes whose human counterpart had expanded showed no excess (1.7%).

Part of the mouse PT's positional organisation therefore lives in gene copies that have no single human counterpart,
and a one-to-one comparison silently drops it.

Thus, COAST shows that the human PT is not a scaled copy of the mouse PT. Transporters for glucose, drugs and urate
occupy different stretches, and their switch points sit later, at a resolution that segment-level clustering largely
misses. Mouse findings about where a transporter works in the PT should therefore be checked by position before they
are carried to human. All human results come from one donor and two cortical sections, and transcript position need
not match protein position. Antibody staining for SGLT2, OAT1, OCT2 and URAT1 along the human PT in more donors would
test these conclusions directly.

---

### COAST resolves where acute injury reshapes the proximal tubule

Acute kidney injury damages the PT unevenly. The S3 segment, the last part of the PT, is usually hit hardest
[Ransick 2019]. We compared two control mouse kidneys with two kidneys after ischemia–reperfusion injury (IRI) on one
shared axis. To check that injury did not distort the axis, we rebuilt it without 3,812 injury-responsive genes.
Injured tubules kept their order nearly as well as control tubules (Methods).

**Healthy gradients flatten inside S2.** Within S2, expression of the injury biomarker KIM-1 (HAVCR1) climbed along
the tubule after injury (+1.0 log2 from the start to the end of S2, against −0.3 in controls) (Fig. 6a). At the same time, healthy S2
gradients of the vitamin D enzyme CYP24A1, angiotensinogen (AGT) and CYP4A14 flattened. These single-gene changes
did not meet our false-call threshold (Methods), so we show them as examples. To ask whether flattening is general,
we chose gradients in one control mouse and measured how much of each the other mice kept. The second control kept
most of each gradient (0.84–1.16 of it in S2), while each injured mouse kept about half (0.49–0.56). Injured mice
kept less in all 12 comparisons across the three segments (Fig. 6b).

**Proliferation rises toward S3.** Cell-cycle genes rose after injury in every segment, and most in S3: CDK4 rose
+0.9 log2 in S1 and +1.9 in S3 (Fig. 6d). In all, the injury response of 268 genes differed by at least 1 log2 between S3 and S1, against 0 and 4 genes
when control and injured sections were swapped. Within S2, the DNA replication programme climbed
along the tubule (Fig. 6c). Among 122 pathways whose injury response changed within a segment (none in either swap, after a
correction described in Methods), DNA replication, G1/S transcription and mismatch repair rose through S2 (Fig. 6e). CDK4/6 inhibitors given
before injury reduce early AKI in mice [DiRocco 2014; Pabla 2015], although longer treatment promoted fibrosis
[Zhou 2024].

**Mitochondrial loss concentrates in S3.** Oxidative phosphorylation fell along S3 after injury (Fig. 6f), and the
respiratory-chain pathways were among the S3 pathway calls (Fig. 6h). Other drug targets also responded most in S3
(Fig. 6g):
- xanthine oxidase (XDH, the target of allopurinol and febuxostat) rose +1.3 log2 in S1 and +2.6 in S3;
- 11β-HSD1 fell −3.1 in S3;
- the cystine transporter rBAT fell −1.8 in S3;
- HNF1B, a transcription factor that maintains PT identity, fell −1.5 in S3.

**Switch points move upstream.** Of 52 transporter pairs, 11 switched at a different position after injury, and
none did when sections were swapped (post hoc). Most moved toward the start of the PT (Fig. 6i). They included:
- OAT1 → OCT2, the main uptake transporters for many renally cleared drugs (−0.14 of the PT length);
- 4F2hc → rBAT, where rBAT is mutated in cystinuria (−0.20);
- B0AT1 → GAT2 (−0.23);
- the folate pair SLC46A3 → PCFT (−0.22).

The zinc pair ZIP5 → ZIP8 moved the other way (+0.17); loss of ZIP8 protects mice from ischemic injury [Cai 2025].
The cystine pair b0,+AT → AGT1 kept its position, although AGT1 fell about 20-fold. This shows that a switch point
can stay put while expression levels change. Here, unlike between species, straight lines through the cluster means
recovered the direction of every shift. The sizes differed: clustering gave −0.07 for OAT1 → OCT2 (Fig. 6j).

Thus, COAST shows that injury does more than raise or lower genes in each segment. It flattens healthy gradients,
shifts transporter switch points toward the start of the PT, and concentrates proliferation and mitochondrial loss
in S3. With two mice per group, these comparisons are descriptive and need confirmation in larger cohorts.
