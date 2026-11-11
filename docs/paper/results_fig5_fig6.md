# Results: Figs. 5 and 6 (draft, 2026-11-09)

Written in the style of the Fan lab's Results sections (STalign, STdeconvolve, STcompare; read in full 2026-11-11, notes in
`citation_checks/H_results_figure_style.md`). Headings name what COAST enables, paragraphs run without bullets, genes are
italic symbols with the protein named once, and each subsection closes with a "Taken together" or "Thus" sentence. Numbers come from
`fig5-human-vs-mouse-along-pt.md`, `fig6-aki-along-pt.md` and `fig5-6-highlights.md`. Citations are given as
[Author year], with sources in `citation_checks/G_fig56_clinical_hooks.md` and `D_results_literature.md`. Post-hoc
analyses are marked as such in the text.

---

### COAST enables comparison of positional gene expression along the mouse and human proximal tubule

The proximal tubule (PT) is the main site of solute and drug uptake in the kidney, and much of what is known about how
it divides this work along its length comes from mice. However, human kidney atlases resolve the PT into only two or
three discrete segments [Lake 2023], and cross-species comparisons have quantified how much of each transporter is
present rather than where along the PT it is expressed [human cortex proteomics, PMC5118637; cross-species proteomics,
PMID 31123036]. As such, whether transporters occupy the same positions along the mouse and human PT remains an open
question. We therefore applied COAST to two mouse kidney sections and two cortical sections from one human kidney
assayed by Visium HD, placing PT cross-sections from both species onto a shared pseudostructure coordinate (Fig. 5a,
Methods). Human PT tubules spanned the full coordinate, although about half lay in the range that mouse segment labels
assign to S1, compared to a third of mouse tubules (Extended Data Fig. X). Because human gradients were consistently shallower than mouse gradients, for reasons we could not separate in
this cohort (Methods), we min-max scaled each gene's fitted expression within each section in order to compare where
genes are expressed rather than how much.

As a first example, we focused on *Slc5a2*, which encodes the sodium–glucose cotransporter SGLT2 targeted by
gliflozins. In mouse, *Slc5a2* expression was restricted to PT tubules early along the pseudostructure coordinate,
whereas in human it extended across most PT tubules (Fig. 5b). Quantifying this along pseudostructure, *Slc5a2*
remained above half of its maximum through 0.38–0.39 of the mouse PT but through 0.66–0.71 of the human PT (Fig. 5c).
In contrast, *Slc2a2*, which encodes GLUT2, the basolateral transporter that carries the same glucose out of the cell,
ended at a similar position in both species (0.55–0.58) (Fig. 5d). To evaluate whether such differences could arise
from variability between sections, we swapped one mouse and one human section between the groups. The difference for
*Slc5a2* was more than four times larger than under either swap (Methods).

We next applied the same analysis to transporters that mediate renal drug secretion. *Slc22a6* and *Slc22a8*, which
encode the organic anion transporters OAT1 and OAT3, peaked at the start of the human PT (0.06–0.16) but in the
early-to-mid PT in mouse (0.24–0.47) (Fig. 5d). Likewise, *Abcc2*, which encodes the efflux pump MRP2, was expressed
from the start of the human PT but only past the middle of the mouse PT (Extended Data Fig. X). In contrast, *Slc22a2*,
which encodes the organic cation transporter OCT2, rose later along the human PT (0.51–0.57 versus 0.32–0.36), and
*Slc22a12*, which encodes the urate transporter URAT1, peaked late in the human PT (0.84) but mid-PT in mouse (0.54)
(Fig. 5d). Each of these differences again exceeded both section swaps at least four-fold. Because OAT1, OAT3 and OCT2
mediate the basolateral uptake of many renally cleared drugs [OCT/MATE review, PMC3675737; to verify], we speculate
that these positional differences may change which PT cells first encounter a given drug in human compared to mouse.

Many PT transporters act in pairs, in which one family member is expressed early along the PT and is replaced by a
related member later on. We reasoned that the position at which the two cross, which we refer to as the switch point,
could differ between species even when both genes are expressed in both. Of 53 transporter pairs with a switch point
in all four sections, 13 (25%) switched at a different position in human than in mouse, whereas none did when sections
were swapped between species (Methods). Eight of these 13 pairs switched later in human, including *Slc5a2* →
*Slc5a10* (SGLT2 → SGLT5; +0.21 of the PT length), *Slc22a23* → *Slc22a12* (+0.25), *Slc22a23* → *Slc22a2* (+0.24),
*Slc7a7* → *Slc7a13* (y+LAT1 → AGT1; +0.24) and *Slc6a19* → *Slc6a18* (B0AT1 → B0AT3; +0.18) (Fig. 5e). In each, the early partner was maintained further along the human PT. All seven later-in-human pairs that
could be evaluated in an independent human single-cell atlas [Lake 2023] reproduced. These included *Slc5a2* →
*Slc5a1* (SGLT2 → SGLT1), whose *Slc5a1* partner is detected in too few of our human tubules to be shown here.

We next sought to evaluate whether these shifts could be recovered by clustering analysis alone. To this end, for each
section we fit straight lines through the S1, S2 and S3 cluster means, as a segment-level analysis would, and located
where each pair crossed (Methods). If the shifts were apparent at segment resolution, we would expect the cluster means
to give shifts similar to those from COAST. However, the cluster means recovered only +0.07 of the +0.25 shift for
*Slc22a23* → *Slc22a12*, +0.06 of the +0.24 shift for *Slc22a23* → *Slc22a2* and +0.02 of the +0.21 shift for
*Slc5a2* → *Slc5a10* (Fig. 5f; post hoc). They gave no switch point for three pairs and moved one monocarboxylate
transporter pair in the opposite direction. Several of these switch points lie within S2, which three cluster means
cannot resolve. We note that straight lines through three cluster means are one simple segment-level baseline, and
finer clustering may recover more of each shift. Thus, these species differences are visible only along a continuous coordinate.

Finally, cross-species comparisons typically restrict analysis to genes with one-to-one orthologs. However, several PT
gene families, such as *Cyp4a*, *Cyp2d*, *Ugt2b* and *Slco1a*, expanded in mouse. We therefore asked whether genes
without a one-to-one human counterpart carry positional information of their own. In mouse, 9.9% of genes from
families with more mouse than human copies and 14.2% of genes from families expanded in both species were strongly
zonated along the PT, compared to 3.8% of one-to-one orthologs (Mantel–Haenszel odds ratio 2.24, P = 8 × 10⁻⁷)
(Fig. 5g). In contrast, mouse genes whose human counterpart had expanded showed no such excess (1.7%). As such, part of
the positional organization of the mouse PT resides in gene copies without a single human counterpart and would be
omitted by a one-to-one comparison.

Taken together, these results demonstrate that COAST can compare positional gene expression across species on a shared
continuous coordinate. Applied to the PT, it reveals that transporters for glucose, drugs and urate occupy different
positions in mouse and human, and that their switch points shift in a manner largely missed by segment-level
clustering. We note that all human data derive from two cortical sections of one donor and that transcript position
need not reflect protein position. Antibody staining for SGLT2, OAT1, OCT2 and URAT1 along the human PT in additional
donors would be needed to confirm these differences.

**Fig. 5 | COAST enables comparison of positional gene expression along the mouse and human proximal tubule.**
**a** One mouse and one human kidney section assayed by Visium HD. PT tubules are coloured by their position on the
shared pseudostructure coordinate; other tubules are grey. **b** *Slc5a2* (SGLT2) expression in the same PT tubules.
**c** *Slc5a2* expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each mouse (blue) and human
(magenta) section over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Levels
are not comparable between species because the two probe panels differ; positions are. Dashed lines mark the S1/S2 and
S2/S3 boundaries of the mouse segment labels.
**d** As in c for *Slc2a2* (GLUT2), *Slc22a6* (OAT1), *Slc22a8* (OAT3), *Slc22a2* (OCT2) and *Slc22a12* (URAT1).
**e** Five transporter switches that move later in human. For each species, the early (light) and late (dark) partner
are shown as the mean of the two sections' min-max scaled curves (thick, shaded) and each section's curve (thin), over
the range used to locate the switch. Background shading marks S1, S2 and S3. Dashed blue and dotted magenta lines mark
the mouse and human switch points. Right, each section's switch point with its bootstrap 95% interval (mouse blue,
human magenta). **f** Shift in switch point (human −
mouse) for the 12 pairs that switched at a different position in human and whose partners are both measured, from
COAST curves (filled) and from straight lines through the S1, S2 and S3 cluster means (open; × where the cluster means
give no switch). **g** Percentage of PT-expressed genes strongly zonated along pseudostructure in mouse and human, by
orthology class (human : mouse gene copies). Error bars are Wilson 95% confidence intervals.

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
