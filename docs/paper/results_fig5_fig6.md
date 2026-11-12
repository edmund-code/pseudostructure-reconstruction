# Results: Figs. 5 and 6 (draft, 2026-11-09)

Written in the style of the Fan lab's Results sections (STalign, STdeconvolve, STcompare; read in full 2026-11-11, notes in
`citation_checks/H_results_figure_style.md`). Headings name what COAST enables, paragraphs run without bullets, genes are
italic symbols with the protein named once, and each subsection closes with a "Taken together" or "Thus" sentence. Numbers come from
`fig5-human-vs-mouse-along-pt.md`, `fig6-aki-along-pt.md` and `fig5-6-highlights.md`. Citations are given as
[Author year], with sources in `citation_checks/G_fig56_clinical_hooks.md` and `D_results_literature.md`. Post-hoc
analyses are marked as such in the text.

---

### COAST enables comparison of positional gene expression along the mouse and human proximal tubule

Much of what is known about how the proximal tubule (PT) divides its work along its length comes from mice, whereas
human kidney atlases resolve the PT into only two or three segments [Lake 2023]. As such, whether drug targets occupy
the same positions along the human PT as along the mouse PT remains largely unknown. We therefore applied COAST to two
mouse kidney sections and two cortical sections from one human kidney assayed by Visium HD, placing PT tubules from both
species on a shared pseudostructure coordinate (Fig. 5a, Methods). Human gradients were about half as steep as mouse
gradients, a difference also seen between independent human and mouse single-nucleus atlases (Methods).
Because this could reflect species, donor age or tissue source, we compared where genes are expressed rather than how
much, and evaluated each difference against swaps of one mouse and one human section between the groups (Methods).

As a first example, we focused on *Slc5a2*, which encodes the sodium–glucose cotransporter SGLT2 targeted by gliflozins.
In mouse, *Slc5a2* expression was restricted to PT tubules early along the pseudostructure coordinate, whereas in human
it extended across most PT tubules (Fig. 5b). Quantifying this along pseudostructure, *Slc5a2* remained above half of
its maximum through 0.38–0.39 of the mouse PT but through 0.66–0.71 of the human PT (Fig. 5c). In contrast, *Slc2a2*,
which encodes GLUT2, the basolateral transporter that carries the same glucose out of the cell, ended at a similar
position in both species (0.55–0.58) (Fig. 5d). The difference for *Slc5a2* was more than four times larger than under
either section swap (Methods).

Where SGLT2 ends matters for how gliflozins act. Empagliflozin and dapagliflozin are approved for type 2 diabetes, heart
failure and chronic kidney disease [Tran 2024]. SGLT2 reabsorbs about 97% of filtered glucose, and when it is
inhibited, SGLT1 further along the PT reabsorbs part of the glucose that escapes, so that only 50–60% of filtered glucose
is excreted in mice [Rieg 2014]. Modelling in the rat further suggests that SGLT2 inhibition shifts sodium reabsorption
downstream and roughly doubles oxygen consumption in S3 [Layton 2015]. Both effects depend on how much of the PT lies
beyond SGLT2. Because SGLT2 extends further along the human PT, we speculate that the stretch available for SGLT1
compensation, and the segment that takes on the extra workload, may be shorter in human than in mouse. As such, mouse
estimates of SGLT1 compensation and of late-PT workload under gliflozins may need to be adjusted for the position of
SGLT2 in human, which may also bear on dual SGLT1 and SGLT2 inhibition with sotagliflozin, approved for heart failure [Tran 2024].

We next applied the same analysis to transporters that carry drugs and urate across the PT. *Slc22a8*, which encodes
the organic anion transporter OAT3, peaked at the start of the human PT (0.06–0.07) but further along the mouse PT
(0.24–0.28). *Slc22a12*, which encodes the urate transporter URAT1, peaked late in the human PT (0.83) but mid-PT in
mouse (0.54) (Fig. 5d). Both differences exceeded either section swap at least four-fold. *Slc22a6*, which encodes
OAT1, varied 48- to 289-fold along the mouse PT, peaking at 0.47, but only 1.6- to 2.3-fold along the human PT. In
contrast, *Slc22a2* (OCT2) and *Abcc2* (MRP2) varied less than 3.5-fold along the PT of either species, too little to
assign them a position (Methods).

These positions bear on where drugs enter PT cells. OAT1 and OAT3 sit on the basolateral membrane of PT cells, where
they take up many anionic drugs, among them antibiotics, non-steroidal anti-inflammatory drugs, diuretics and
antivirals, from the blood for secretion into urine; regulatory agencies highlight both as key transporters for drug
excretion and drug–drug interactions [Nigam 2015]. Uptake on the basolateral side is loosely coupled to efflux into
urine by the apical pumps MRP2 and MRP4 [Nigam 2015]. The antiviral tenofovir is one such substrate: in mice given
tenofovir, PT mitochondria took on an irregular shape with sparse cristae, and mice lacking OAT1 were spared, placing
OAT1 uptake upstream of the change [Kohler 2011]. In mouse, OAT1 was concentrated around the middle of the PT, whereas
in human it was spread nearly evenly and OAT3 was concentrated at the start. We speculate that the PT cells that take
up OAT substrates may therefore sit at different positions in mouse and in human. Similarly, URAT1 sits on the apical
membrane of PT cells and reabsorbs filtered urate; it is the target of the uricosuric drugs benzbromarone and
lesinurad, which are approved to lower serum urate in gout [Guo 2025]. Its later position in human suggests that these
drugs act on a later stretch of the PT in human than in mouse. Thus, COAST can test whether the position of a drug
target along the PT transfers from mouse to human.

### COAST identifies genes whose positional gradient reverses between species

A shallower gradient in human could reflect donor or tissue rather than species. In contrast, a gradient that rises
along the PT in one species and falls in the other cannot arise from attenuation alone. We therefore screened 3,872
PT-expressed genes for those whose S3-versus-S1 difference had one sign in both mouse sections and the opposite sign in
both human sections, with at least a 1.4-fold difference in each (Methods). We found 28 such genes, and none under
either section swap. To evaluate whether these reversals reproduce in independent data, we compared each species with
microdissected mouse segments from both sexes [Chen 2021; Chen 2023] and with human single-nucleus profiles of PT
segments [Lake 2023]. For 11 of the 28 genes, both directions reproduced (Fig. 5e). Of the remainder, four had a mouse
gradient that weakens or reverses in female mice [Chen 2023; Xiong 2023], and the others either did not reproduce in one species or were
detected by few probes.

For example, *Mme*, which encodes the peptidase neprilysin, rose 2.7–3.5 log2 from S1 to S3 in
mouse, peaking at 0.82–0.83, but fell 0.8–1.0 log2 in human, peaking at 0.10–0.11 (Fig. 5e). Conversely, *Dcxr*, which
encodes the uronate-pathway enzyme dicarbonyl/L-xylulose reductase, peaked at the start of the mouse PT (0.06) and at the
end of the human PT (0.90–0.95), and *Slc17a3*, which encodes the urate exporter NPT4, followed *Mme*. The mouse direction
agrees with rat, in which neprilysin protein is most concentrated in the S3 segment [Vio 2018]; the human direction
has not been described.

Neprilysin is the main enzyme that degrades natriuretic peptides, and it also processes bradykinin, substance P,
angiotensin II and endothelin; of all tissues, it is most abundant on the brush border of PT cells
[Judge 2015]. Its inhibitor sacubitril, given with valsartan, lowered cardiovascular death and heart failure
hospitalization in patients with reduced ejection fraction [McMurray 2014; Jhund 2015]. In mouse, neprilysin was
concentrated at the end of the PT, whereas in human it was concentrated at the start. We speculate that the PT cells in
which sacubitril acts may therefore differ between species, so that renal effects of neprilysin inhibition measured in
mice may not arise in the same cells in patients. More generally, a gene whose gradient reverses cannot be placed on
the human PT by transferring its mouse segment, and each such gene calls for positional measurement in human tissue.

### COAST locates where paired transporters hand over along the proximal tubule

Many PT transporters act in pairs, in which one family member is expressed early along the PT and is replaced by a
related member later on. We reasoned that the position at which the two cross, which we refer to as the switch point,
could differ between species even when both genes are expressed in both. We located the switch point on each gene's
curve scaled to its own range, so that it does not depend on how efficiently each probe detects its gene. Because this
scaling stretches a nearly flat curve to the full range, we counted a switch point only where both partners changed at
least two-fold along the section (Methods; post hoc). Of 12 transporter pairs that met this criterion in all four
sections, 7 switched later in human and none earlier, whereas none did when sections were swapped between species
(Methods). Most of the other pairs failed the criterion in human, where gradients are shallower. The seven included
*Slc5a2* → *Slc5a10* (SGLT2 → SGLT5; +0.21 of the PT length), *Slc6a19* → *Slc6a18* (B0AT1 → B0AT3; +0.18) and
*Slc7a7* → *Slc7a13* (y+LAT1 → AGT1; +0.24) (Fig. 5f,g). In each, the early partner was maintained
further along the human PT; human *SLC6A19* remained the more abundant transcript through S3 (Fig. 5f), as it did in an
independent human single-nucleus atlas [Lake 2023]. All six later-in-human pairs that could be evaluated in that atlas
reproduced. These included *Slc5a2* → *Slc5a1* (SGLT2 → SGLT1; +0.30), whose *Slc5a1* partner is detected in too few of
our human tubules to be shown here, and whose later switch agrees with the longer SGLT2 domain above.

A later switch point means that the early partner serves more of the PT in human, which matters when the early partner
is a drug target. B0AT1 reabsorbs neutral amino acids, including phenylalanine, and is being pursued as a target in
phenylketonuria, in which phenylalanine accumulates in the blood. Mice lacking *Slc6a19* excrete over 100-fold more
phenylalanine in urine, and loss of *Slc6a19* improves a mouse model of phenylketonuria [Belanger 2018]. In mouse,
B0AT3 takes over from B0AT1 in S3, whereas in human B0AT1 remains dominant through S3. We speculate that inhibiting
B0AT1 may therefore block a larger share of neutral amino acid reabsorption in human than in mouse, with less of the PT
left to B0AT3. As such, the urinary loss achieved by a B0AT1 inhibitor in patients may differ from that predicted from
mice.

### Segment-level clustering misses half of the switch-point shifts

We next sought to evaluate whether these shifts could be recovered by clustering analysis alone. To this end, for each
section we fit straight lines through the S1, S2 and S3 cluster means, as a segment-level analysis would, and located
where each pair crossed (Methods). If the shifts were apparent at segment resolution, we would expect the cluster means
to give shifts similar to those from COAST. The cluster means recovered most of the shift for the three amino acid
pairs (Fig. 5g; post hoc). However, they recovered only +0.02 of the +0.21 shift for *Slc5a2* → *Slc5a10* and +0.07
of the +0.25 shift for *Slc22a23* → *Slc22a12*, and gave no switch point for *Slc5a2* → *Slc5a11* (SGLT2 → SMIT2).
Every human switch point lies within S2 (0.53–0.60), where three cluster means place a single value. We note that
straight lines through three cluster means are one simple segment-level baseline, and finer clustering may recover
more of each shift. Thus, a map of drug targets built from segment-level clusters would miss half of these
differences, including the later handover from SGLT2 to its downstream partners.

Taken together, these results demonstrate that COAST can compare positional gene expression across species on a shared
continuous coordinate. Applied to the PT, it reveals that the targets of gliflozins, uricosurics and B0AT1 inhibitors
occupy different positions in mouse and human, that neprilysin, the target of sacubitril, reverses its gradient, and
that transporter switch points move later in human, half of which segment-level clustering misses. We note that all
human data derive from two cortical sections of one donor and that transcript position need not reflect protein
position. Antibody staining for SGLT2, OAT3, URAT1, B0AT1 and neprilysin along the human PT in additional donors would
be needed to confirm these differences.

**Fig. 5 | COAST enables comparison of positional gene expression along the mouse and human proximal tubule.** **a** One
mouse and one human kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared
pseudostructure coordinate; other tubules are grey. **b** *Slc5a2* (SGLT2) expression in the same PT tubules. **c**
*Slc5a2* expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each mouse (blue) and human (magenta)
section over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Levels are not
comparable between species because the two probe panels differ; positions are. Dashed lines mark the S1/S2 and S2/S3
boundaries of the mouse segment labels. **d** As in c for *Slc2a2* (GLUT2), *Slc22a6* (OAT1), *Slc22a8* (OAT3) and
*Slc22a12* (URAT1). **e** Genes whose gradient reverses between species, as in c: *Mme* (neprilysin), *Slc17a3* (NPT4)
and *Dcxr*. **f** Three transporter pairs whose switch point moves later in human, in log2(1 + counts per 10,000). For
each species, the early (light) and late (dark) partner are shown as the mean of the two sections' curves (thick) and
each section's curve (thin), over the range used to locate the switch; each species has its own y axis, so heights
compare the two partners within a species, not mouse with human. Background shading marks S1, S2 and S3. Switch points
are located on scaled curves (Methods). **g** Switch point in mouse (start of arrow) and human (arrowhead)
for the six shifted pairs whose partners are both in the screened gene set, from COAST curves (black) and from straight
lines through the S1, S2 and S3 cluster means (grey, dashed; post hoc). Background shading marks S1, S2 and S3.


> Fig. 5 section restructured 2026-11-12 by capability (setup; positions of drug targets; reversed gradients; handovers; what clustering misses), each with one drug story. Figure changes this text assumes: d drops OCT2 (changes < 3.5-fold in both species); new e = reversals (*Mme*, *Slc17a3*, *Dcxr*); old e, f become f, g. Drug sources read in full 2026-11-12 (citation_ledger.md); McMurray 2014 still unread (paywalled).

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
