# Results: Figs. 5 and 6 (draft, 2026-11-12, v3)

Results report what the data show; clinical background and interpretation are in `discussion_fig5_fig6.md`
(user, 2026-11-12). Discovery comes from one screen applied to both comparisons: each section's curve is reduced to its
shape (zero mean, unit variance), shapes are differenced between conditions, genes are called against balanced
section swaps, and the called difference curves are clustered (`analysis/paper/shape_clusters.py`). Pathways are
tested on the same difference curves. The previous version, with drug background in each section and the transporter
handover panels, is `results_fig5_fig6_v2_drug_stories.md`. Style: Fan-lab Results (`citation_checks/H_results_figure_style.md`).

---

### COAST enables comparison of positional gene expression along the mouse and human proximal tubule

Human kidney atlases resolve the proximal tubule (PT) into two or three segments [Lake 2023], so whether a gene occupies
the same stretch of the human PT as of the mouse PT is largely unknown. We applied COAST to two mouse kidney sections
and two cortical sections from one human kidney assayed by Visium HD, placing PT tubules from both species on a shared
pseudostructure coordinate (Fig. 5a). Human gradients were about half as steep as mouse gradients, a ratio also seen
between independent human and mouse single-nucleus atlases (Methods). Because the two species are measured on different
probe panels, we compared the shape of each gene's profile along the PT rather than its level or amplitude.

*Slc5a2*, which encodes SGLT2, illustrates the comparison. In mouse its expression was restricted to PT tubules early
along the pseudostructure coordinate, whereas in human it extended across most PT tubules (Fig. 5b). Along
pseudostructure, *Slc5a2* remained above half of its maximum through 0.38–0.39 of the mouse PT but through 0.66–0.71 of
the human PT (Fig. 5c).

### Clustering of shape differences identifies four classes of positional change between species

To find such differences without specifying their form in advance, we reduced each section's curve to its shape by
centring it and scaling it to unit variance, and computed for each gene the difference between the mean human shape
and the mean mouse shape (Methods). Of 3,872 PT-expressed genes, 892 changed at least 1.4-fold along the PT in both
sections of at least one species and so had a defined shape. We called a gene when its difference curve was large
(root mean square ≥ 0.5 in units of standard deviation) and at least twice as large as under either balanced swap of
one mouse and one human section between the groups. This called 434 genes, and none under either swap.

Hierarchical clustering of the 434 difference curves gave four classes (Fig. 5f). In the first (63 genes), expression
peaked at the start of the mouse PT (0.05) but in the middle of the human PT (0.52); it included *Slc5a2* and
*Slc7a7* (y+LAT1), whose peaks stayed early in both species while the human decline was more gradual (Fig. 5d). In the
second (86 genes), mouse expression rose to a mid-PT peak (0.56) whereas human expression peaked early (0.14); it
included *Slc22a6* (OAT1), which peaked at 0.47 in mouse and 0.15 in human. The third and largest class (158 genes)
rose toward the end of the mouse PT (peak 0.86) and fell along the human PT (peak 0.12); it included *Slc10a2*
(ASBT), which peaked at 0.78 in mouse and 0.15 in human. The fourth (127 genes) was the converse, mid-PT in mouse
(0.46) and late in human (0.88); it included *Slc22a12* (URAT1), which peaked at 0.54 in mouse and 0.83 in human.

The third and fourth classes held the genes whose gradient reversed direction between species (Fig. 5e). *Mme*
(neprilysin) rose 2.7–3.5 log2 from S1 to S3 in mouse, peaking at 0.82–0.83, but fell 0.8–1.0 log2 in human, peaking
at 0.10–0.11; *Slc17a3* (NPT4) followed the same pattern; and *Dcxr* peaked at the start of the mouse PT (0.06) and at
the end of the human PT (0.90–0.95). For each of these three genes, both directions agree with independent
microdissected mouse segments [Chen 2021] and human single-nucleus PT profiles [Lake 2023]. The mouse direction of
*Mme* also agrees with rat, in which neprilysin protein is most concentrated in S3 [Vio 2018].

### Pathway-level shape differences between species are as large as gene-level ones

We next asked whether genes of a pathway change shape together. For each of 122 Reactome, KEGG and Hallmark sets with
at least ten genes with a defined shape, we averaged the members' difference curves, so that the average is large only
when members change in the same direction at the same positions (Methods). Forty-one sets changed by the gene-level
criterion, and none under either section swap. None, however, changed more than random sets of the same size drawn from
the 892 genes (lowest q = 0.12). Four representative sets, one from each cluster of the 41 set-level difference curves,
are shown in Fig. 5g. Lysosome genes (23) peaked at the end of the mouse PT (0.89) and at the start of the human PT
(0.08). Bile acid metabolism (27 genes) and peroxisomal protein import (21) peaked mid-PT in mouse (0.61 and 0.57) and
late in human (0.85 for both), whereas fatty acid metabolism (53) peaked mid-PT in mouse (0.58) and early in human
(0.12).

Taken together, these results demonstrate that COAST, combined with clustering of shape differences, identifies where
along the PT gene expression differs between species without a predefined pattern. Applied to mouse and human, it
finds 434 genes in four classes, including early-to-late and late-to-early reversals, and pathways whose members move
together, although no pathway moves more than the genome-wide background. All human data derive from two cortical
sections of one donor.

**Fig. 5 | COAST enables comparison of positional gene expression along the mouse and human proximal tubule.** **a** One
mouse and one human kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared
pseudostructure coordinate; other tubules are grey. **b** *Slc5a2* (SGLT2) expression in the same PT tubules. **c**
*Slc5a2* expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each mouse (blue) and human
(magenta) section over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Levels
are not comparable between species because the two probe panels differ. Dashed lines mark the S1/S2 and S2/S3
boundaries of the mouse segment labels. **d** As in c for one gene from each shape class in f: *Slc7a7* (y+LAT1),
*Slc22a6* (OAT1), *Slc10a2* (ASBT) and *Slc22a12* (URAT1). **e** Genes whose gradient reverses between species, as in
c: *Mme* (neprilysin), *Slc17a3* (NPT4) and *Dcxr*. **f** The four classes of shape difference. Each section's curve is
centred and scaled to unit variance; lines show the mean over the class's genes for each section (thin) and each species
(thick). Background shading marks S1, S2 and S3. Gene counts above each panel; two member genes are named. **g** As in
f for one pathway from each group of pathway-level shape differences.

---

### COAST resolves where acute injury reshapes the proximal tubule

We applied COAST to two control mouse kidneys and two kidneys after ischaemia–reperfusion injury (IRI), placing all
four on one shared pseudostructure coordinate (Fig. 6a). We rebuilt the axis without the 3,812 injury-responsive genes;
injured tubules kept their order nearly as well as control tubules did (Spearman ρ 0.90 and 0.91 against 0.94 and 0.97;
Methods). Because both conditions are measured on one probe panel, levels can be compared directly here.

*Havcr1*, which encodes KIM-1, rose by 10.6, 12.3 and 9.3 log2 in S1, S2 and S3 after injury (Fig. 6b,c). Within S2 it
climbed along the tubule after injury (+1.0 log2 from the start to the end of S2, against −0.2 in controls) (Fig. 6c).

### Clustering of shape differences identifies four classes of positional change after injury

We applied the same shape-difference screen as for the species comparison. Of 5,521 PT-expressed genes, 1,744 changed
at least 1.4-fold along the PT in both sections of at least one condition. The screen called 685 genes, against one
and none under the two balanced swaps of one control and one injured section.

Clustering the 685 difference curves gave four classes (Fig. 6f). In the first (193 genes), expression rose toward the
end of the injured PT (peak 0.89) but not of the control PT (peak 0.39); it included *Hmox1* (heme oxygenase 1), flat
in controls and rising into S3 after injury, and the DNA synthesis gene *Tk1* (Fig. 6d). In the second (236 genes),
control expression rose into S3 (peak 0.84), whereas after injury it peaked mid-PT (0.47) and fell in S3; it included
*Slc6a13* (GAT2), which peaked at 0.87 in controls and 0.43 after injury. The third (194 genes) followed the same
direction with a deeper control trough in S2 (control peak 0.84, injured peak 0.49); it included *Slc22a2* (OCT2),
which peaked at 0.82 in controls and 0.52 after injury, *Acadm* and *Cyp4a14*. The fourth and smallest (62 genes)
peaked early in controls (0.10) and mid-PT after injury (0.53), and included *Havcr1* and *Cyp24a1*. *Slc6a8* (CRT1),
also in this class, kept its early peak but held its level through S2 after injury before falling in S3. In the
second and third classes, which together hold 430 of the 685 genes, the position of peak expression moved from S3
toward the middle of the PT after injury.

Several genes graded within a segment in controls were nearly flat there after injury (Fig. 6e). *Cyp4a14*, *Agt* and
*Cyp24a1* each changed by more than 2-fold from the start to the end of S2 in controls (+1.3, +1.6 and −2.1 log2),
but by −0.4, +0.4 and +0.1 log2 after injury. To ask whether this is general, we chose gradients in one control mouse
and measured how much of each the other mice kept. The second control kept 0.84–1.16 of each S2 gradient, whereas each
injured mouse kept 0.49–0.56, and the injured mice kept less in all 12 comparisons across the three segments. Genes
whose overall level barely changed after injury kept 0.66–0.68 of their S2 gradient, against 0.45–0.54 for
injury-responsive genes, so about a quarter of the loss is attributable to noisier placement of injured tubules
(Methods).

### Pathway-level shape differences after injury exceed the genome-wide background

Of 428 sets tested, 66 changed shape by the gene-level criterion and none under either swap; 26 of these also changed
more than random sets of the same size (q ≤ 0.05). Clustering the 66 set-level difference curves gave four groups, and
the set with the strongest competitive evidence in each is shown in Fig. 6g. E2F targets (39 genes; q = 0.007) and
Myc targets (55; q = 0.007) stayed at their lowest through the first half of the injured PT and peaked in S3 (0.84
and 0.92); in controls they were also highest at the end of the PT but higher in S1 than in mid-PT. Mitochondrial fatty acid β-oxidation (20; q = 0.013) peaked at 0.61 in
controls and 0.49 after injury and fell steeply into S3 after injury. Lysosome genes (38; q = 0.007) peaked at the
start of the control PT (0.13) and mid-PT after injury (0.44), as did the citric acid cycle and respiratory electron
transport (36; q = 0.013; peak 0.14 to 0.47).

Taken together, these results demonstrate that COAST resolves injury along the PT. After IRI, the peak of 430 genes
moves from S3 toward the middle of the PT, fatty acid oxidation, respiration and lysosome programmes are lost from S3,
proliferation programmes rise in S3, and gradients within segments are about half erased. With two mice per
condition these comparisons are descriptive, and they describe one early time point (about 12 h–2 d after
reperfusion).

**Fig. 6 | COAST resolves where acute injury reshapes the proximal tubule.** **a** One control and one injured (IRI)
mouse kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared pseudostructure
coordinate; other tubules are grey. **b** *Havcr1* (KIM-1) expression in the same PT tubules. **c** *Havcr1*
expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each control (blue) and injured (orange)
mouse over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Dashed lines mark
the S1/S2 and S2/S3 boundaries of the control segment labels. **d** As in c for one gene from each shape class in f:
*Hmox1*, *Slc6a13* (GAT2), *Slc22a2* (OCT2) and *Slc6a8* (CRT1). **e** Three genes whose within-segment gradient is
lost after injury: *Cyp4a14*, *Agt* and *Cyp24a1*. **f** The four classes of shape difference. Each mouse's curve is
centred and scaled to unit variance; lines show the mean over the class's genes for each mouse (thin) and each
condition (thick). Background shading marks S1, S2 and S3. Gene counts above each panel; two member genes are named.
**g** As in f for one pathway from each group of pathway-level shape differences.

> Methods sentences still to write: shape screen (eligibility 1.4-fold in both sections of one condition; RMS ≥ 0.5;
> both swaps ≤ 0.5 × contrast; Ward k = 4), pathway step (≥ 10 shaped members; same rule on the mean difference
> curve; competitive p from 10,000 size-matched random sets, BH), and the 4-cluster choice (fixed, not tuned).
