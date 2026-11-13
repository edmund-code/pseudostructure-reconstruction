# Results: Figs. 5 and 6 (draft, 2026-11-12, v4: COAST fits)

Results report what the data show; clinical background and interpretation are in `discussion_fig5_fig6.md`
(user, 2026-11-12). Discovery comes from one screen applied to both comparisons: each section's curve is reduced to its
shape (zero mean, unit variance), shapes are differenced between conditions, genes are called against balanced
section swaps, and the called difference curves are clustered (`analysis/paper/shape_clusters.py`). Pathways are
tested on the same difference curves. Every number below comes from the COAST fits of the v4fix segmentations
(`analysis/paper/coast_fit.py`, 2026-11-12); v3, on the package fits of the earlier segmentations, is
`results_fig5_fig6_v3_package.md`, and the version with drug background in each section is
`results_fig5_fig6_v2_drug_stories.md`. Style: Fan-lab Results (`citation_checks/H_results_figure_style.md`).

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
pseudostructure, *Slc5a2* remained above half of its maximum through 0.26–0.29 of the mouse PT but through 0.58–0.62 of
the human PT (Fig. 5c).

### Clustering of shape differences identifies four classes of positional change between species

To find such differences without specifying their form in advance, we reduced each section's curve to its shape by
centring it and scaling it to unit variance, and computed for each gene the difference between the mean human shape
and the mean mouse shape (Methods). Of 3,868 PT-expressed genes, 1,593 changed at least 1.4-fold along the PT in both
sections of at least one species and so had a defined shape. We called a gene when its difference curve was large
(root mean square ≥ 0.5 in units of standard deviation) and at least twice as large as under either balanced swap of
one mouse and one human section between the groups. This called 666 genes, against four under each swap.

Hierarchical clustering of the 666 difference curves gave four classes (Fig. 5f). In the first (110 genes), the mean
mouse shape peaked at the end of the PT (0.90) and the human shape in the middle (0.52); it included *Slc5a2* and
*Slc7a7* (y+LAT1), whose peaks stayed early in both species (mouse 0.05–0.14, human 0.23–0.24) while the human decline
was more gradual (Fig. 5d). In the second (200 genes), expression peaked at the start of the mouse PT (0.09) and at the
end of the human PT (0.88); it included *Slc25a15* (ORNT1), which peaked at 0.16–0.30 in mouse and 0.90–0.94 in human.
The third and largest class (256 genes) was the converse, late in mouse (0.86) and early in human (0.05); it included
*Slc10a2* (ASBT), which peaked at 0.78–0.80 in mouse and at 0.05–0.32 in human. The fourth (100 genes) was mid-PT in
mouse (0.53) and late in human (0.88); it included *Slc22a12* (URAT1), which peaked at 0.53–0.54 in mouse and 0.88–0.94
in human.

The second and third classes held the genes whose gradient reversed direction between species (Fig. 5e). *Mme*
(neprilysin) rose 2.9–3.6 log2 from S1 to S3 in mouse, peaking at 0.83–0.84, but fell 1.0–1.4 log2 in human, peaking
at 0.18–0.24; *Slc17a3* (NPT4) rose 1.1–1.4 log2 in mouse and fell 1.3–1.7 log2 in human; and *Dcxr* peaked at the
start of the mouse PT (0.05–0.06) and at the end of the human PT (0.90–0.92). For each of these three genes, both
directions agree with independent microdissected mouse segments [Chen 2021] and human single-nucleus PT profiles
[Lake 2023]. The mouse direction of *Mme* also agrees with rat, in which neprilysin protein is most concentrated in S3
[Vio 2018].

### Pathway-level shape differences between species are as large as gene-level ones

We next asked whether genes of a pathway change shape together. For each of 364 Reactome, KEGG and Hallmark sets with
at least ten genes with a defined shape, we averaged the members' difference curves, so that the average is large only
when members change in the same direction at the same positions (Methods). Ninety-seven sets changed by the gene-level
criterion, and none under either section swap. Two of these also changed more than random sets of the same size drawn
from the 1,593 genes: protein secretion (13 genes) and protein processing in the endoplasmic reticulum (34; q = 0.04 for
both). Four representative sets, one from each cluster of the 97 set-level difference curves, are shown in Fig. 5g.
Oxidative phosphorylation genes (22) peaked at the start of both the mouse and human PT (0.05 and 0.21) but fell more
steeply along the human PT. Ubiquitin-specific processing proteases (16) and protein secretion (13) peaked at the end of
the mouse PT (0.86 and 0.90) and at the start of the human PT (0.05 for both), whereas glyoxylate metabolism and glycine
degradation (16) peaked mid-PT in mouse (0.51) and later in human (0.78).

Taken together, these results demonstrate that COAST, combined with clustering of shape differences, identifies where
along the PT gene expression differs between species without a predefined pattern. Applied to mouse and human, it
finds 666 genes in four classes, including early-to-late and late-to-early reversals, and pathways whose members move
together; two secretory-pathway sets move more than the genome-wide background. All human data derive from two cortical
sections of one donor.

**Fig. 5 | COAST enables comparison of positional gene expression along the mouse and human proximal tubule.** **a** One
mouse and one human kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared
pseudostructure coordinate; other tubules are grey. **b** *Slc5a2* (SGLT2) expression in the same PT tubules. **c**
*Slc5a2* expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each mouse (blue) and human
(magenta) section over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Levels
are not comparable between species because the two probe panels differ. Dashed lines mark the S1/S2 and S2/S3
boundaries of the mouse segment labels. **d** As in c for one gene from each shape class in f: *Slc7a7* (y+LAT1),
*Slc25a15* (ORNT1), *Slc10a2* (ASBT) and *Slc22a12* (URAT1). **e** Genes whose gradient reverses between species, as in
c: *Mme* (neprilysin), *Slc17a3* (NPT4) and *Dcxr*. **f** The four classes of shape difference. Each section's curve is
centred and scaled to unit variance; lines show the mean over the class's genes for each section (thin) and each species
(thick). Background shading marks S1, S2 and S3. Gene counts above each panel; two member genes are named. **g** As in
f for one pathway from each group of pathway-level shape differences.

---

### COAST resolves where acute injury reshapes the proximal tubule

We applied COAST to two control mouse kidneys and two kidneys after ischaemia–reperfusion injury (IRI), placing all
four on one shared pseudostructure coordinate (Fig. 6a). We rebuilt the axis without the 3,635 injury-responsive genes;
injured tubules kept their order as well as control tubules did (Spearman ρ 0.91 and 0.92 against 0.90 and 0.95;
Methods). Because both conditions are measured on one probe panel, levels can be compared directly here.

*Havcr1*, which encodes KIM-1, rose by 10.3, 12.6 and 9.3 log2 in S1, S2 and S3 after injury (Fig. 6b,c).

### Clustering of shape differences identifies four classes of positional change after injury

We applied the same shape-difference screen as for the species comparison. Of 5,496 PT-expressed genes, 1,788 changed
at least 1.4-fold along the PT in both sections of at least one condition. The screen called 668 genes, and none under
either balanced swap of one control and one injured section.

Clustering the 668 difference curves gave four classes (Fig. 6f). In the first and smallest (37 genes), the mean control
shape peaked mid-PT (0.49) and the injured shape earlier (0.24); it included *Khk* (ketohexokinase), which peaked at
0.48 in controls and 0.26–0.29 after injury, and *Slc26a1* (SAT1), which kept its early peak but no longer fell in S3
after injury (Fig. 6d). In the second (180 genes), expression rose toward the end of the injured PT (peak 0.88) but not
of the control PT (peak 0.34); it included *Hmox1* (heme oxygenase 1), flat in controls and peaking at 0.84–0.88 after
injury, and the DNA synthesis gene *Tk1*. In the third (206 genes), control expression rose into S3 (peak 0.83),
whereas after injury it peaked mid-PT (0.46) and fell in S3; it included *Slc6a13* (GAT2), which peaked at 0.86–0.89 in
controls and 0.42–0.44 after injury, *Slc22a2* (OCT2; 0.80–0.82 against 0.50–0.58) and *Cyp4a14*. The fourth and
largest (245 genes) followed the same direction (control peak 0.88, injured peak 0.50); it included *Acadm*, which
peaked at 0.84–0.85 in controls and 0.51 after injury, *Havcr1* and *Cyp24a1*. *Slc6a8* (CRT1), also in this class,
kept its early peak but held its level through S2 after injury before falling in S3. In the third and fourth classes,
which together hold 451 of the 668 genes, the position of peak expression moved from S3 toward the middle of the PT
after injury.

Several genes graded within a segment in controls were nearly flat there after injury (Fig. 6e). *Cyp4a14*, *Agt* and
*Cyp24a1* each changed by more than 2-fold from the start to the end of S2 in both controls (+1.1, +1.35 and −1.5 to
−1.9 log2), but by −0.5 to −0.4, +0.1 to +0.3 and +0.1 log2 after injury. To ask whether this is general, we chose
gradients in one control mouse and measured how much of each the other mice kept. The second control kept 0.78–1.33 of
each S2 gradient, whereas each injured mouse kept 0.39–0.56, and the injured mice kept less in all 12 comparisons across
the three segments. Only 4 and 7 S2 gradients belonged to genes whose overall level barely changed after injury, too
few to separate the loss due to injury from noisier placement of injured tubules (Methods).

### Pathway-level shape differences after injury exceed the genome-wide background

Of 447 sets tested, 70 changed shape by the gene-level criterion and none under either swap; 24 of these also changed
more than random sets of the same size (q ≤ 0.05). Clustering the 70 set-level difference curves gave four groups, and
the set with the strongest competitive evidence in each is shown in Fig. 6g. E2F targets (38 genes; q = 0.01) and
Myc targets (56; q = 0.01) stayed at their lowest through the first half of the injured PT and peaked in S3 (0.81
and 0.89); in controls they were also highest at the end of the PT but higher in S1 than in mid-PT. The citric acid
cycle and respiratory electron transport (43; q = 0.01) and lysosome genes (38; q = 0.01) both peaked at the start of
the control PT (0.12 and 0.13) and mid-PT after injury (0.42 for both), and fell steeply into S3 after injury.

Taken together, these results demonstrate that COAST resolves injury along the PT. After IRI, the peak of 451 genes
moves from S3 toward the middle of the PT, respiration and lysosome programmes are lost from S3, proliferation
programmes rise in S3, and gradients within segments are about half erased. With two mice per condition these
comparisons are descriptive, and they describe one early time point (about 12 h–2 d after reperfusion).

**Fig. 6 | COAST resolves where acute injury reshapes the proximal tubule.** **a** One control and one injured (IRI)
mouse kidney section assayed by Visium HD. PT tubules are coloured by their position on the shared pseudostructure
coordinate; other tubules are grey. **b** *Havcr1* (KIM-1) expression in the same PT tubules. **c** *Havcr1*
expression, in log2(counts per 10,000 + 0.5), against pseudostructure for each control (blue) and injured (orange)
mouse over pseudostructure 0.05–0.95. Points are binned means of tubules; lines are fitted curves. Dashed lines mark
the S1/S2 and S2/S3 boundaries of the control segment labels. **d** As in c for one gene from each shape class in f:
*Slc26a1* (SAT1), *Hmox1*, *Slc6a13* (GAT2) and *Slc6a8* (CRT1). **e** Three genes whose within-segment gradient is
lost after injury: *Cyp4a14*, *Agt* and *Cyp24a1*. **f** The four classes of shape difference. Each mouse's curve is
centred and scaled to unit variance; lines show the mean over the class's genes for each mouse (thin) and each
condition (thick). Background shading marks S1, S2 and S3. Gene counts above each panel; two member genes are named.
**g** As in f for one pathway from each group of pathway-level shape differences.

> Methods sentences still to write: shape screen (eligibility 1.4-fold in both sections of one condition; RMS ≥ 0.5;
> both swaps ≤ 0.5 × contrast; Ward k = 4), pathway step (≥ 10 shaped members; same rule on the mean difference
> curve; competitive p from 10,000 size-matched random sets, BH), and the 4-cluster choice (fixed, not tuned).
>
> Not recomputed on the COAST fits: the human/mouse gradient ratio of the first paragraph ("about half as steep") and
> its atlas comparison (`results/human_gradient_attenuation/`, from `analysis/scripts/attenuation_external_check.py`).
