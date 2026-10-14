# Open issues after the switch to the package pseudospace

Recorded on 2026-10-12, after notebooks 70–77 reran every coordinate- and label-dependent result on Caleb's package
(`pseudospace_reconstruction` 0.1.0, stages 01–04). Each item names its source. None is resolved; the paper outline,
figure plan and draft still quote the old numbers, and notebook 78 (figures) is not built.

## Paper framing

1. **"Continuous beats discrete" is weakly supported.** The segment F-test (the positional F-test with S1/S2/S3 steps
   in place of the spline) passes the swap rule with 76 human vs mouse calls, 71 of them shared with the positional
   F-test. The spline adds 12 redundancy clusters, 10 of them near misses; the steps add 5 (notebooks 74, 75). The
   supported message is narrower: testing shape against specimen shape noise finds differences that level-based
   clustering and conventional DE miss or confuse with probe-panel differences (notebooks 71, 75).
2. **The package coordinate is noisy within segments.** Not met (notebook 72): human gene-fold agreement 0.49, equal
   depth 0.73 (human S1 0.55), species registration gap 0.37, 10 of 12 cells noise-dominated, held-out-gene gains.
   One fold refit reverses the mouse PT. Every within-segment result is therefore descriptive (notebook 74).
3. **Our landmark coordinates beat the package within segments.** LCP and LRS-Z: gene-fold 0.93–0.97 and 0 of 12
   noise-dominated cells on the same units (notebook 72). Their place in the paper (comparator, dropped, or a
   within-segment add-on) is undecided.
4. **Atlas replication is not set-specific.** The 83 calls replicate in both atlases against random sets and noise
   calls, but not against random genes of the same positional strength (p 0.35, 0.64; old 0.0013, 0.0057;
   notebook 77). Replication supports genes, not gene sets.
5. **The human vs mouse F-test passes calibration narrowly:** 0.071 against the 0.075 limit (notebook 71).

## Labels and units

6. **Human PT-S2 shrinks to 632 units** (373 + 259; old 1,813). The package labels 47% and 60% of the old human S2
   units as S1 (notebook 70). S2−S1 classified genes fall from 9,013 to 6,269, and part of the S2−S1 agreement is
   undefined (notebook 73). HUK1_MED1's segment transitions cross on the curve (notebook 72, C1).
7. **The package drops much of the injured AKI S3.** An Unassigned cluster of 2,425 units (72% IR) is dropped; 569 of
   its units sit on the PT stretch, 556 in S3, 501 from IR mice. An 80-unit injured-PT cluster (Havcr1, Spp1, Krt20,
   Lcn2) is labelled cTAL. Deep-S3 units: IR 251 and 32 against control 783 and 700 (notebooks 70, 76). The labels
   were confirmed without overrides; this is a limitation of every AKI result.

## AKI

8. **Only two AKI claims hold** (notebook 76): the injury-marker positive control and a specific cell-cycle program.
   Transport loss, its progression, late-PT depletion and the decoupling of the two programs are not claimable. The
   outline's "progressive, segment-dependent decoupling" sentence is unsupported.
9. **The AKI screens fail calibration:** 0.120 (notebook 71) and 0.144 (notebook 76's universe). No AKI call list can
   be presented as a discovery list.
10. **The IRI time point is unknown.** The cell-cycle program reads as "undetermined" between division and arrest.

## The package itself (worth telling Caleb)

11. **Stage 03 is chaotic.** The BLAS/OpenMP thread count, or a lognorm layer with unsorted indices, changes its
    clusters and moves PT units (notebook 72). Runs must keep the default threads and file layout.
12. From the scratch trial (2026-10-11): `plain_log2fc_level` subtracts the typical-gene curve, unlike its
    documentation; stage 01 forks one worker per sample (about 22 GB RAM for four samples); the 150-gene cap on gene
    sets drops the mTORC1 and OXPHOS Hallmarks; the QuPath feature `id` is ignored and `Nphs1` is skipped silently.
13. **The pass-2 PCA uses 631 genes, not the 763 it records.** Pass-1 clustering writes `var['highly_variable']`
    (its 2,000 genes) onto the unit object; the pass-2 object inherits it, and scanpy 1.11.5's `sc.tl.pca` then uses
    that column as its gene mask. Pass 2 therefore runs PCA on the 631 genes that are both integration genes and
    pass-1 clustering genes, while `uns['pca']` says 763. Checked on the human vs mouse run (2026-10-12): a refit on
    the 631 genes reproduces the stored PCs exactly; a refit on all 763 correlates 0.989–0.998 over PCs 1–5. The
    effect on the embedding is small, but `X_harmony`, the curve and DPT all inherit it, and it should hold for the
    AKI run and the reference run too (not checked). Not patched.

## Recorded deviations

14. Local `*_v4.geojson` mouse segmentations instead of the reference `v4fix` (192 more polygons); raw
    `HUK1_*_v2.geojson` for human; a nephron-wide curve with the PT as a stretch (notebook 70).
15. The notebook 74 placebo refit lowered `integration.hvg_min_disp` in a config copy, because only 45 non-zonal genes
    passed the variable-gene filter.
16. Notebook 76 loosened a reuse tolerance from 1e-6 to 1e-4 after seeing two dispersion values differ by 1.4e-5.

## Not yet done

17. The paper outline, figure plan, draft v0.6 and `pseudostructure-current-state.md` still quote SCF13-era numbers.
18. Notebook 78 (paper figures) is not built.
19. The method-comparison page shows notebook 71 only; the segment F-test and conventional DE rows of notebook 75 are
    not on it.
20. The literature check of the 83 calls (untracked scratch) names candidate new species differences: sterol
    synthesis, glyoxylate/glycine cleavage and glutathione synthesis. Its Nt5e groups (G08, G37) were classed
    inconsistently.

## From notebooks 81 and 82 (published findings re-read along the PT)

20. **"No published table compares the S1→S3 profile across species" needs qualifying.** Cheval et al. 2012 (PMID
    23056504) microdissected S1 and S3 in both human and mouse (SAGE). They report global divergence but no per-gene or
    per-program profile comparison. Our species × position interaction replicates in their data (notebook 82, B2b).
21. **Human zonation is the weakest measurement.** Against microdissected human S1/S3, our cortex-section human S3 − S1
    agrees at ρ 0.14 (mouse 0.46).
22. **Within-segment positions.**
    - Raw within-species program curves carry a shape that random gene sets share; notebook 81 centres by the median
      gene.
    - After that, only two mouse gradients are confirmed, and no human position. Human mid-S1 turns may mark the old
      S2 units the package labels as S1 (issue 6).
23. **Two `test_program_localization.py` tests** fail instead of skipping when R cannot find limma (an unactivated
    shell).
24. **Mouse S1 units carry neighbour signal that human S1 units barely carry.** In both mice, S3 − S1 (centred) is
    −3.3 to −5.0 for DCT genes (S100g, Slc12a3, Calb1, Pvalb, Slc8a1), −3.7 to −4.3 for podocyte Nphs2 and −2.0 to
    −2.3 for Ren1. Our human sections show at most −1.2 for the same genes. TAL genes (Umod, Slc12a1) rise toward S3
    in both species. This fits S1 at the vascular pole next to the DCT and glomerulus, and the 10 µm S1 odds ratio
    (notebook 42).
    - **Consequence:** S100g's mouse S1 signal (−5.0, −4.7) is absent in Census mouse PT (−0.07; notebook 82,
      section 10b). Raghubar's S100g relocation is therefore not reproduced as PT biology.
    - **Bound:** a gene's mouse-only S1 bias can be neighbour signal. B2b replicates in microdissected tissue, so
      neighbour signal does not drive the program-level interaction. Single-gene positional claims should be checked
      against an atlas, as section 10b does.
25. **Within-segment pathway positions (notebook 84) are not supported on the elastic curve.**
    - Chen's segment-level pathways do not separate inside S1 or S2 beyond random gene groups (p 0.12, 0.13).
    - PTS3 enriches too few pathways to test.
    - Separation appears only when non-PT-specific genes enter, which is the neighbour-signal route of issue 24.
    - Testing S3 needs another discrete source of PTS3 pathways, for example a single-nucleus atlas's S3 cluster or
      the clustering method's S3 calls. That would be a new protocol.
    - Peak-sorted heatmaps show a diagonal by construction, so they need the reportable marks.

## From notebook 85 (sharp or gradual transitions)

26. **Pooled segment units misplace mouse positions.**
    - Notebooks 79–84 read positions in units pooled over all four specimens (cuts 0.473, 0.605). The mouse's own
      S1→S2 label transition is at 0.341, which is pooled S1 0.72.
    - So 83's "switches inside S1" (NaPi-IIc S1 0.69, U-M3 S1 0.72) are at the mouse's own S1/S2 boundary (mouse units
      S1 0.96, S1 1.00).
    - 84's pooled "S2 thirds" were the second half of mouse S2.
    - Within-species statements should use the species' own units. Human units are undefined, because the human
      transitions cross (issue 6).
27. **The last stretch of the mouse curve is thin-limb units labelled PT-S3.**
    - From mouse S3 ≈ 0.3, the package's `primary_identity` turns Thin_Limb: 33–62% at S3 0.4 and 52–95% afterwards.
    - In total, 425 and 547 mouse "PT" units carry it.
    - 18 sharp-looking falls of PT pathways (S3 0.38–0.66) sit at this switch, and so does 83's NaPi-IIa switch (S3
      0.36; 29% thin-limb units).
    - Late-S3 positions are not within-PT statements until these units are separated. Like issue 7, this is a
      labelling limit of the package.
28. **71's 6-df spline cannot show a sharp switch near S2.**
    - A perfect step reads 1.04 pooled units wide at the mouse S2/S3 boundary and 0.81 in S2, with up to 24% overshoot.
    - 81's "no sharp switch" and every width read from 71's curves (81, 83, 84) carry this floor.
    - Notebook 85 measures widths on 40 equal-count bins instead (floor 0.05).

## From notebook 86 (genes outside one-to-one orthology)

29. **The Cyp4a family disagrees with male microdissection.** Our curve and our labels put Cyp4a12a in early S2 (1.19
    mouse units) and Cyp4a10, Cyp4a14 and Cyp4a31 at the S2/S3 boundary (2.08–2.20). Chen 2023's male microdissection has
    Cyp4a10 and Cyp4a14 highest in PTS2 and Cyp4a12a in PTS3. 62 found the same order on the retired coordinate, so the
    coordinate is not the cause. Cyp4a12a's probe status is "near", and its microdissection level is low (9/26/55 TPM).
    Any Cyp4a positional statement needs an orthogonal check (ISH).
30. **Genes with no included probe are structural zeros.** In our panels, mouse Gsta1, Nat8f6, Nat8f7 and Apol9b and
    human UGT2B15 and AKR1C1 are absent. Gsta1 is the mouse partner of the canonical human PT marker GSTA1. A "human-only
    alpha-GST" result would be a panel artefact.
31. **Family sums fail the 1:1 test's calibration** (relabeled share 0.089 against the 0.075 limit; notebook 86, H7), and
    the 468 AKI genes outside 76's pool fail gene-level calibration (0.114). Neither yields claimable gene or family
    calls; only the class-level comparisons are reported.

## From notebook 89 (positional programs against a cluster-level view)

32. **The S1-high pathway gradient (P1) may be partly neighbour signal.** Its sets are generic: ECM organisation, EMT,
    focal adhesion, interferon-γ and TNF-α responses, hypoxia, mitotic spindle. Issue 24 found mouse S1 units carrying
    neighbour signal. P1 falls at S2/S3 and inside S3 rather than at S1/S2, which argues against a pure S1 artefact, but
    no positional claim for P1's immune and ECM members should be made without a neighbour-signal check.
33. **E2's cluster-view threshold was mis-specified.** The protocol called a segment when its centred value reached
    0.25 log2, the same number as the curve's amplitude (a range). A centred segment value is about half a range, so
    the cluster view "missed" 120 of 123 sets largely by construction. Any future cluster-versus-curve comparison should
    use a statistical cluster test (the segment F-test) or a range-matched threshold.

