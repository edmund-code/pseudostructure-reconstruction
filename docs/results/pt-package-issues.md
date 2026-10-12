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

## Recorded deviations

13. Local `*_v4.geojson` mouse segmentations instead of the reference `v4fix` (192 more polygons); raw
    `HUK1_*_v2.geojson` for human; a nephron-wide curve with the PT as a stretch (notebook 70).
14. The notebook 74 placebo refit lowered `integration.hvg_min_disp` in a config copy, because only 45 non-zonal genes
    passed the variable-gene filter.
15. Notebook 76 loosened a reuse tolerance from 1e-6 to 1e-4 after seeing two dispersion values differ by 1.4e-5.

## Not yet done

16. The paper outline, figure plan, draft v0.6 and `pseudostructure-current-state.md` still quote SCF13-era numbers.
17. Notebook 78 (paper figures) is not built.
18. The method-comparison page shows notebook 71 only; the segment F-test and conventional DE rows of notebook 75 are
    not on it.
19. The literature check of the 83 calls (untracked scratch) names candidate new species differences: sterol
    synthesis, glyoxylate/glycine cleavage and glutathione synthesis. Its Nt5e groups (G08, G37) were classed
    inconsistently.
