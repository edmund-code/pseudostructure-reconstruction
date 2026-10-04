# Open items for the coordinator (workstream 4: Introduction and Methods)

Ordered by how much they block the manuscript. Each item says what is missing, where I looked, and
who can resolve it.

## A. Missing study facts (needs lab records; cannot be recovered from the repository)

1. **Mouse specimens.** Strain (the brief says C57BL/6-type; nothing in the repository records
   it), age, supplier, housing, euthanasia, and animal-protocol approval. Sex is inferred only from
   Y-linked transcripts (notebook 31 statement); confirm from records.
2. **Human donor.** Tissue source (nephrectomy margin, deceased donor, or other), age, sex from
   records, kidney function, consent/IRB. State whether `HUK1_COR1` and `HUK1_MED1` are serial
   sections or come from separate blocks, and why one is named "MED". The data show it is cortex,
   and the drafts say so.
3. **Tissue processing and imaging.** Fixation (FFPE or fixed frozen), section thickness, H&E
   protocol, scanner, and Visium HD workflow. A CytAssist image exists for the human sections;
   the mouse workflow is unknown. Also needed: library preparation, sequencer and read depth.
4. **Space Ranger differs between species.** Mouse used 3.1.1 (mm10, Mouse Probe Set v2.0) and human
   used 4.0.1 (GRCh38-2024-A, Human Probe Set v2.1.0). The Methods report this; decide whether
   reprocessing the mouse data with 4.0.1 is needed or whether reporting suffices. (I read the
   versions from h5 attributes and web summaries, which are metadata, not raw data.)
5. **Data availability.** Deposition plan for Visium HD data and polygons, code archive DOI and
   licence (`docs/publication-readiness.md` lists these as open too).

## B. Segmentation provenance (unverifiable from code)

6. **Human segmentations.** `HUK1_*_v2.geojson` are described as "re-segmentations". Nothing records
   how they were made (kidney_panoptic or not, which checkpoint, whether D4 test-time augmentation
   was used, manual edits), or whether they were validated. The model was trained only on mouse
   H&E. A reviewer will ask.
7. **Upstream mouse QC.** The rules and code that produced `*_kept_tubules_labeled_fine.geojson`
   (`kept_tubule`, `doublet_tier1_flag`, removal of low-support structures) are not in the
   repository. About 7.6k of 13.6k v4 instances were kept for `Ctrl1A2`. The Methods need this
   procedure. Also confirm that the v4 + D4 test-time-augmentation inference in
   `segmentation/model_summary_v2.txt` §7c is the run that produced these files.
8. **Coordinate frame.** Notebook 01 assigns bins with Space Ranger full-resolution pixel
   coordinates (scale factor 1), so the GeoJSON must be in that frame. Confirm how the H&E used for
   segmentation was registered to the Visium HD image. This matters most for human, where 2,331 of
   8,747 (`COR1`) and 6,963 of 14,186 (`MED1`) polygons receive no bins. Explain why: capture area
   smaller than the segmented area, or misregistration.
9. **Segmentation performance numbers.** Decide whether to report them: held-out PQ 0.880 with D4
   test-time augmentation on 98 objects from one nephrectomy slide; patch-level PQ 0.727 on 1,298
   objects. `Ctrl1A2` was partly training data, and human segmentation is unvalidated.
10. **Segmentation software versions** used in the actual inference run (PyTorch 2.5.1 and
    Python 3.10 are pinned; OpenCV and OpenSlide versions in the run are unknown). The OpenMidnight
    encoder has no peer-reviewed citation, only a model card and blog post; choose the citation.

## C. Inconsistencies found in the pipeline that the paper must resolve or disclose

11. **Notebook 13 and notebook 03 use different partitions.**
    - Both have 26,839 structures and 11 clusters, but different membership hashes (13:
      `41fd77bb98da`; 03: `4ebd19739331`) and different cluster numbering.
    - Sizes differ slightly: AL 859 vs 858; Unassigned 1,890 vs 1,897; pass-2 nephron 24,335 vs
      24,328; PT 12,866 vs 12,871.
    - Notebook 13 says its labels were "confirmed against the reviewed notebook 03 memberships", but
      no concordance is saved.
    - Recommendation: describe only notebook 13's partition in the paper and report its agreement
      with 03 (adjusted Rand index and a cross-tabulation) in a supplement.
12. **The sampled-region check uses notebook 03's glomerulus labels.** Notebook 33 reads glomeruli
    from `human_vs_healthy_mouse/cross_species_harmony_pass1.h5ad` (03), while PT labels come from
    13. Either switch to 13's glomerulus cluster or state the provenance.
13. **The DPT trajectory-method comparison is not reproducible from committed code.** It uses an
    earlier run of notebook 12 (logic `12-pathway-remodeling-v5-bh-by-statistic`, 11,106 genes,
    1,516 pathways, 126 calls) on 03's DPT coordinate (`results/pt_pathway_remodeling/`). Current
    notebook 12 reads only 13. Fix by adding a coordinate switch or a committed notebook, and rerun
    under whatever coordinate workstream 1 picks.
14. **The Census extraction script is not committed** (`SCRATCH/census_pseudobulk.py`). It
    hard-codes machine-specific absolute paths and reads `results/` for the gene list. It must become
    a committed script or notebook with root flags. Do the same for the GEO and Xiong downloads:
    record URLs, dates and checksums in `docs/`, not only in the gitignored
    `data/external/*/SOURCE.txt`. Record the Xiong repository commit; the local files were renamed
    to `PT-S{1,2,3}_fWTvsmWT.csv`.
15. **The specificity rule was amended after the first run.** The second condition (relabeled
    calls < 5% of tested pathways) was added after a pathway-score model passed the NCDR alone. The
    Methods disclose this. Workstream 2 should decide whether to keep both conditions.
16. **Two kinds of q-value.** Specificity counts use normal-approximation q-values, while discovery
    uses empirical matched q-values. The relabeling step of the funnel also uses the normal q. The
    Methods state this; workstream 2 should confirm or unify.
17. **Docs say "13-cluster" but there are 11.** AGENTS.md and `docs/workflows/pseudospace.md`
    refer to a "13-cluster reference fingerprint". The cross-species partitions have 11 clusters;
    the 13 probably refers to the mouse-only notebook 02. Clarify in the docs.
18. **Uncommitted machine path.** The working copy of notebook 01 contains a commented machine
    path (a `%env PSEUDOSPACE_VISIUM_ROOT=...` line with a local absolute path). HEAD is clean. Do not stage that cell.

## D. Scientific points for workstream 1 (coordinate), which the Methods depend on

19. **Orientation markers overlap the biology.** The marker axis uses Gatm, Slc5a2 and Slc5a12
    (early) and Slc22a7, Slc7a13 and Cyp7b1 (late). Gatm is a headline biology gene (creatine
    zonation shift) and Cyp7b1 is a male-biased mouse S3 gene. The axis only chooses the root tip,
    but say so explicitly, or show that the coordinate is unchanged without these genes.
20. **Projection uncertainty is not quantified.** scFates pseudotime uses one projection
    (`n_map = 1`). Consider more projections to quantify uncertainty.
21. **The coordinate check reuses the labels.** The same reviewed S1/S2/S3 labels both check the
    coordinate and anchor the registration warp. External datasets (notebook 34) are the
    independent check; keep that framing.
22. **Common support truncates the human tail.** Support is [0.018, 0.864]; human 99th percentiles
    are about 0.87, mouse reaches 1.0. Late mouse S3 beyond 0.864 is excluded from the gene
    models. Every count in the Methods (12,272 structures, 11,071 genes, 1,513 pathways, segment
    tables, transition positions) must be refreshed if the coordinate changes.

## E. Dependencies on workstreams 2 and 3

23. **Workstream 2: final reporting set.** Which list is primary: confident (82), robust (66) or
    core (42)? Confirm the funnel order and the definition of "specimen-sensitive". The
    *Sensitivity runs and the final pathway list* section is a placeholder.
24. **Workstream 2: two-relabeling ceiling.** Only two balanced relabelings exist. Methods and
    Discussion should keep "screens for gross nonspecificity; not a calibrated FDR".
25. **Workstream 3: literature-check Methods.** Methods should follow WS3's final procedure.
    Confirm the novelty sentence in the Introduction ("no study has compared human and mouse PT
    gene programs position by position").
26. **Workstream 3: CYP2E1 conflict.** Literature conflicts on human renal CYP2E1 (PMID 38042273
    localises it to human PT). I left CYP2E1 out of the Introduction.

## F. Introduction claims still marked [CITE]

27. Continuous (rather than stepwise) PT zonation within and across S1–S3. I could not find an
    abstract-level source; WS3 may know a segment-resolved study that shows this.
28. A classic or review reference for S3 ischaemic vulnerability. Currently Kim et al. 2011 (mouse,
    PMID 22025970). Witzgall et al. 1994 (PMID 7910173) is an alternative, but its abstract does
    not name the species.
29. A review stating that the PT is the main site of drug-induced kidney injury. Currently a 2026
    in-vitro paper (PMID 41699311) whose abstract says so.

## G. Small items to fill

30. **HCOP table.** Download date and version; local sha256 prefix `0cfb78e4eb273751`.
31. **Pathway libraries.** Source and date for the Enrichr-format files (`Reactome_2022`,
    `MSigDB_Hallmark_2020`, `KEGG_2019_Mouse`). `GO_Biological_Process_2023` is present but unused.
32. **Census.** Confirm that `2025-11-08` is the release used and whether it is an LTS release. The
    human dataset's collection DOI is a preprint (10.1101/2025.03.06.637075); the published
    version appears to be Genome Biol 2026 (PMID 41821037), with the same first author and the same
    dataset description. Confirm.
33. **Rat symbol matching.** Rat genes were matched to mouse symbols case-insensitively, which is
    approximate. Consider an ortholog table, or state the limitation (the Methods state it).
