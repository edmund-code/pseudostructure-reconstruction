# Open items for the paper (after draft v0.3)

This file replaces `docs/paper/open_items.md`. Items are ordered by who must act and by how much they
block submission. Each says what is missing, where it appears in `draft_v03.md`, and what would close
it. Resolved items from the earlier list are summarised at the end.

## A. The user must supply (lab records; cannot be recovered from the repository)

1. **Human pixel scale (item 13).**
   - Notebook 42 reads the pixel size from each GeoJSON feature's `Area px` and `Area µm²` fields:
     0.441 µm/px for mouse, but exactly 0.500 µm/px for both human files (`HUK1_*_v2.geojson`).
   - The segmentation model works at 0.44068 µm/px. Human bin coverage (0.765) matches the 0.757
     expected if the human µm fields are mis-scaled.
   - If the true human scale is 0.44 µm/px, human distances and lengths are overstated by 13.5%
     (corrected glomerulus distances 131 / 148 / 228 µm) and areas by 29%.
   - Unaffected: the bin-to-polygon join (pixel coordinates) and the structure filters (counts).
   - Affected: notebook 42's human distances (R1, Figure S10) and notebook 11's extent analysis.
   - **Needed:** the scanner metadata (objective, µm per pixel) of the human H&E images, and whether
     the human `v2` segmentation ran at 0.44068 µm/px. Draft: R1, Methods 6.3 and 6.25, Figure S10.
2. **Mouse specimens.** Strain (the plan says C57BL/6-type; nothing records it), age, supplier,
   housing, euthanasia. Sex from records: all four specimens express Y-linked transcripts, which is the
   only current evidence (Methods 6.1).
3. **Human donor.**
   - Tissue source (nephrectomy margin, deceased donor or other), age, kidney function, comorbidity
     and medication, and pre-operative fasting.
   - Ischaemia times; sex from records; consent.
   - Whether `HUK1_COR1` and `HUK1_MED1` are serial sections or separate blocks, and why one is named
     "MED".
   - This matters for the Limitations "Shared confounders" paragraph and the M4 response (Methods 6.1;
     Discussion).
4. **Ethics.** Animal-protocol and human-tissue (IRB) approvals, with numbers (Methods 6.1).
5. **Tissue processing, imaging and sequencing.** Fixation and embedding, section thickness, H&E
   protocol, scanner, the Visium HD workflow for each species (a CytAssist image exists for human),
   library preparation, sequencer and depth (Methods 6.1–6.2).
6. **Segmentation provenance** (Methods 6.3–6.4; Table S8).
   - How the human `*_v2.geojson` segmentations were made: model and checkpoint, D4 test-time
     augmentation, manual edits, pixel scale (item 1), and whether they were validated. The model
     was trained on mouse H&E only.
   - Rules and code of the upstream mouse QC that produced `*_kept_tubules_labeled_fine.geojson`.
   - How the H&E images were registered to the Visium HD frame. Notebook 42 shows exact centroid
     agreement and that zero-bin polygons lie outside the capture area, but not the method.
   - Whether to report segmentation PQ: held-out PQ 0.880 with D4 test-time augmentation on 98 objects
     from one nephrectomy slide. `Ctrl1A2` was partly training data; human segmentation is
     unvalidated.
   - Software versions of the inference run (PyTorch 2.5.1 and Python 3.10 are pinned; OpenCV and
     OpenSlide are unknown).
7. **Space Ranger versions differ** (mouse 3.1.1 with mm10 and probe set v2.0; human 4.0.1 with
   GRCh38-2024-A and v2.1.0). Decide whether reporting suffices or the mouse data should be reprocessed
   (Methods 6.2; Limitations).
8. **Data and code availability.** Deposition of the Visium HD data and polygons (controlled access
   for human tissue), code archive DOI and licence (Methods 6.31).

## B. Reference inputs and wording (user or coordinator)

9. **HCOP table.** Download date and version (local sha256 prefix `0cfb78e4eb273751`; Methods 6.5).
10. **Pathway libraries.** Source and date of the Enrichr-format files (`Reactome_2022`,
    `MSigDB_Hallmark_2020`, `KEGG_2019_Mouse`; Methods 6.13). The round 3 legends list a
    WikiPathways (W) tag, but Table S2's 1,513 pathways contain none (Reactome 1,181, KEGG 282,
    Hallmark 50). The tag was dropped from the v0.3 legend; the figures agent can drop it from
    notebook 41's label helper.
11. **Harmony θ = 6** (referee B29; Methods 6.7).
    - No current notebook records why θ = 6 was chosen.
    - A stale diagnostic, `human_vs_healthy_mouse/diagnostics/pass1_clustering_sensitivity.csv`
      (2026-09-08), is no longer written by notebook 03. Its reference partition matches θ = 2
      (ARI 1.0), and θ = 6 scores 0.42–0.46, so it cannot support the choice.
    - Needed: a rationale, or a current sensitivity table.
12. **Literature check wording** (referee B34; Methods 6.22). State plainly that AI agents searched
    Europe PMC with defined queries, and archive the queries.
13. **OpenMidnight citation** (Kaplan et al. 2025). Only a model card and blog post exist; choose the
    citation form.
14. **GSE267280 and GSE277302** (notebook 44 fasting and time-of-day data). Find the associated
    publications, if any (Methods 6.27 [CITE]).
15. **UGT1A9 probes** (referee B25; R5). Show that both panels' UGT1A9 probes target the unique exon 1.
    Otherwise the calls may measure the UGT1A locus.
16. **S3 vulnerability source** (referee B35). Kim et al. 2011 is a label-retaining-cell study; a
    dedicated source would be better.
17. **Census release.** Confirm that 2025-11-08 is the release used and whether it is LTS. Confirm that
    the human dataset's preprint (10.1101/2025.03.06.637075) is the Genome Biol 2026 paper (PMID
    41821037).

## C. Coordinator checks and consistency across documents

18. **Transporter count.** 33 of 43 zonated transporters conserved under symmetric classes, against
    196 of 388 zonated genes. Counted by workstream 4 from `gene_classes_symmetric.csv`, not by a
    notebook. Confirm, or add the count to notebook 43. Item 9's 30 of 54 is the notebook 40 count;
    R3 gives both.
19. **Abedini PMID.** `data/external/abedini2024_visium/SOURCE.txt` cites PMID 38514613, which is a
    different paper. The correct PMID is 39048792 (Nat Genet 56:1712–1724). The file is local
    (gitignored), and any committed doc that copied it needs the same fix.
20. **Documents that still carry superseded statements.**
    - `docs/paper/r2_revision_results.md`: 72% (table 0.715 → 71%).
    - `docs/paper/revision_state_citations.md`: rat OAT1 S3/S2 −6.6 (table −6.55 → −6.5); Gamt "first
      protein-level support".
    - `docs/paper/r3_revision_wording.md`: Pkhd1 human-only (indeterminate under the symmetric union).
    - `docs/paper/r5_draft.md` and `docs/results/pt-pathway-literature-novelty.md`: old reversal labels
      and the OAT "mouse-only" module. Under symmetric classes Slc13a3 is conserved and Slc22a6
      indeterminate.
    - `docs/paper/r1_draft.md`: 19 / 47 (now 10 / 26).
21. **Notebook 45 anchor.** Notebook 45's re-implementation of the step screen gives 50 species calls
    in our cohort, against 16 in notebook 37; relabelings are 0 and 0 in both. If Figure 2g shows our
    cohort, say which implementation the star uses.
22. **Analyses requested by the referee and not run.** Decide whether to run each or to state it as
    a limitation:
    - species-separate clustering with marker review, and R3 without integration genes (M1 fixes 2
      and 3). Cross-fitted transfer was run instead;
    - a program-level or permutation comparison to replace pathway-level Mann–Whitney p-values (M7
      fix 4);
    - class-level enrichments under the symmetric classes. The draft cites notebook 40's only, as
      superseded;
    - medullary-ray location of human S3 (M1 fix 1).
23. **Table S2 coordinate columns.** Merge the `coordinate_robustness.csv` columns ("DPT13-robust"
    and the others) into `table_S2_pathways.csv`, or ship them as a second sheet. Assemble Table 1
    from notebook 43, and Tables S9–S11 from notebooks 45, 42/46 and 44.
24. **Orientation markers overlap the biology** (earlier item 19). Gatm, Slc5a2 and Slc5a12 (early)
    and Slc22a7, Slc7a13 and Cyp7b1 (late) set the root tip. Methods 6.9 says that the axis only
    chooses the root. Keep.
25. **Unrecorded downloads.** Record URLs, dates and checksums for the GEO, KTEA and Xiong downloads
    in `docs/`, not only in gitignored `SOURCE.txt` files.

## D. Figures (round 3 and later)

26. **Figure 5.** The coordinator's gene set is drawn. Recommendations:
    - show **Gclm** instead of Gclc. Gclm is a notebook 38 "headline" gene with balanced probes; Gclc
      is "supporting" with unequal probes (3 against 6). Both are indeterminate;
    - consider adding **Rbp4**, the strongest human-specific signal (expressed only in human; rise in
      7 of 8 Visium samples). Its card would show the human side only;
    - Cyp2e1 shows only the rodent side (not expressed in human). Say so in its card.
27. **Figure 4a.** Replace or supplement the historical "DPT (notebooks 03, 12)" row with a DPT13 row
    from `coordinate_robustness.csv` (item 15).
28. **Figure S11** (notebook 44 state figure). Assemble in notebook 41, or drop the references to it in
    R3.
29. **Figure 2f.** Put our cohort's relabelings on the same rows as their strategies, with black stars
    (item 15).
30. **Figure S10a.** Draw the human distances at the confirmed pixel scale once item 1 is resolved, or
    state the units.

## E. Resolved since the earlier list

- **Notebook partitions.** Concordance between notebooks 03 and 13 is in Methods 6.7 (ARI 0.979 / 0.981;
  98.5% of S1/S2/S3 labels).
- **Glomerulus provenance** is in Methods 6.16 (22 of 1,581).
- **The DPT comparison** on notebook 03's coordinate is labelled historical, and DPT13 is run through
  notebook 37 (D4).
- **Census extraction scripts** are committed (`analysis/scripts/fetch_census_pt_segments.py`,
  `fetch_census_pt_cells.py`).
- **Notebooks 35–46** are committed.
- **The amended specificity rule** (5% condition) and the two kinds of q-value are disclosed (Methods
  6.14, Discussion decisions).
- **[PENDING D4] and [PENDING D8]** are resolved in v0.3.
- **Lake injury-marker list** is in Methods 6.19.
- **Zero-bin human polygons** are explained (Methods 6.4).
- **Mouse snRNA confirmation** uses 12 male donors (r3_draft corrected).
- **Previously open [CITE] items** in the Introduction are closed (Maass et al. 2019; Heyman et al.
  2010).

Still from the earlier list:
- the "13-cluster" wording in AGENTS.md and `docs/workflows/pseudospace.md`; the cross-species
  partitions have 11 clusters;
- the notebook 01 working copy's commented machine path, which must not be staged.
