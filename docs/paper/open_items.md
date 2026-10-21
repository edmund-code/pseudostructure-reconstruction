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

18. **Transporter count — CLOSED 2026-10-20.** 33 of 43 zonated transporters conserved (0.767) against
    196 of 388 zonated genes (0.505). It is **not** a hand count: notebook 47 writes it to
    `results/pt_revision_addenda/tables/transporter_class_shares.csv` (row `Slc/Abc/Aqp`), with
    `transporter_direction_agreement.csv` adding 40 zonated in both species and direction agreement
    0.825 against 0.723 for all genes, and it reproduces exactly from notebook 43's
    `gene_classes_symmetric.csv`. `draft_v03_changelog.md` §6.1 and `referee_response_v01.md` M8 still
    say no notebook saves it and are stale. Item 9's 30 of 54 is the notebook 40 count; R3 gives both.
19. **Abedini PMID.** Fixed 2026-10-04: the local `SOURCE.txt` now cites PMID 39048792 (Nat Genet
    56:1712–1724); no committed document carried the wrong PMID.
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

## C′. Added 2026-10-20 by the continuous-axis species round

Full record in `docs/paper/hvm-continuous-discovery-round.md`. These are scientific blockers, not bookkeeping.

31. **The human arm's gradient attenuation — the highest-value item in this file. Narrowed 2026-10-21 from three
    candidate causes to two.** Median human/mouse curve amplitude 0.633 on the 95-gene core and **0.453** on a
    1,429-gene clean set genome-wide, 75.8% of core genes flatter, and `criteria_hvm.csv` C12 ("A3: no
    attenuation flag") recorded **not met** at fraction 0.829. A blurred or mis-registered human axis is now
    **excluded**: a positional smear can only lower an amplitude and must attenuate rough curves more, and
    `analysis/scripts/attenuation_mechanism_human_vs_mouse.py` shows no smear width reaches the observed median,
    every width gives the opposite-sign roughness dependence (−0.57 to −0.69 against an observed +0.302, +0.283
    with independent noise in the stratifier), and the simulated sharpened share is 0.1–4.8% against 13.1%
    observed. The flatness also reproduces across the two human sections (flatter in all four pairings for 0.997
    of the most-attenuated quartile). What remains is **a per-gene capture difference between the two species'
    probe panels** versus **genuinely shallower human zonation** — and note that the centre of mass's algebraic
    offset-invariance does not remove a per-gene per-species *scale* difference, which is the surviving
    candidate. Biology and measurement are still perfectly confounded with species in a 2-vs-2 design, and both
    arms of both balanced relabelling refits carry the asymmetry identically, so no internal control separates
    the remaining two. **Closed by:**
    RNAscope on one human block measuring S1-vs-S3 fold change for Nt5e, Miox and Slc34a3 against the mouse values.
    Until it is closed, no species magnitude should be published — and the AKI axis shares the human-side
    machinery, so it is not a species-only item.
32. **No positional positive control in the human arm.** Megalin and cubilin have essentially no axial gradient in
    our human sections (amplitudes 0.170 and 0.197, the lowest checked; LRP2's two human sections anti-correlate at
    ρ −0.168), against a literature expectation of S1/S2-weighted endocytosis. **Closed by:** megalin IHC on the
    same block — which also answers item 31 in one stain.
33. **The reversal wording must drop "monotone".** Only 1 of the 14 called reversals has both limbs at monotonicity
    ≥ 0.80, and it is the gene to drop. Affects any draft text and figure legend naming monotone gradient reversal.
34. **Two coordinate-construction genes were carried as reversal evidence.** Miox and Cox7a1 are
    `fold-protected = False` in `d_end_contrasts.csv.gz` (only 5.7% of genes are), i.e. they helped fit the axis
    they are being positioned on. Decide whether to exclude unprotected genes from every species claim by rule, and
    state the rule in Methods.
35. **No sentence may name a human PT segment.** Human PT-S2 has segment length 0.0609 in HUK1_COR1 and **−0.0769**
    in HUK1_MED1, with S2→S3 (0.5895) occurring before S1→S2 (0.6665); human S1 spans 0.543–0.666 of the axis
    against mouse 0.342/0.339. Sweep the drafts for "human S2" and "S3-ward in human".
36. **The collaborators' SGLT2→SGLT1 question cannot be answered with this panel.** Slc5a1 has detection_human
    0.061 and is a mouse neighbour marker (ratio 7.12); Slc5a2–Slc2a2 is flat. Tell them, and scope any SGLT2i work
    to an assay that detects SGLT1 in human.

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
