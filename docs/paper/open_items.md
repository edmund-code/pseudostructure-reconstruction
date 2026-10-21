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

## C″. Added 2026-10-21 by the AKI-axis audit and the full-text literature pass

Sources: `aki-axis-thesis-audit.md`, `literature-novelty-ledger.md`. Items 37–40 are the ones that change what the
paper may claim; 41–46 are corrections and chores.

37. **Three sentences are now forbidden, on published prior art.** "First to order unordered cross-sections"
    (Chakravorty 2026, Cell — scPrisma on 337 per-tubule aggregates, no x/y). "Generalises to any repeated
    structure such as crypts and lobules" — that is a 2022 SpaceFold quotation (Niec 2022), four years older than
    anything the project's earlier method survey found, and the word "pseudospace" is itself taken at least three
    times. "Aggregate per structure then order the aggregates." The novelty must move onto the **identification
    contract** with all five qualifiers attached (partial cross-sections at unknown position; unknown endpoints;
    specimen-dependent coverage; a single-cell-type axis; specimen-level calibration with cross-species transfer),
    and the defence is that **no published open-axis method can run on a single-cell-type epithelium** — SpaceFold's
    axis is built from cell-type composition and needs a deconvolution reference.
38. **Two experiments the field will expect, neither yet run.** (i) scPrisma's **linear** mode on our tubule
    aggregates as a named comparator to the elastic curve. (ii) The physical validation this class of method is now
    judged by: inferred position against each cross-section's **own measured diameter or area** — not cortical
    depth, which this project has already ruled out as a selector. Chakravorty validates at ρ 0.94 against a
    measured physical dimension; a reviewer who knows that paper will ask.
39. **The highest-value single experiment from the whole round: stain DPP4 in matched human and mouse kidney
    cortex.** Dpp4 is the top gene of the sharper-in-human set, it survives the admixture screen as an anti-rider
    (ρ_distal −0.579/−0.598), and its species asymmetry is already published in human-versus-mouse **liver** with
    immunofluorescence validation (Yakubovsky 2026, Nature — verified in full text). Antibody staining carries no
    probe panel, so a positive result removes the per-gene capture explanation for the entire 121-gene set in one
    slide, and the prediction is pre-registered in someone else's figure. This now outranks item 31's RNAscope
    panel on cost-to-information.
40. **A third candidate cause for the attenuation, testable on our own data.** Item 31 records two surviving
    causes. The literature supplies a third that is distinct from the positional smearing already excluded:
    a **composition** difference, i.e. fewer pure-S3 units in the human arm. Test internally — condition the
    amplitude ratio on per-specimen axial coverage and pure-S3 unit counts from `composition.csv` and report
    whether it moves. Until then item 31 should read *three* causes, not two.
41. **Reframe DCXR and SLC17A3 as validations, in those words.** Both have published human limbs that agree with
    us (McEvoy 2022 Fig. 1e for DCXR in 19 donors; Sakaguchi 2024 Table 1 for SLC17A3 at 6.02/3.82/3.56). DCXR is
    the paper's strongest validation — two independent designs, neither carrying our probe panels or our axis, get
    the same answer in both species — and only its within-pure-S1 and within-pure-S3 half is new. Published as a
    discovery it is the claim that costs the paper credibility.
42. **Delete the PMID 18981171 attribution wherever it appears.** It is a human cloning and enzymology paper with
    no mouse data and no intra-renal localisation; it cannot support the UGT3A1 species framing. Already struck in
    `hvm-continuous-discovery-round.md`; check the figure legends and any draft text.
43. **Stop calling ABHD14A a reversal** — the published rodent record has mouse flat, so the claim is "no mouse
    gradient, a strong human gradient", the wording already used for NT5E. And **drop the ABHD17B/ASCT2
    mechanism**: it requires apical ASCT2 in the mouse straight segment and the published rodent data have none.
44. **Three abstract-only sources are load-bearing and must go through the library** before the species sections
    are written: **Vrhovac 2015** (the only human measurement behind the S1/S2-versus-S3 SGLT2/SGLT1 scheme),
    **Zhai 2003** (mouse PT ultrastructure, the classical claim that runs against claim B), **Guder & Ross 1984**.
45. **Cite the literature that runs against claim B in the same breath as claim B.** Albadry 2024 measures
    cross-species gradient intensity in liver and finds human similar to **mouse** with rat as the flat outlier —
    the opposite direction for our species pair. Claim B's novelty is intact (no publication in any organ reports a
    species-systematic axial amplitude difference, and the best-placed paper compares direction only, R = 0.15),
    but it must not be presented as unopposed.
46. **Corrections found by the AKI audit.** The constant-within-label control at 4.6%, on which the surviving
    change-density claim rests, has **not** been re-derived — do that before publishing the 6–20× ratio.
    `fold-protected` has **opposite polarity** between `a_gene_classes.csv.gz` (True for 93.8%) and
    `d_end_contrasts.csv.gz` (unprotected ~5.7%); rename one. `pt-continuous-axis-experiments.md`'s "1,009 /
    141 replicated / 115 fold-protected" should read **145 / 118**; and it names Slc34a1 and Slc5a2 as
    "continuous shapes favoured" when the lowest-BIC model of six for both is a **step** model.

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

47. **M6's convergence ladder needs a new primary null, and the AKI half of it needs withdrawing.**
    Re-run 2026-10-21 with the project's own `stratified_line_null` and one change — strata of expression decile
    **× |amplitude| quartile** instead of decile alone, because a large-effect gene is more likely to be called by
    every line and the published null does not model that. The AKI ladder collapses (5.08× → **1.63×** at k ≥ 6;
    the frozen k ≥ 3 bar sits at **1.04×**), so withdraw the AKI multiples. The species ladder survives with
    deflated multiples (k ≥ 4 **1.91×**, k ≥ 5 **3.21×**, k ≥ 6 **16.7×**) and should be reported with the
    amplitude-preserving null as primary and the decile-only as a sensitivity. Both atlases are *below* 1.0 at
    k ≥ 1 and k ≥ 2 under either null — lines overlap less than chance at the bottom — which deserves its own
    sentence, being the opposite of what a shared-artefact story predicts.
48. **Two lines cannot be both evidence and anchor; fix the bookkeeping.** Verified as exact column identities:
    `call S7 handoff` and `is handoff` are the **same 58-gene column**, so an "anchored share" on a core counting
    S7 is partly tautological (S7 carries 7 of the 24 species genes at k ≥ 5); all 195 `call A9 androgen` calls lie
    inside the 566-gene `is androgen set`; and `supported` is **exactly** `call A5 | call A6`, two of the nine
    counted AKI lines. Drop S7 from the counted lines or stop calling handoff membership an anchor, and say that
    `supported` overlaps the line count.
49. **A package-era species ladder has six lines, not seven.** The atlas-level line **S4** is an external-atlas
    contrast and cannot be rebuilt package-era, which is the answer to the open question recorded in
    `human-mouse-paper-ready.md` §6. Replacements exist on disk for the other gaps (level line in
    `s3_level_gene_table.csv`; F-test in `gene_statistics_human_vs_mouse.csv`, 3,584 of 11,055 against notebook
    67's SCF13-era 3,889 at Jaccard 0.62; zonation class in `pt_package_label_comparison/.../43/`).

## C‴. Added 2026-10-21 by hypothesis round 2

50. **Build the PT-identity restriction for the human arm.** It is the filter that decides the gene list, and it is
    the one instrument that separates "the compartment changed" from "the cell type changed": it flipped Hao2's
    mouse sign (retention 0.047, so its mouse proximal weighting was thin-limb units inside the PT-S3 label) and
    it is the reason to trust Pck1 (0.808), Dcxr (0.861) and Kif12 (0.943). `post_hoc_unit_identity_by_bin.csv`
    carries it for the two mice only (PT identity 0.87–1.00 per bin); the human side needs building. Until it
    exists, no sharper-in-human or reversal list is final. This supersedes the marker-score admixture screen, whose
    rank correlation partly re-encodes axial direction (corr with ρ_bin +0.686/+0.697) and whose admixture-specific
    enrichment is **2.13×**, not the 6.3× first recorded.
51. **Ugt3a1 comes off the reversal list.** It is the only many:many ortholog in the set, its mouse paralog Ugt3a2
    is absent from the shared panel so paralog substitution cannot be excluded internally, 94% of its human effect
    is thick-limb-attributable, the sign flips to +0.392 extrapolated to zero differential, and its human distal
    bins hold 5–24 counts. The project had it ranked first. Also off: Hao2 (and the oxalate/glyoxylate trio),
    Abhd14a (detection 0.256–0.261), Abhd17b (0.351/0.351), Slc26a1/Sat1, the Acy1⇄Acy3 baton pass, and
    Slc6a18→Slc6a19 as a species handoff.
52. **The reversal screen's operating point is τ = 0.20, not 0.15.** Re-derived independently: at τ = 0.15 the
    human-section within-arm null still fires twice on the raw statistic, and is empty only after a filter stack
    containing an external-atlas filter. At **τ = 0.20 both within-arm nulls are empty with no filtering**, giving
    **47 genes** on internal-only filters, **47 of 47** passing an orientation-free re-check (Spearman against bin
    index over all 20 bins, no purity selection, opposite-signed in all four specimens). Use that pair of numbers.
53. **Put the human-replication fact in the limitations verbatim.** The two human sections are two sections of one
    block and they agree **worse than two separate animals**: 44 gene pairs differ between them at |Δz| ≥ 0.3
    against 4 between the two mice. Every "4-of-4 agreement" in this project is 2 biological mouse replicates plus
    2 technical human sections.
54. **Tell the collaborators the SGLT2i answer, because a powered null is more useful than silence.** The testable
    unit SGLT2 → GLUT2 is **identically positioned** in the two species: matching index −0.393 mouse against
    −0.396 human, ΔD = −0.003, against an amplitude surrogate generating a median |ΔD| of 0.383 with 48.6% of
    4,000 draws exceeding the observation; the two amplitude ratios are 0.533 and 0.531. So no positional
    explanation for species differences in SGLT2i response is supported here, the distal escape route is
    unmeasured, and the assay that would answer it must detect SGLT1 in human PT (RNAscope or IHC, not this panel),
    scored on straight versus convoluted profiles.
55. **Resolve PCK1's mouse limb before publishing it as a reversal.** Our mouse Pck1 falls toward S3
    (E −0.354/−0.255) while the published rodent record reports PEPCK *highest* in S3 (Drewnowski 2002,
    **abstract only**, so load-bearing-on-abstract). Human rises strongly (+0.743/+0.864) at detection 0.98/0.99,
    the best-detected gene in the list, with PT-identity retention 0.808 — so the claim is worth the slide, and the
    stain already proposed resolves both questions at once.
56. **Examine the NAD branch point; it is the round's most interesting unexamined result.** QPRT reverses
    independently (mouse +0.127/+0.185 → human −0.296/−0.434 on PT-identity units, and also on the committed
    fitted statistic, detection 0.656/0.567, balanced probes) and neither track promoted it. With HAAO conserved
    distal in both arms and ACMSD reversed, three genes of the branch point now point the same way.
57. **KIF12 is the best new reversal and has a named upstream cause.** Mouse −0.471/−0.503 against human
    +0.384/+0.435, detection 0.87/0.86, PT-identity retention 0.943, and a direct HNF1B transcriptional target —
    which is testable, unlike most of the list.

## C⁗. Added 2026-10-21 by the drug-geography and figure round

58. **Repair the G6pc / G6PC1 pairing — the terminal gluconeogenic step is invisible to the project.** Mouse
    `G6pc` and human `G6PC1` are the same 1:1 ortholog pair in the same inventory component (6169), both well
    detected (mouse 0.555/0.608, human 0.987/0.952), yet `G6pc` is **absent from the 11,055-gene shared space**
    under either symbol, so glucose-6-phosphatase cannot appear in any cross-species table. It is also the
    steepest gluconeogenic gene in the mouse arm (−3.2/−3.5). Fix before the glucose story is written. Not
    repairable, and to be stated as panel asymmetry: SLC9A3R1, LDHA, ALDOA, TPI1 are mouse-only.
59. **Correct three detection facts in the project record.** `Slc2a9` fails in the **mouse** arm (0.082; human
    0.802) — the record had it as a human failure. `Abcc4` fails in mouse (0.089; human 0.533). `Slc5a8` fails in
    human (0.071; mouse 0.851), so the apparently conserved SMCT1 contrast must not be used. Also `Oxct1` human
    detection is **0.157** with the committed neighbour flag True (ratio 5.03): report its human limb as **not
    measured**, not as flat, and drop the SLC9A3/NHE3 reversal for the same reason.
60. **Report 18, not 16, for the opposite-sign share of the 29 movers.** 18 is the count on the **pure-compartment**
    contrast per specimen (the retained statistic); 16 is the same count on the **curve-derived** arm contrasts
    (the retracted one). They agree on 25 of 29; the four that differ are Hao2, Chp1, Msrb1 and Txnip. Classes on
    the retained statistic: 18 opposite sign / 6 same sign > 1.5× apart / 4 same sign within 1.5× / 1 specimens
    disagree within an arm.
61. **The reversal phenomenon does not need the reconstructed axis, and that should be said.** An axis-free
    contrast — `log2((S3+0.1)/(S1+0.1))` from the committed per-specimen `segment_means`, i.e. H&E segment labels
    and raw counts only, no coordinate, spline, registration or purity bins — gives within-arm nulls of **0 of
    2,946** (mouse) and **14 of 2,664** (human sections) against **380 of 915** cross-species over 10,669 genes.
    The axis is how the genes were found and how shape is quantified; it is not what makes them visible. This is
    the strongest available answer to "is the coordinate inventing the biology".
62. **Publish the endocytic negative.** LRP2, CUBN, AMN and PDZK1 have **no resolvable gradient in either arm**.
    A nephrotoxicity reader will ask where aminoglycoside uptake sits, and the honest answer here is that this
    cohort cannot place it. Say so rather than omitting the panel.
63. **The PCK1 dissociation is the round's strongest biology and needs its own paragraph.** Mouse gluconeogenesis
    is one proximal block co-located with SGLT2 (every step S1-weighted); in human PCK1 inverts (axis-free
    +1.21/+1.43 against mouse −0.73/−0.64) while FBP1, ALDOB, GPT2, GOT1 and GOT2 stay proximal. A dissociation
    inside one pathway in one assay cannot be produced by global flattening or by a registration error. Its
    hazards are already recorded as items 55 (the rodent PEPCK direction) and the human-distal admixture gate,
    which PCK1 clears with an atlas ratio of 0.19.
64. **Fix the SMCT naming in nb92's `h1_gene_shapes.csv`:** SLC5A8 is SMCT1, SLC5A12 is SMCT2; the table has them
    swapped.
