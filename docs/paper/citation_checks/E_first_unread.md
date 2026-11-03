# E: full-text check of the five unread references

Manuscript: `docs/paper/manuscript_v1_intro_methods.md` (not edited).
Downloads: `scratchpad/fulltext/dl_E/`. Extracted text: `txt_E/`. Page images: `img_E/`.
No personal identifier was sent to any service. Unpaywall was not used. No proof-of-work or anti-bot page was solved.
Every copy that was read came from an Internet Archive (Wayback Machine) snapshot of a publisher or author file. The one exception is a direct repository download that turned out to be a preprint.

## Summary

| Ref | Full text read | Verdict | Fix |
|---|---|---|---|
| [6] Moor 2018 *Cell* | **Y**, published version (28-page Cell PDF incl. STAR Methods and supplementary figure legends) | SUPPORTED | None to the sentence. Change the tag to "[published PDF read]". |
| [60] Kanehisa & Goto 2000 *NAR* | **Y**, published version (OUP full-text HTML: all sections, both figure legends, both tables; figure images not viewed) | SUPPORTED, as the canonical KEGG reference. The 2019 library content rests on [61]. | None. Keep [60,61] together. Change the tag. |
| [66] Mao 2015 SimplePPT | **Y**, published proceedings PDF (9 pp.). The appendices A–D it mentions are not part of the proceedings paper. | SUPPORTED. scFates source confirms `method='ppt'` is SimplePPT. | None required. Optional: one-line convergence note (below). |
| [45] Hastie & Stuetzle 1989 *JASA* | **Y**, all 15 pp. of the publisher PDF, read as page images | SUPPORTED for the principal-curve concept and the arc-length projection index. PARTLY if [45] is read as the source of the elastic energy or of the end extension. | Move [45] onto the projection sentence (wording below). |
| [33] Haghverdi 2016 *Nat. Methods* | **Y for the accepted author manuscript**, peer-reviewed and pre-copyedit, with Online Methods and Supplementary Notes. The copyedited version of record was **not** read (paywalled). bioRxiv v1 was also read. | SUPPORTED (Intro and Methods) | None required. Optional precision about Scanpy's approximation (below). Change the tag. |

---

## [6] Moor, A. E. et al. *Cell* 175, 1156–1167.e15 (2018)

**Access route.** OpenAlex marks the article as bronze open access at `cell.com/article/S0092867418311644/pdf`, but curl gets HTTP 403 there. The copy that was read is the Wayback snapshot of `www.cell.com/cell/pdf/S0092-8674(18)31164-4.pdf` (2022-12-09), fetched raw with `id_` (`dl_E/wb_moor_cell.pdf`). It is the published PDF: 28 pp. with Highlights, main text, STAR Methods (e1–e5) and supplementary figure legends S1–S7.
Routes that failed: bioRxiv `api.biorxiv.org` returned the metadata (preprint 10.1101/261529, CC BY-NC), but the JATS and PDF on `www.biorxiv.org` gave HTTP 429. Europe PMC has no PMCID, and its preprint record PPR19042 `fullTextXML` gave HTTP 500. Semantic Scholar gave 429. The lab site and CORE were not needed.

**Full text read:** Y, the published version.

**Citing sentence (Intro):** "Many organs are built from repeated units whose function changes continuously along an axis, such as the liver lobule[5], the intestinal villus[6] and the kidney nephron[7,8]."

**Verdict: SUPPORTED** for the villus.
- Summary: "The intestinal epithelium is a highly structured tissue composed of repeating crypt-villus units."
- Summary: "… utilized it for single-cell spatial reconstruction, uncovering broad zonation of enterocyte function along the villus."
- Discussion: "enterocytes are not terminally differentiated cells but, rather, continuously transdifferentiate as they migrate along the villus axis"

**Notes worth knowing**
- The coordinate depends on landmark genes. Each cell gets x = (summed top landmarks) / (summed top + bottom landmarks) (STAR Methods, "Zonation reconstruction algorithm"). The landmarks come from laser-capture microdissection of five villus segments, and the cells are then binned into six zones. Introduction: "This approach is infeasible, however, when no prior knowledge exists regarding zonated landmark genes."
- On the villus, position also tracks cell age. Introduction: "The positions of enterocytes along the villus axis correlate with their age".
- Discussion: "the single-cell spatial reconstruction assigns every cell a continuous coordinate (Figure S1E)", although genes were summarized in six zones.

**Proposed fix:** none to the sentence. Change the reference tag from "[inaccessible (abstract only)]" to "[published PDF read]".

**Writing ideas**
- "Landmark-based reconstructions place each cell by the ratio of zonated landmark genes measured by microdissection[6]; pseudostructure reconstruction needs no landmark panel, only an orientation."
- Discussion contrast: on the villus, position along the axis is also cell age[6]. Proximal-tubule cells do not migrate along the nephron, so t is a purely spatial coordinate. This needs a kidney reference before it is used.

---

## [60] Kanehisa, M. & Goto, S. *Nucleic Acids Res.* 28, 27–30 (2000)

**Access route.** The copy that was read is the Wayback snapshot of `academic.oup.com/nar/article/28/1/27/2384332` (2026-09-24, gzip HTML, `dl_E/wb_kegg_oup.html`). It contains the full article text: Abstract, Introduction, Genomic information, Pathway information, Availability, the Figure 1–2 legends and Tables 1–2.
Routes that failed:
- Europe PMC lists PMC102409. Its `fullTextXML` returned HTTP 500.
- NCBI efetch has no body, only a `pmc-pdf` link (`gkd027.pdf`).
- The PMC PDF URL serves a JavaScript proof-of-work page, which was not solved.
- europepmc.org PDF render and `ptpmcrender`: HTTP 403.
- OUP article-pdf: 403.
- The PMC OA service URL returned 404.

**Full text read:** Y, the published version. The figure images were not viewed. The two figures are schematics whose legends were read.

**Citing sentence (Methods, Pathways):** "Gene sets came from MSigDB Hallmark 2020[58], Reactome 2022[59] and KEGG 2019 Mouse[60,61], as distributed by Enrichr[62,63] …"

**Verdict: SUPPORTED** as the canonical database citation. It cannot document the 2019 content, so [61] stays.
- Abstract: "KEGG (Kyoto Encyclopedia of Genes and Genomes) is a knowledge base for systematic analysis of gene functions, linking genomic information with higher order functional information."
- Abstract: "… the PATHWAY database, which contains graphical representations of cellular processes, such as metabolism, membrane transport, signal transduction and cell cycle."
- Table 1 title: "The summary of KEGG release 12.0 (October 1999)". The paper therefore describes release 12.0, not the 2019 library.

**Detail.** In 2000, organism-specific metabolic pathways were "constructed computationally by correlating genes in the genome with gene products (enzymes) in the reference pathways according to the matching EC numbers". Regulatory pathways were drawn per organism. Ortholog identifiers were still being introduced.

**Proposed fix:** none. Keep [60,61] together and change the tag from "[inaccessible (abstract only)]" to "[published full text read]".

**Writing ideas**
- "KEGG pathways for mouse (Enrichr library KEGG_2019_Mouse[60–63])". Naming the library file makes the release explicit.
- Note that KEGG mouse pathways are organism-specific maps derived from reference pathways[60]. This is why no ortholog mapping was applied to that library, unlike the human-symbol libraries.

---

## [66] Mao, Q., Yang, L., Wang, L., Goodison, S. & Sun, Y. SimplePPT. *Proc. SIAM SDM 2015*, 792–800

**Access route.** The copy that was read is the Wayback snapshot of `epubs.siam.org/doi/pdf/10.1137/1.9781611974010.89` (2023-03-15, `dl_E/wb_simpleppt.pdf`). It is the SIAM publisher PDF, 9 pp., with Abstract, §§1–5 and the references. The appendices A–D it cites (proofs, gap statistic) are not in the proceedings paper and were not read.
Routes that failed: OpenAlex lists it as closed with no repository copy. Semantic Scholar gave 429. The author pages and CiteSeerX were not needed.

**Full text read:** Y, the published version, without the separate appendices.

**Citing sentence (Methods, other methods):** "… scFates 1.2.5[65] (SimplePPT principal tree[66]) …"

**Verdict: SUPPORTED.**
- §3.2: "The proposed methods is named as SimplePPT, the pseudo-code of which is given in Algorithm 1."
- Abstract: "As a special case, a principal tree model is proposed and a new algorithm is developed that learns a tree structure automatically from data."
- §4.2: "λ is an important parameter that controls the tradeoff between the curve fitting error and the length of a principal tree."

**scFates check.** The installed scFates is in `kidney-pseudospace/site-packages/scFates/tools/graph_fitting.py`.
- Its docstring reads: "If `method=='ppt'`, uses simpleppt implementation from [Soldatov19]_."
- It also reads: "`ppt_lambda` … penalty for the tree length [Mao15]_" and "`ppt_sigma` Regularization parameter for simpleppt [Mao15]_".
- The `simpleppt` 1.1.7 loop (`simpleppt/ppt.py`) follows Algorithm 1: soft assignment R ∝ exp(−d/σ), a minimum spanning tree over the principal points, L = D − B, and F = XR(λL + Λ)⁻¹.
- One small difference from the paper: the code uses the plain distance (the `metric`, Euclidean here) inside the exponent, while the paper uses ‖x − f‖².

**Notes**
- The paper calls the algorithm "parameter-free with guaranteed convergence", yet §4.1 lists "four parameters that need to be specified" (M, the initialization, σ, λ).
- Our run used 50 nodes, σ from `explore_sigma`, λ = 100 and 50 steps. `methods_notes.md` records "not converged" on both datasets.

**Proposed fix:** none required; [65] + [66] is right. Optional Supplementary line: "SimplePPT did not reach its convergence tolerance within 50 iterations (relative node change 0.035 and 0.082)."

**Writing ideas**
- "SimplePPT re-estimates the tree as a minimum spanning tree over its principal points at every iteration[66], so it can branch on a single path. We used scFates' pseudotime along the tree as given."
- Contrast for our method: "the elastic curve fixes a chain topology; SimplePPT learns topology and node positions jointly[66]."

---

## [45] Hastie, T. & Stuetzle, W. Principal curves. *J. Am. Stat. Assoc.* 84, 502–516 (1989)

**Access route.** The copy that was read is the Wayback snapshot of `web.stanford.edu/~hastie/Papers/Principal_Curves.pdf` from 2019-01-07 (`dl_E/wb_hastie_2019.pdf`). It is the JSTOR/Atypon publisher PDF, 15 pp., with an OCR text layer that scrambles the two columns. All 15 pages were rendered at 110 dpi (`img_E/hastie-01…15.png`) and read as images.
Routes that failed:
- The 2020–2021 snapshots are truncated at about 1 MiB; only the 2019 snapshot is complete (2.2 MB).
- The live Stanford URL now redirects (302).
- OpenAlex lists it as closed. Semantic Scholar returned 404 for the DOI.
- The UW scan used earlier (`dl_C/hastie_uw.pdf`) has no text layer and was not needed.

**Full text read:** Y, all pages: definitions, algorithm, the data-set algorithm, bias, examples, principal surfaces and the Appendix.

**Citing sentences**
- Intro: "… fits a single elastic principal curve[45–47] through all cross-sections, placing each at a position t …"
- Methods: "An elastic principal curve[45–47] was fitted … Each unit was projected onto its nearest edge, and its position is the path length from the start of the curve to its projection; units beyond an end node were placed on the straight continuation of the end edge, up to twice that edge's length."

**Verdict: SUPPORTED** for the principal-curve concept and for position as arc length to the projection.
- Abstract: "Principal curves are smooth one-dimensional curves that pass through the middle of a p-dimensional data set, providing a nonlinear summary of the data."
- §2.2: "The projection index λf(x) of x is the value of λ for which f(λ) is closest to x."
- §5 (data sets): "we use the arc-length parameterization, for which λ1 = 0 and λi is the arc length along the polygon from f1 to fi."
- §5.1 (projection step): the paper defines d_ik as the distance from x_i to its closest point on each line segment of the polygon, then "we replace λi by the arc length from f1(j) to fi(j)". This is our "nearest edge, then path length".
- §7.1 (SLAC): coordinates included "the arc length from the beginning of the arc till the point of projection onto the ideal curve".

**PARTLY SUPPORTED** for two other readings:
1. **The elastic energy (fit + stretching + bending) is not from [45].** Hastie and Stuetzle define a principal curve by self-consistency, "E(X | λf(X) = λ) = f(λ) for a.e. λ" (Definition 1), and fit it by local averaging or smoothing splines (§5.2, §5.5). The energy U comes from [46,47].
2. **Extending units past the end nodes is not in [45].** Their data algorithm projects only onto segments between consecutive vertices, and they warn that "Discontinuities can certainly occur at the endpoints of a curve." (§4).

**Further caveat for the identification-limits text.** §2.2 and §6 show that a principal curve is biased for a curved generating path: "the principal curve is biased for the functional model". The bias goes to 0 as the noise variance becomes small relative to the radius of curvature.

**Proposed fix (Methods):** "An elastic principal curve[46,47], a regularized form of the principal curve[45], was fitted … Each unit was projected onto its nearest edge, and its position is the arc length from the start of the curve to that projection (the projection index of ref. 45); units beyond an end node were placed, by our convention, on the straight continuation of the end edge …". The Intro "[45–47]" can stay.

**Writing ideas**
- "t is each cross-section's projection index[45]: the arc length along the curve from the root end to the cross-section's nearest point on it."
- Identification-limits sentence: "a principal curve through noisy data around a curved path is biased toward the outside of the bend[45]. Positions are therefore defined up to this smoothing bias, which shrinks as noise falls relative to curvature."

---

## [33] Haghverdi, L., Büttner, M., Wolf, F. A., Buettner, F. & Theis, F. J. *Nat. Methods* 13, 845–848 (2016)

**Access routes**
- **Read:** Laleh Haghverdi's TUM dissertation (2016, mediaTUM 1325451). The Wayback snapshot of `mediatum.ub.tum.de/doc/1325451/1325451.pdf` (2025-08-09, 138 pp.) was saved as `dl_E/haghverdi_thesis.pdf`.
  - Its **Appendix B** is headed: "This is a pre-copyedited, author-produced PDF of an article accepted for publication in Nature Methods following peer review."
  - It contains the main text, figure legends, Online Methods and the Supplementary Notes 1–7.
  - The extract is `txt_E/dpt_accepted_main.txt`; the full text is in `txt_E/thesis.txt`.
- **Also read:** Helmholtz Munich PuSH `deliver.php?id=21330`. It is labelled "Postprint", but the file is the **bioRxiv v1 preprint** (2016-02-29) (`dl_E/dpt_hmgu_postprint.bin`).
- **Failed:**
  - nature.com is paywalled (`"isAccessibleForFree": false`); the page shows only the abstract and the figure and supplementary-figure legends.
  - MDC `edoc 19027oa.pdf` (accepted version) returned HTTP 401. MDC `edoc 19019` holds metadata only.
  - Europe PMC has no PMCID or author manuscript.
  - Semantic Scholar returned 429. The fatcat API gave no response. scholar.archive.org shows an anti-scrape check, and mediaTUM is behind an Anubis proof-of-work page; neither was solved.

**Full text read:** Y for the accepted author manuscript, which is peer-reviewed but not copyedited. The **version of record was not read.**

**Citing sentences**
- Intro: "Trajectory inference orders cells along continuous processes from their expression profiles, usually given a start cell or cluster that sets the direction[32–36]"
- Methods: "ordered by diffusion pseudotime[33] as implemented in Scanpy 1.11.5[55] (max(30, √n) neighbours, 15 diffusion components, rooted at the same root unit)"

**Verdict: SUPPORTED** for both sentences. DPT requires a root cell.
- Main text (accepted): "Provided with a known `root cell`, we measure the developmental progression of each cell in the dataset by computing its DPT with respect to the root cell."
- Online Methods (accepted), "Overview of algorithmic steps": the "user-provided parameters" include "b. The index or indices of one or several root cells."
- Online Methods (accepted): "Consequently, DPT is based on the full rank 𝑻 rather than a low rank approximation, i.e. we do not use diffusion maps as a dimensionality reduction method …"

**Scanpy check.** In Scanpy 1.11.5, `neighbors/__init__.py::_get_dpt_row` sums (λ/(1−λ))²-weighted differences over the computed eigenpairs, i.e. the leading `n_dcs` diffusion components. `tl.dpt` stops with an error unless `adata.uns['iroot']` is set. So "15 diffusion components" truncates the paper's full-rank definition, and "as implemented in Scanpy" already makes that clear.

**Differences between the preprint and the accepted manuscript**
1. **Abstract.**
   - bioRxiv v1: "This allows us to identify cells that undergo branching decisions or are in metastable states, and thereby genes differentially regulated at these states."
   - Accepted: "This allows us to reconstruct the cells' developmental progression and to identify transient or metastable states, branching decisions or differentiation endpoints."
   - Published abstract (nature.com): "Our DPT software implementations make it possible to reconstruct the developmental progression of cells and identify transient or metastable states, branching decisions and differentiation endpoints."
2. **Box 1 ("Universal time") was removed from the main text.** Its content moved to an Online Methods section, "Relation of pseudotime in snapshot data to actual time", and to Supplementary Note 6. The sentence "Pseudotimes are proxies for universal time" survives there, so the writing idea in `B_refs_22_40.md` can still be quoted from the accepted version.
3. **The main text was rewritten.**
   - The DPT equation moved from the main text to the Online Methods.
   - The diffusion map is now introduced as "a non-linear method for recovering the low-dimensional structure underlying high-dimensional observations".
   - The root-cell sentence was added to the main text.
   - "DPT does not rely on dimension reduction" was added.
4. **Root cells.**
   - The preprint v2 phrasing was "One (or several) root cell(s)" plus "If several roots are defined, DPT averages the pseudotime … over these roots".
   - The accepted overview says "The index or indices of one or several root cells". The averaging sentence was not found in the accepted Online Methods overview.
   - This makes no difference here, because we use one root.

**Proposed fix:** none required. Optional precision for the Methods: "… diffusion pseudotime[33] as implemented in Scanpy 1.11.5[55], which computes it from the leading diffusion components (15 here; max(30, √n) neighbours), rooted at the same root unit". Change the tag to "[accepted author manuscript read (TUM thesis, Appendix B); version of record paywalled]".

**Writing ideas**
- "DPT orders each unit by its random-walk distance to the root unit[33]; like our curve, it needs only the root to set the direction."
- Borrow the accepted-version framing: "pseudotimes are proxies for universal time"[33]. Likewise, t is a proxy for position along the shared structure, not a measured distance.
