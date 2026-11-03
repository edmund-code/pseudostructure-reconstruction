# C · Methods references: full-text check (refs 41–47 and the Methods placeholders)

Checked 2026-11-03 against `docs/paper/manuscript_v1_intro_methods.md` (Methods) and `claims.md`. The downloads are
in `fulltext/dl_C/`, the extracted text in `fulltext/txt_C/`, and the page renders in `fulltext/img_C/`. Parsing
scripts are in `fulltext/scripts_C/` and were run with `python -I`. Where a sentence turned on what we actually ran,
I checked the code: the installed `pseudospace_reconstruction` (elastic curve, Harmony call), `segmentation/src`
(encoder, PQ, watershed), `pseudospace/nb_gam.py` and `stats_gam.py` (gene models), `pseudospace/cross_species.py`
(orthologs), `data/external/hgnc/SOURCE.txt` and `data/external/census_pt_segments/SOURCE.txt`.

**Verdicts:** SUPPORTED means the full text backs the specific technical statement. PARTLY means the citation is
right but one stated detail is not in that paper or differs from it. A row marked NEW is a reference I propose
adding.

## Summary table

| Ref | Access route | Full text | Verdict | Proposed fix | Full citation (DOI) |
|---|---|---|---|---|---|
| [41] DINOv2 | arXiv 2304.07193v2 PDF (TMLR camera-ready, header "Published in TMLR (01/2024)") | Y | PARTLY | ViT-g/14, 1.1B and the 14-px patch are supported. **Registers are not in DINOv2**: add Darcet 2024 next to "with registers". | Oquab, M. et al. DINOv2: learning robust visual features without supervision. *Trans. Mach. Learn. Res.* (2024). https://openreview.net/forum?id=a68SUt6zFt (arXiv:2304.07193) |
| NEW Darcet 2024 (registers) | arXiv 2309.16588v2 PDF (ICLR 2024) | Y | SUPPORTED | Cite for "with registers". | Darcet, T., Oquab, M., Mairal, J. & Bojanowski, P. Vision transformers need registers. In *Int. Conf. Learn. Represent. (ICLR)* (2024). arXiv:2309.16588 |
| [42] Midnight | arXiv 2504.05186 PDF | Y | SUPPORTED | None. Midnight itself does not mention registers, so do not cite [42] for that word. | Karasikov, M. et al. Training state-of-the-art pathology foundation models with orders of magnitude less data. In *MICCAI 2025*, LNCS 15967, 573–583 (Springer, 2025). doi:10.1007/978-3-032-04984-1_55 |
| [43] OpenMidnight | Zenodo API record + Zenodo PDF (SOPHONT-TR-2025-001 v1.0) + GitHub README | Y | SUPPORTED | "Open replication of Midnight" is the authors' own wording. Registers and 224-px crops are stated only in the README. The report says "1.1B"; our 1.13 billion is the count of the released checkpoint as loaded (1,134,777,344), so either write "about 1.1 billion" or add "(1.13 billion in the released checkpoint)". Consider adding the training data (about 12,000 TCGA H&E slides). | Kaplan, D., Grandhi, R. S., Lane, C., Warner, B., Abraham, T. M. & Scotti, P. S. How to train a state-of-the-art pathology foundation model with $1.6k. Technical Report SOPHONT-TR-2025-001 (Sophont, 2025). doi:10.5281/zenodo.20711012 |
| [44] Harmony | NCBI efetch JATS of PMC6884693 (Europe PMC fullTextXML returned HTTP 500) + R harmony 2.0.5 Rd | Y | SUPPORTED | θ and λ are the paper's names. "30 iterations" should be "at most 30 rounds" (`max_iter`, with early stopping on). λ = 1 is fixed explicitly: harmony 2.0.5's default `lambda = NULL` estimates λ instead. | Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* 16, 1289–1296 (2019). doi:10.1038/s41592-019-0619-0 |
| [45] Gorban & Zinovyev 2005 | arXiv cond-mat/0405648v2 PDF (journal ref. on p. 1) | Y | PARTLY | Same three-term energy, but its bending term is the second difference, 4× our μ term, and its curve grows by bisecting the **most loaded** edge, not the lowest-energy one. Keep it paired with [46] for the energy; cite [46] alone for the growth rule. | Gorban, A. N. & Zinovyev, A. Y. Elastic principal graphs and manifolds and their practical applications. *Computing* 75, 359–379 (2005). doi:10.1007/s00607-005-0122-6 |
| [46] ElPiGraph | Europe PMC fullTextXML PMC7516753 | Y | SUPPORTED | Energy, bisect-edge grammar, lowest-energy selection and the λ/μ guidance all match. The paper's default start is **two** nodes on PC1; our three-node start is elpigraph-python's `InitNodes=3` default, so attribute it to the software. The extension past end nodes is ours (ElPiGraph has its own leaf-extrapolation function). | Albergante, L. et al. Robust and scalable learning of complex intrinsic dataset geometry via ElPiGraph. *Entropy* 22, 296 (2020). doi:10.3390/e22030296 |
| Hastie & Stuetzle 1989 | Scanned JSTOR copy on the second author's UW page; pp. 502 and 504 read as rendered images (no text layer) | Partial (2 of 15 pp.) | SUPPORTED (for what is cited) | Add it at "elastic principal curve" for the principal-curve concept and the arc-length projection index behind "position = path length to the projection". | Hastie, T. & Stuetzle, W. Principal curves. *J. Am. Stat. Assoc.* 84, 502–516 (1989). doi:10.1080/01621459.1989.10478797 |
| [47] tradeSeq | Europe PMC fullTextXML PMC7058077 | Y | PARTLY | It supports NB, log(N_i) offset, K = 6 shared cubic basis and quantile knots. It differs in being penalized (GCV smoothing parameter, mgcv joint dispersion, Wald tests); ours is an **unpenalized** regression spline with ML dispersion held fixed and LR tests. Write "as in tradeSeq[47], but without a smoothing penalty", and add McCarthy 2012 for NB GLMs with an offset compared by LR between nested fits. The Introduction's "negative-binomial generalized additive models[47]" should match Methods. | Van den Berge, K. et al. Trajectory-based differential expression analysis for single-cell sequencing data. *Nat. Commun.* 11, 1201 (2020). doi:10.1038/s41467-020-14766-3 |
| NEW McCarthy 2012 | Europe PMC fullTextXML PMC3378882 | Y | SUPPORTED | Cite beside [47] in "Models". Our dispersion is plain ML (no Cox–Reid), so cite it for the GLM and LR framework, not the dispersion estimator. | McCarthy, D. J., Chen, Y. & Smyth, G. K. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Res.* 40, 4288–4297 (2012). doi:10.1093/nar/gks042 |
| U-Net (placeholder) | arXiv 1505.04597v1 PDF | Y | SUPPORTED | Fill "[ref]". | Ronneberger, O., Fischer, P. & Brox, T. U-Net: convolutional networks for biomedical image segmentation. In *MICCAI 2015*, LNCS 9351, 234–241 (Springer, 2015). doi:10.1007/978-3-319-24574-4_28 |
| Kirillov 2019 PQ (placeholder) | CVF open-access PDF (CVF pp. 9404–9413) | Y | PARTLY | PQ = SQ × RQ and unique matching at IoU > 0.5 are supported. Kirillov names the second factor **recognition quality (RQ)**; "detection quality (DQ)" is HoVer-Net's term. Kirillov averages PQ **over classes**; we select on pooled class-agnostic PQ. Write "recognition quality" (or cite Graham 2019 for DQ) and say "class-agnostic PQ". | Kirillov, A., He, K., Girshick, R., Rother, C. & Dollár, P. Panoptic segmentation. In *Proc. IEEE/CVF CVPR* 9396–9405 (2019). doi:10.1109/CVPR.2019.00963 |
| NEW (optional) HoVer-Net | arXiv 1812.06499v5 PDF | Y | SUPPORTED | Only if "detection quality" is kept. | Graham, S. et al. Hover-Net: simultaneous segmentation and classification of nuclei in multi-tissue histology images. *Med. Image Anal.* 58, 101563 (2019). doi:10.1016/j.media.2019.101563 |
| Meyer & Beucher 1990 (placeholder) | INACCESSIBLE — abstract only (ScienceDirect 403; Semantic Scholar lists it as CLOSED; abstract from the earlier notes' BibBase mirror) | N | Unverified | Prefer the 1993 chapter below, whose full text was read. | Meyer, F. & Beucher, S. Morphological segmentation. *J. Vis. Commun. Image Represent.* 1, 21–46 (1990). doi:10.1016/1047-3203(90)90014-M |
| Beucher & Meyer 1993 (alternative) | Author's lab PDF (CMM Mines Paris), pp. 433–481 | Y | SUPPORTED | Use for "marker-controlled watershed". | Beucher, S. & Meyer, F. The morphological approach to segmentation: the watershed transformation. In *Mathematical Morphology in Image Processing* (ed. Dougherty, E. R.) 433–481 (Marcel Dekker, 1993). doi:10.1201/9781482277234-12 (CRC reissue) |
| HCOP (placeholder) | Europe PMC fullTextXML PMC8574622 | Y | PARTLY | Right source, but HCOP gives support-counted assertions and does not itself label pairs one-to-one. State our rule: pairs supported by ≥ 3 HCOP resources, kept when reciprocal best by support (`cross_species.py`, `min_support=3`). | Yates, B., Gray, K. A., Jones, T. E. M. & Bruford, E. A. Updates to HCOP: the HGNC comparison of orthology predictions tool. *Brief. Bioinform.* 22, bbab155 (2021). doi:10.1093/bib/bbab155 |
| Scanpy (placeholder) | Europe PMC fullTextXML PMC5802054 | Y | SUPPORTED | Fill "[ref]". Our HVG thresholds are Scanpy 1.11.5's `seurat`-flavour defaults (its docstring cites Satija 2015; not read). | Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* 19, 15 (2018). doi:10.1186/s13059-017-1382-0 |
| scFates (placeholder) | Europe PMC fullTextXML PMC9805561 | Y | PARTLY | scFates is right, but the baseline ran `tl.tree(method='ppt')`, i.e. SimplePPT. Write "scFates 1.2.5 (SimplePPT principal tree)" and cite Mao 2015 too (Crossref only, not read). | Faure, L., Soldatov, R., Kharchenko, P. V. & Adameyko, I. scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single-cell data. *Bioinformatics* 39, btac746 (2023). doi:10.1093/bioinformatics/btac746. Add: Mao, Q., Yang, L., Wang, L., Goodison, S. & Sun, Y. SimplePPT: a simple principal tree algorithm. In *Proc. 2015 SIAM Int. Conf. Data Mining* 792–800 (2015). doi:10.1137/1.9781611974010.89 |
| [22] Leiden (Methods use) | Europe PMC fullTextXML PMC6435756 | Y | SUPPORTED | None. | Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. *Sci. Rep.* 9, 5233 (2019). doi:10.1038/s41598-019-41695-z |
| Benjamini & Hochberg 1995 | Text-layer PDF of the JRSS-B article (dcscience.net copy) | Y | SUPPORTED | None needed. The proof assumes independence; our set statistics share genes (see section). | Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* 57, 289–300 (1995). doi:10.1111/j.2517-6161.1995.tb02031.x |
| HGNC gene groups (placeholder) | Europe PMC fullTextXML PMC9825485 (Seal 2023) and PMC4739092 (Gray 2016); group membership downloaded from genenames.org | Y | SUPPORTED | Cite Seal 2023 (the current resource paper, and the one our `SOURCE.txt` already names); Gray 2016 is optional, for the group hierarchy. Name the groups: 752 "Solute carriers" (76 subgroups) and 417 "ATP binding cassette transporters" (7 subgroups); a family is an HGNC subgroup. Give the download date (2026-10-12). | Seal, R. L. et al. Genenames.org: the HGNC resources in 2023. *Nucleic Acids Res.* 51, D1003–D1009 (2023). doi:10.1093/nar/gkac888. Optional: Gray, K. A., Seal, R. L., Tweedie, S., Wright, M. W. & Bruford, E. A. A review of the new HGNC gene family resource. *Hum. Genomics* 10, 6 (2016). doi:10.1186/s40246-016-0062-6 |
| Enrichr (placeholder) | Kuleshov 2016: Europe PMC fullTextXML PMC4987924; Xie 2021: NCBI efetch PMC8152575; Enrichr site citation text | Y | SUPPORTED | Cite Xie 2021 (+ Kuleshov 2016) as the library source. Enrichr asks for Chen 2013 too (not read). Also cite each source database (below). | Xie, Z. et al. Gene set knowledge discovery with Enrichr. *Curr. Protoc.* 1, e90 (2021). doi:10.1002/cpz1.90; Kuleshov, M. V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res.* 44, W90–W97 (2016). doi:10.1093/nar/gkw377 |
| MSigDB Hallmark (placeholder) | NCBI efetch PMC4707969 | Y | SUPPORTED | Fill "[refs]". | Liberzon, A. et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst.* 1, 417–425 (2015). doi:10.1016/j.cels.2015.12.004 |
| Reactome (placeholder) | Europe PMC fullTextXML PMC8689983 | Y | SUPPORTED | Best match to Enrichr's `Reactome_2022` (1,818 terms). Reactome is human-curated; say the sets were mapped to mouse through the ortholog map. | Gillespie, M. et al. The reactome pathway knowledgebase 2022. *Nucleic Acids Res.* 50, D687–D692 (2022). doi:10.1093/nar/gkab1028 |
| KEGG (placeholder) | Kanehisa 2019: Europe PMC fullTextXML PMC6324070. Kanehisa & Goto 2000: INACCESSIBLE — abstract only (efetch had no body; OUP and Europe PMC PDF links returned HTML) | 2019 Y / 2000 N | SUPPORTED (2019) | Cite Kanehisa & Goto 2000 (canonical; on KEGG's own list) + Kanehisa 2019 NAR (matches `KEGG_2019_Mouse`). | Kanehisa, M. & Goto, S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res.* 28, 27–30 (2000). doi:10.1093/nar/28.1.27; Kanehisa, M., Sato, Y., Furumichi, M., Morishima, K. & Tanabe, M. New approach for understanding genome variations in KEGG. *Nucleic Acids Res.* 47, D590–D595 (2019). doi:10.1093/nar/gky962 |
| CZ CELLxGENE Census (placeholder) | Europe PMC fullTextXML PMC11701654 + Census docs ("Citing Census") + CELLxGENE curation API for the source dataset | Y | SUPPORTED | Cite the NAR 2025 paper **and** the source study, as the Census docs ask: Acera-Mateos et al. 2025 (bioRxiv; collection metadata only, not read). Give the Census release (2025-11-08). "Six donors" is right once a filter is stated: the dataset has 10 donors, and 6 have ≥ 50 cells in both the PCT and S3 classes. | CZI Cell Science Program et al. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *Nucleic Acids Res.* 53, D886–D900 (2025). doi:10.1093/nar/gkae1142; source data: Acera-Mateos, M. et al. Systematic evaluation of single-cell multimodal data integration for comprehensive human reference atlas. Preprint at *bioRxiv* doi:10.1101/2025.03.06.637075 (2025) |

---

## [41] Oquab et al. 2024, DINOv2 (TMLR)

**Access:** arXiv 2304.07193v2 PDF (`dl_C/dinov2.pdf`). Its first line reads "Published in Transactions on Machine Learning Research (01/2024)", so this is the TMLR version.

**Our sentences:**
- Intro: "segments tubule cross-sections from H&E histology with a frozen pathology foundation model[41–43]".
- Methods: "OpenMidnight (an open replication of Midnight, a DINOv2 ViT-g/14 with registers and 1.13 billion parameters)[41–43] … supplies detail below the 14-pixel patch size".

**Verdict: PARTLY.** DINOv2 supports ViT-g/14, about 1.1B parameters, the 14-px patch, and frozen features. The word "register" does not occur anywhere in the paper. Register tokens come from Darcet et al. 2024 (next entry).
- §5 Efficient implementation: "we use an embedding dimension of 1536 with 24 heads (64 dim/head), rather than 1408 with 16 heads (88 dim/head)."
- §5: "our ViT-g backbone counts 1.1B parameters."
- §1 Introduction: "self-supervised pretraining alone is a good candidate for learning transferable frozen features".

**Writing ideas:**
- DINOv2 states what is frozen in a single clause: "training a simple classifier over a frozen backbone, and do not perform finetuning of the backbone weights" (§7.1). We could say the same: encoder frozen, only the stem, decoder and heads trained.
- It names the crop sizes the model saw: "forwarding both large crops (at resolution 224) and small crops (resolution 98)". This is the natural place to state that our 518-px input differs from the 224-px pretraining crops.

## NEW: Darcet et al. 2024, Vision transformers need registers (ICLR 2024)

**Access:** arXiv 2309.16588v2 PDF (`dl_C/registers.pdf`), header "Published as a conference paper at ICLR 2024".

**Our sentence:** "a DINOv2 ViT-g/14 with registers" (Methods).

**Verdict: SUPPORTED.** This is the source of the register tokens that `dinov2_vitg14_reg` carries. The OpenMidnight README confirms 4 registers.
- §2.2 Hypothesis and remediation: "we explicitly add new tokens to the sequence, that the model can learn to use as registers."
- §2.2: "At the end of the vision transformer, these tokens are discarded, and the [CLS] token and patch tokens are used as image representations".

**Writing ideas:**
- Registers can be explained in one phrase. The paper says "We add these tokens after the patch embedding layer, with a learnable value", and they are discarded at the output. A decoder reads only patch tokens, so registers need no further mention.

## [42] Karasikov et al. 2025, Midnight (MICCAI 2025)

**Access:** arXiv 2504.05186 PDF (`dl_C/midnight.pdf`); Crossref confirms LNCS pp. 573–583.

**Our sentence:** as for [41].

**Verdict: SUPPORTED.** It is the recipe OpenMidnight replicates: ViT-g/14, 1.1B parameters, started from the DINOv2 checkpoint, 12k TCGA slides, KDE regularizer. It does not mention registers.
- §2 Methods, "Self-supervised training with DINOv2": "We use the DINOv2 self-distillation framework to train ViT-g14 models with 1.1B parameters".
- §2, "Training data": "TCGA contains 12k FFPE slides from 32 cancer types".
- §2, "Extraction of training tiles": "We trained our FMs on tiles of size 256×256 cropped from the original WSIs at magnifications of 2, 1, 0.5, and 0.25 µm/px."

**Writing ideas:**
- Midnight gives the pixel scale of its training tiles (0.25–2 µm/px). Our images at 0.44 and 0.5 µm/px fall inside that range, which is worth one clause, because it is why no resampling was needed. (OpenMidnight does not state its own scale; see [43].)
- Midnight describes its model as a delta from a base recipe: "Our algorithm is based on the original DINOv2 algorithm [20] with several modifications". We can do the same: OpenMidnight = Midnight recipe on DINOv2 ViT-g/14-reg.

## [43] Kaplan et al. 2025, OpenMidnight (Sophont technical report)

**Access:** Zenodo REST record 20711012 (DataCite "Report", 2025-11-14, CC BY 4.0, six creators). Also the Zenodo-hosted PDF (`dl_C/openmidnight_zenodo.pdf`, "Sophont Technical Report SOPHONT-TR-2025-001. Published November 14, 2025. Version 1.0") and the GitHub README (`dl_C/openmidnight_readme.md`). The Zenodo PDF matches the blog PDF except for one sentence.

**Our sentence:** "OpenMidnight (an open replication of Midnight, a DINOv2 ViT-g/14 with registers and 1.13 billion parameters)".

**Verdict: SUPPORTED, with two notes.**
- Abstract: "We present OpenMidnight, a replication and improvement of the Midnight pathology foundation model". README: "Fully open-source and improved replication of Kaiko.AI's pathology foundation model Midnight." Calling it "an open replication of Midnight" is right; "open replication and improvement" is closer still.
- "How we did it": "We start with a pretrained 1.1B parameter ViT-G DINOv2 [19]." README, Methods/Training Recipe: "Student/teacher are ViT‑G/14 with 4 register tokens." Registers appear only in the README, not the report.
- **Data:** "We train OpenMidnight only on one publicly available dataset: TCGA, [22] which contains 12K Formalin-Fixed, Paraffin-Embedded (FFPE) slides from 32 different cancer types." Also: "We only train on images from tissue slides stained with H&E". README: "Multi‑crop: 2 global crops (224) and multiple local crops (98)." Neither the report nor the README gives a µm/px scale.
- **Parameter count:** I built `dinov2_vitg14_reg` from the cached hub code without weights and loaded the cached `SophontAI/OpenMidnight/teacher_checkpoint_load.pt`, which is the file `encoder.py` loads. It has 40 blocks, a positional embedding of 257 tokens (16×16 + 1, i.e. 224-px training) and **1,134,777,344** parameters. "1.13 billion" is therefore correct for what we load; the sources say "1.1B". If Methods quotes 1.13, it is our count; "about 1.1 billion" matches all three sources.
- The report asks to be cited as "Sophont Technical Report SOPHONT-TR-2025-001, 2025. doi:10.5281/zenodo.20711012". The README asks users to cite Kaplan, Karasikov and Oquab together, which is exactly our [41–43].

**Writing ideas:**
- One sentence can carry provenance and scale: "a replication of Midnight trained on ~12,000 TCGA H&E slides (96M patches seen)". The report's own figures are "96M image patches seen, and 29M unique image patches."
- Name the checkpoint the way the README does ("teacher checkpoint"), plus the Hugging Face repository, in Code availability.

## [44] Korsunsky et al. 2019, Harmony

**Access:** NCBI E-utilities efetch JATS of PMC6884693, an author manuscript with full Online Methods (Europe PMC fullTextXML returned HTTP 500). Also the Rd documentation of R harmony 2.0.5 in the analysis environment.

**Our sentences:** Intro "harmonizes sections[44]". Methods: "R harmony 2.0.5[44] … corrected for section with Harmony (θ = 6, λ = 1, 30 iterations)".

**Verdict: SUPPORTED.** The paper uses θ and λ with the meanings we imply.
- Online Methods §2.2: "For each batch variable, we add a new parameter θ. θ decides the degree of penalty for dependence between batch membership and cluster assignment."
- §3.1: "λ is the ridge penalty hyperparameter. Larger values of λ will shrink Wk more towards 0."
- §1.1: "The Harmony algorithm inputs a PCA embedding (Z) of cells, along with their batch assignments (ϕ), and returns a batch corrected embedding (Z^)."

**Details to fix:**
- "30 iterations" should be a maximum. R Rd: "max_iter: Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step." `early_stop` is TRUE by default and our call leaves it on.
- θ = 6 is three times the default. R Rd: "Default theta=2". It is worth saying so, as `two_pass.py` does.
- In harmony 2.0.5 the signature default is `lambda = NULL`, which "will start lambda estimation mode to determine lambdas automatically". We pass λ = 1 explicitly, so write "λ fixed at 1".

**Writing ideas:**
- Name the two penalties by role in one clause: "θ (diversity penalty on clusters dominated by one section) and λ (ridge penalty on the section correction)".
- Harmony states what it corrects in one line ("inputs a PCA embedding … returns a batch corrected embedding"). That makes clear the curve is fit in corrected PCs, not genes.

## [45] Gorban & Zinovyev 2005, elastic principal graphs (Computing)

**Access:** arXiv cond-mat/0405648v2 PDF. Its p. 1 footnote reads "Computing 75, 359–379 (2005), DOI: 10.1007/s00607-005-0122-6". I checked the equations on a rendered page (p. 4).

**Our sentences:** "An elastic principal curve[45,46] was fitted"; the energy "U = (1/N) Σ‖xᵢ − n(xᵢ)‖² + λ Σ_edges … + μ Σ_inner nodes ‖nⱼ − (nⱼ₋₁ + nⱼ₊₁)/2‖² … the terms measure fit, stretching and bending"; "grown by bisecting the edge whose bisection gave the lowest energy".

**Verdict: PARTLY.**
- §2, eq. 1: "U = U(Y) + U(E) + U(R)". Also: "U(E) is the analogue of energy of elastic stretching and U(R) is the analogue of energy of elastic bending of the net."
- The fit and stretching terms match ours. The bending term, eq. 4, is the second difference ‖R(1) + R(2) − 2R(0)‖² (the arXiv text prints R(i)(0) twice, a typo). That is 4× our μ term, so μ values are not comparable across the two papers.
- §4 Adaptive strategies: the curve grows by "Determining the edge which has the largest load, by summing the number of data points (or the sum of their weights) projected to both ends of every edge". This is not our lowest-energy rule, which comes from ElPiGraph.
- It also proposes λ = λ0·s^((2−d)/d) (for a curve, λ ∝ number of nodes). ElPiGraph and our port do not scale λ.

**Fix:** keep [45,46] on the sentence that introduces the curve and its energy. Attach the growth rule ("bisecting the edge whose bisection gave the lowest energy") to [46] alone.

**Writing ideas:**
- Their physical gloss is compact: nodes are tied "by elastic bonds to the closest data points and simultaneously to the adjacent nodes". It suits a reader who does not know principal graphs.
- They state what each penalty prevents: "U(R) is a smoother term, preventing the grid from folding and twisting".

## [46] Albergante et al. 2020, ElPiGraph (Entropy)

**Access:** Europe PMC fullTextXML PMC7516753 (full Materials and Methods). I compared the energy with our port, `trajectory/elastic_curve.py`, and found it unscaled: `MSE + lam·Σ edges + mu·Σ stars`.

**Our sentences:** as for [45], plus "initialized as three nodes on the first principal component … until k = 50 (λ = 0.01, μ = 0.1)" and "The curve code is a port of elpigraph-python 0.3.4 reduced to the principal-curve case."

**Verdict: SUPPORTED.**
- §2.4, eq. 7: "We define the optimization functional for fitting a graph to data as a sum of the truncated approximation error and the elastic energy of the graph". With unit weights and no trimming radius this is our 1/N fit term.
- §2.3: "a sum of squared edge lengths (weighted by the elasticity moduli λi) and the sum of squared deviations from harmonicity for each star (weighted by the μj)". A 2-star's harmonic deviation is exactly our ‖nⱼ − (nⱼ₋₁ + nⱼ₊₁)/2‖².
- §2.5: "with a graph grammar containing only the 'bisect an edge' operation, and a graph composed by two nodes connected by a single edge as initial condition, produces an elastic principal curve". Algorithm 2, step 4: "select an injection map with the minimum energy".
- §2.9: "Empirically, it is shown that μ should be made at least one order of magnitude larger than λ. Typical starting values for λ are 10−2". Our λ = 0.01 and μ = 0.1 follow this exactly.

**Notes:**
- §2.7 gives the default start as "two nodes oriented along the first principal component". Our three-node start is elpigraph-python's `computeElasticPrincipalCurve` default (`InitNodes=3`, per the port's docstring). Say "as in elpigraph-python" rather than citing the paper for three nodes.
- §2.8 notes that leaf nodes "might not capture the extreme points of the data distribution" and offers a leaf-extrapolation function. Our "straight continuation of the end edge, up to twice that edge's length" is our own rule (`curve.py`); keep it described as ours.

**Writing ideas:**
- Use their parameter glossary verbatim in spirit: "λ—coefficient of 'stretching' elasticity; μ—coefficient of 'bending' elasticity". For curves: "For linear branches this means deviation from linearity".
- Their justification of the bending term fits in one clause: "a simple surrogate for minimizing the local curvature".
- Their framing of the fit as alternating steps is a useful one-liner: "alternate between 1) computing the partitioning K … and 2) computing new map φ provided the partitioning K".

## Hastie & Stuetzle 1989, Principal curves (JASA) (Methods addition)

**Access:** scanned JSTOR copy on the second author's University of Washington page (`dl_C/hastie_uw.pdf`). It has no text layer, so I rendered pp. 502 (abstract, Introduction) and 504 (§2.1–2.2 definitions) and read them as images. The rest of the paper was not read. The author's Stanford PDF failed on an SSL certificate error.

**Our sentence (proposed site):** "An elastic principal curve … Each unit was projected onto its nearest edge, and its position is the path length from the start of the curve to its projection".

**Verdict: SUPPORTED for the concept and the projection index.**
- Abstract: "Principal curves are smooth one-dimensional curves that pass through the middle of a p-dimensional data set, providing a nonlinear summary of the data."
- §2.2: "The projection index λf(x) of x is the value of λ for which f(λ) is closest to x."
- §2.1: "There is a natural parameterization for curves in terms of the arc length."

Our elastic curve is not self-consistent in their sense (Definition 1), so cite them for the idea and the projection, not for the fitting criterion.

**Writing ideas:**
- Position = arc-length projection index can be one sentence borrowing their terms: a unit's position is "the value of λ for which f(λ) is closest" to it, with λ in arc length.
- Abstract: "The algorithm for constructing principal curves starts with some prior summary, such as the usual principal-component line". This is a precedent for our PC1 initialization.

## [47] Van den Berge et al. 2020, tradeSeq

**Access:** Europe PMC fullTextXML PMC7058077 (Results and Methods). I compared it with our gene-model code, `pseudospace/nb_gam.py` and `stats_gam.py`. Ours is a patsy `bs(..., degree=3)` with 6 df and internal knots at quantiles, an NB GLM with log-library offset, per-gene ML dispersion "without shrinkage … no Cox-Reid adjustment" held fixed across nested fits, IRLS, and deviance LR. There is **no penalty**.

**Our sentences:** Intro: "Genes are then modelled along t with negative-binomial generalized additive models[47]". Methods: "negative-binomial generalized linear models with a log(library size) offset … a curve (cubic B-spline in x with 6 degrees of freedom, knots at reference x quantiles, the same basis in every section) … [47]".

**Verdict: PARTLY.**
- **Supported:** Methods: "Differences in sequencing depth or capture efficiency between cells are accounted for by cell-specific offsets Ni." Also: "Our default computational implementation sets K = 6." And: "The knots are by default positioned according to the quantiles of the pseudotime values." Basis functions are "enforced to be the same for all genes and lineages".
- **Differs:** "the coefficients βglk are shrunken by substracting a penalty λgβgTSβg from the log-likelihood function … selected using generalized cross-validation". Dispersion comes from "the joint estimation of the NB regression parameters in μgi and dispersion parameter ϕg". Tests are Wald statistics. Ours has no penalty, a fixed ML dispersion, and LR tests.

**Fix:**
- Methods: "…as in tradeSeq[47] (six cubic basis functions with knots at quantiles and a log-library offset) but without a smoothing penalty; nested models were compared by likelihood ratio at fixed dispersion[McCarthy 2012]."
- Intro: change "negative-binomial generalized additive models[47]" to "negative-binomial spline regression[47]" (or keep "GAM" loosely, but then Methods should say "an unpenalized GAM").

**Writing ideas:**
- tradeSeq's one-line model statement is worth copying in form: "Ygi~NB(μgi,ϕg), log(μgi)=ηgi" plus offset. A displayed equation would make our four nested designs unambiguous.
- They call knots "relative markers of progress along the trajectory". This is a good caution for our x-quantile knots.

## NEW: McCarthy, Chen & Smyth 2012 (edgeR GLM)

**Access:** Europe PMC fullTextXML PMC3378882.

**Our sentence:** "The dispersion was estimated per gene and section … held fixed for the nested fits, which were fitted by iteratively reweighted least squares", and "The likelihood ratio of clusters + curve against clusters".

**Verdict: SUPPORTED** for NB log-linear GLMs with a library offset, fitted at a given dispersion and compared by LR between successive models.
- Results: "The total library size Ni serves as an offset in the linear model predictor, capturing the dependence of counts on sequencing depth."
- Methods, GLMs: "Assuming that an estimate is available for ϕg, so the variance can be evaluated for any value of μgi, GLM theory can be used to fit a log-linear model".
- Table 1 caption: "Differentially expressed genes are detected by likelihood ratio tests between successive models".

Caveat: they estimate dispersion by Cox–Reid adjusted profile likelihood; we use plain ML. Do not cite them for the dispersion estimator.

**Writing ideas:**
- Their nested-model table ("Patient", "Patient + tissue", "Patient × tissue" with interpretations) is a good template for our bulk / clusters / curve / clusters + curve, each with a one-line interpretation.

## U-Net, Ronneberger et al. 2015 (placeholder "U-Net-style decoder[ref]")

**Access:** arXiv 1505.04597v1 PDF.

**Verdict: SUPPORTED.**
- Abstract: "The architecture consists of a contracting path to capture context and a symmetric expanding path that enables precise localization."
- §2 Network Architecture: "Every step in the expansive path consists of an upsampling of the feature map … a concatenation with the correspondingly cropped feature map from the contracting path".

Our decoder fuses a frozen ViT token pyramid with a CNN stem pyramid (`decoder.py`: "U-Net/FPN fusion"), so "U-Net-style" is fair.

**Writing ideas:**
- U-Net is also the precedent for two other steps in our Methods. Tiled whole-slide inference: "seamless segmentation of arbitrarily large images by an overlap-tile strategy". Our edge-facing-another-object class: "a weighted loss, where the separating background labels between touching cells obtain a large weight".

## Kirillov et al. 2019, Panoptic segmentation (placeholder "[ref Kirillov 2019]")

**Access:** CVF open-access PDF (CVF pagination 9404–9413; IEEE/Crossref 9396–9405, use the DOI's).

**Our sentences:** "The checkpoint with the highest validation panoptic quality (PQ)[ref Kirillov 2019] was kept." Also: "matching predictions one-to-one to annotations at intersection over union (IoU) above 0.5. In mouse, PQ was 0.880 (segmentation quality 0.914, detection quality 0.963…)".

**Verdict: PARTLY.**
- §4.1: "a predicted segment and a ground truth segment can match only if their intersection over union (IoU) is strictly greater than 0.5."
- §4.2: "PQ can be seen as the multiplication of a segmentation quality (SQ) term and a recognition quality (RQ) term". Also: "RQ is the familiar F1 score". So PQ = SQ × RQ (= DQ) is supported, but **Kirillov's name is "recognition quality"**. "Detection quality (DQ)" is HoVer-Net's term (Graham 2019: "the detection quality (DQ) is the F1 Score"). Our `instance_metrics.py` notes that both names are in use.
- §4.2: "We calculate PQ for each class independently and average over classes." We select on the class-agnostic overall PQ (`evaluate_per_class` → "overall"; AGENTS.md "pooled class-agnostic PQ"). Say so.

**Fix:** "…panoptic quality (PQ = segmentation quality × recognition quality, with one-to-one matching at IoU > 0.5[Kirillov]), computed class-agnostically". Or keep "detection quality" and add HoVer-Net.

**Writing ideas:**
- Kirillov's void rule is the precedent for scoring only annotated objects under partial annotation. Quote: "all pixels in a predicted segment that are labeled as void in the ground truth are removed from the prediction and do not affect IoU computation". One clause would justify our evaluation on annotated objects only.
- "SQ is simply the average IoU of matched segments" is a ready gloss for the SQ number.

## Watershed: Meyer & Beucher 1990 / Beucher & Meyer 1993 (placeholder "marker-controlled watershed[ref]")

**Access:**
- Meyer & Beucher 1990: **INACCESSIBLE — abstract only.** ScienceDirect returned 403, and Semantic Scholar marks the open-access PDF "CLOSED". The earlier notes read only the abstract, via BibBase ("introduce the notion of markers of the objects to be segmented").
- Beucher & Meyer 1993 chapter: full PDF from the authors' lab site (people.cmm.minesparis.psl.eu/users/beucher/publi/SB_watershed.pdf). Its page numbers are 433–481, which match the Crossref record of the CRC reissue.

**Our sentence:** "Objects were separated by marker-controlled watershed[ref]: seeds of high interior probability were flooded up to the edge probabilities."

**Verdict (1993 chapter): SUPPORTED.** The implementation is `skimage.segmentation.watershed(..., markers=...)` in `postprocess/decode.py`.
- §"B. Flooding from Selected Sources": "suppose we know before flooding which minima correspond to the centers of the nuclei and which to the background." Then: "we will bore a hole only in these minima before immersing the relief."
- Same section: "A dam is constructed only for separating floods originating from different pierced minima." Also: "it is not even necessary to choose sources within the minima. Any region may be chosen." This supports seeds that are thresholded interior regions.

**Fix:** cite the 1993 chapter, which was read in full. Keep 1990 only if the authors prefer the journal paper and accept that it is verified from the abstract only. Citing scikit-image for the implementation is optional (not checked here).

**Writing ideas:**
- Their two-step logic gives a one-line description: "markers fix one source per object; the relief (edge probability) is flooded only from these sources, and dams where floods meet are the object boundaries."

## HCOP, Yates et al. 2021 (placeholder "one-to-one mouse orthologs (HCOP[ref])")

**Access:** Europe PMC fullTextXML PMC8574622. I checked our rule in `pseudospace/cross_species.py`: HCOP pairs with `min_support=3`, a support-weighted reciprocal-best match, and only symbols with exactly one partner on both sides.

**Verdict: PARTLY.** HCOP is the right citation for the source. The paper does not, however, define a one-to-one set; the user derives it from support counts.
- The HCOP pipeline: "It is these combined ortholog predictions, essentially a pair of genes with a list of associated databases that support that orthology assertion, that are displayed as HCOP results."
- HCOP Results: "The higher the number of resources agreeing on a particular assertion the more confidence we can have that this assertion reflects a genuine ortholog."
- Applications: VGNC's "high confidence" rule is "a one-to-one orthology relationship between protein coding genes and must be predicted by at least three out of the following four resources". This is a published precedent for a support threshold of three.

**Fix:** "Human genes were mapped to mouse orthologs from HCOP[ref], keeping pairs supported by at least three of its resources and, within multi-member orthogroups, the reciprocal-best pair by support, so that every gene has one partner." Confirm that the package path used for the paper applies the same `min_support=3` rule.

**Writing ideas:**
- HCOP's wording for confidence is plain and reusable: "the number of resources that support each prediction".

## Scanpy, Wolf et al. 2018 (placeholder "Scanpy 1.11.5[ref]")

**Access:** Europe PMC fullTextXML PMC5802054. I also read the installed Scanpy 1.11.5 docstring and signature of `highly_variable_genes`.

**Verdict: SUPPORTED** as the toolkit citation.
- Results: "SCANPY provides preprocessing comparable to SEURAT [10] and CELL RANGER [6], … clustering similar to PHENOGRAPH". The paper itself uses Louvain; Leiden is cited separately, as we do.
- Our HVG thresholds (mean 0.0125–3, normalized dispersion ≥ 0.5) are exactly the defaults of 1.11.5 (`min_mean=0.0125, max_mean=3, min_disp=0.5, flavor='seurat'`). The docstring says "Annotate highly variable genes :cite:p:`Satija2015,Zheng2017,Stuart2019`". If Methods names the method, "Scanpy's seurat-flavour dispersion criterion (default thresholds)" plus Satija 2015 (not read) would be accurate.

**Writing ideas:**
- Saying "default thresholds" is shorter than listing them, and it signals that nothing was tuned.

## scFates, Faure et al. 2023 (placeholder "scFates 1.2.5 (principal tree)[ref Faure 2023]")

**Access:** Europe PMC fullTextXML PMC9805561 (Application Note, full body). I checked the run settings in `analysis/paper/fig3_baselines.py` and `results/paper_figures/fig3/baselines/methods_notes.md`: `tl.explore_sigma` then `tl.tree(method='ppt', Nodes=50, ppt_lambda=100, ppt_nsteps=50)`. The ppt fits are flagged "not converged" in both cohorts.

**Verdict: PARTLY.**
- Body: "The user has the choice of running either ElPiGraph or SimplePPT algorithms to learn a principal graph." Also: "the initial R version included a tree inference approach inspired by SimplePPT (Mao et al., 2015)". Since our run used SimplePPT, cite Mao 2015 with Faure 2023 (Mao 2015 verified only through Crossref: SIAM SDM 2015, pp. 792–800).
- Body: "scFates additionally implements ElPiGraph … allowing investigators to impose topological constraints on trajectories to fit single curves". This is useful context: scFates can fit the same kind of curve we fit, but our baseline was deliberately a tree.

**Fix:** "scFates 1.2.5 (SimplePPT principal tree, 50 nodes)[Faure 2023; Mao 2015]". The non-convergence flag belongs in the Supplementary settings.

**Writing ideas:**
- scFates' description of pseudotime assignment fits a baseline sentence: cells are placed "by projecting cells according to their position relative to the principal points".

## [22] Leiden, Traag et al. 2019 (Methods use)

**Access:** Europe PMC fullTextXML PMC6435756.

**Our sentence:** "Units were clustered by Leiden[22] (resolution 2, 30 nearest neighbours, 2,000 highly variable genes)".

**Verdict: SUPPORTED.**
- Abstract: "We prove that the Leiden algorithm yields communities that are guaranteed to be connected."
- Leiden Algorithm section: "The Leiden algorithm consists of three phases: (1) local moving of nodes, (2) refinement of the partition and (3) aggregation of the network based on the refined partition".
- Introduction, on modularity with resolution γ: "Higher resolutions lead to more communities, while lower resolutions lead to fewer communities."

Implementation detail (not a citation issue): `sc.tl.leiden(flavor="igraph", n_iterations=2)` optimizes modularity with a resolution parameter. "Resolution 2" therefore means γ = 2 in modularity, not CPM.

**Writing ideas:**
- The resolution gloss can be cut to "(modularity resolution γ = 2; higher γ gives more clusters)".

## Benjamini & Hochberg 1995 (Methods addition: "Benjamini–Hochberg corrected within each contrast")

**Access:** text-layer PDF of the JRSS-B article (dcscience.net copy; pp. 289–300 visible on p. 1). I did not confirm that the copy is authorized.

**Verdict: SUPPORTED.**
- Summary: "It calls for controlling the expected proportion of falsely rejected hypotheses—the false discovery rate."
- §3, the procedure (math symbols re-typed, because the scan's text layer garbles them): "let k be the largest i for which P(i) ≤ (i/m) q*; then reject all H(i) i = 1, 2, . . ., k."
- Theorem 1 (verbatim): "For independent test statistics and for any configuration of false null hypotheses, the above procedure controls the FDR at q*."

**Note for authors:** our gene-set statistics share genes, so they are not independent. BH is still standard here, but a reviewer may ask. Our empirical p values come from 1,000 random sets; state the estimator (e.g. (b + 1)/(n + 1)), because the smallest attainable p bounds what BH can call. Neither point was checked against a paper here.

**Writing ideas:**
- "q ≤ 0.05 (Benjamini–Hochberg, within each contrast)" is enough. Use "FDR" rather than "corrected p".

## HGNC gene groups (placeholder "(HGNC groups 752 and 417)[ref HGNC]")

**Access:** Seal et al. 2023 (Europe PMC fullTextXML PMC9825485) and Gray et al. 2016 (Europe PMC fullTextXML PMC4739092). Group membership came from `genenames.org/cgi-bin/genegroup/download?id=752|417&type=branch` (today's files have the same row counts as our 2026-10-12 copies: 474 and 53 genes). I also read `data/external/hgnc/SOURCE.txt` and notebook 101.

**Verdict: SUPPORTED** (Seal 2023, the current HGNC resource paper).
- Seal 2023, HGNC DATA: "The HGNC data box at the top of the Report presents HGNC-curated nomenclature, HGNC ID, chromosomal location and membership of gene groups." Also: "Related gene groups are curated into hierarchies that support browsing between groups."
- Gray 2016: "Gene families may be arranged in complex hierarchies that feature many different levels and each level has its own family report."

**Facts for the sentence:**
- 752 = "Solute carriers" (474 genes in 76 subgroups such as "Solute carrier family 22"). 417 = "ATP binding cassette transporters" (53 genes in 7 subgroups).
- Notebook 101 pairs genes **within subgroups** (`FAMILY = … ['Group name']`) and keeps a gene's first group when it is in several.

**Fix:** "Within each HGNC transporter family (the subgroups of the 'Solute carriers' and 'ATP binding cassette transporters' groups, IDs 752 and 417; downloaded 12 October 2026)[Seal 2023]…".

**Writing ideas:**
- HGNC releases data nightly ("We release our data nightly"), so a download date is the version.

## Enrichr (placeholder "(Enrichr libraries)[refs]")

**Access:** Kuleshov 2016 (Europe PMC fullTextXML PMC4987924). Xie 2021 (NCBI efetch PMC8152575; Europe PMC fullTextXML returned 500). The Enrichr front page's citation request and `datasetStatistics` were read for the three library names. Chen 2013 was not read.

**Verdict: SUPPORTED** as the source of the library files.
- Kuleshov 2016, abstract: "Enrichr currently contains a large collection of diverse gene set libraries available for analysis and download."
- Xie 2021, Basic Protocol 1: "Enrichr enables users to submit lists of human or mouse genes to compare against numerous gene set libraries of known biological function".
- Enrichr site: "Please acknowledge Enrichr in your publications by citing the following references", followed by Chen 2013, Kuleshov 2016 and Xie 2021.
- Library metadata: `MSigDB_Hallmark_2020` has 50 terms; `Reactome_2022` 1,818 terms (link reactome.org/download-data); `KEGG_2019_Mouse` 303 terms.

**Fix:** "Gene sets were the Enrichr[Xie 2021; Kuleshov 2016] libraries MSigDB_Hallmark_2020[Liberzon], Reactome_2022[Gillespie] and KEGG_2019_Mouse[Kanehisa]". Note that the human-symbol libraries were mapped to mouse through the ortholog map (`coordinate_inputs.py` passes `ortholog_map`).

## MSigDB Hallmark, Liberzon et al. 2015

**Access:** NCBI efetch JATS PMC4707969 (Results and Methods; Europe PMC fullTextXML returned 500).

**Verdict: SUPPORTED.**
- Abstract: "Each hallmark in this collection consists of a "refined" gene set, derived from multiple "founder" sets, that conveys a specific biological state or process and displays coherent expression."
- Results, "Generating the Hallmark Collection": it describes the 50 themes and "Founders for the final set of 50 hallmarks (Table 1)".

**Writing ideas:**
- "50 hallmark sets (refined, low-redundancy)" explains in four words why Hallmark sits beside the much larger Reactome and KEGG.

## Reactome, Gillespie et al. 2022

**Access:** Europe PMC fullTextXML PMC8689983.

**Verdict: SUPPORTED** as the citation matching Enrichr's `Reactome_2022`. Enrichr gives no Reactome release number; this paper describes release 78 (October 2021), the nearest NAR description.
- Introduction: "The Reactome Knowledgebase is distinctive in focusing its manual annotation effort on a single species, Homo sapiens".
- Introduction: "Reactome (version 78, October 2021) has entries for 10 726 (52.5%) of the 20 442 predicted human protein-coding genes".

**Note:** because the curation is human-only, Methods should say the sets were carried to mouse genes through the one-to-one ortholog map.

## KEGG, Kanehisa et al. 2019 and Kanehisa & Goto 2000

**Access:**
- Kanehisa 2019: Europe PMC fullTextXML PMC6324070.
- Kanehisa & Goto 2000: **INACCESSIBLE — abstract only.** NCBI efetch had no body. OUP's article-pdf link and Europe PMC's `?pdf=render` both returned HTML bot pages.
- The KEGG website (kegg.jp/kegg/kegg1.html, read through WebFetch's summarizer) asks users to cite four papers: Kanehisa et al. NAR 2027 (gkag1042), Kanehisa et al. NAR 2025 (53, D672–D677), Kanehisa Protein Sci 2019 (28, 1947–1951) and Kanehisa & Goto 2000.

**Verdict (2019): SUPPORTED.**
- Abstract: "KEGG pathway maps, BRITE hierarchies and KEGG modules have been developed as generic molecular networks with KEGG Orthology nodes of functional orthologs".
- Introduction: "pathways for a specific organism can be automatically reconstructed from the gene set in its genome once genes are annotated with KO identifiers." This is why a mouse KEGG library exists.

**Fix:** cite Kanehisa & Goto 2000 (canonical) + Kanehisa et al. 2019 NAR (matching the 2019 library). If the journal wants KEGG's own current list, add the NAR 2025 paper (not read).

## CZ CELLxGENE Discover / Census (placeholder "single-cell PT segment profiles of six donors in CZ CELLxGENE Census[ref]")

**Access:**
- Europe PMC fullTextXML PMC11701654 (NAR 2025).
- Census documentation home page, "Citing Census".
- CELLxGENE curation API (public) for dataset 09b518f9-da64-44cc-aec8-70a89d55611f → collection 4cbb929b… "Multimodal benchmarking dataset for renal cortex characterization (mBDRC)", DOI 10.1101/2025.03.06.637075. Crossref title: "Systematic evaluation of single-cell multimodal data integration for comprehensive human reference atlas" (Acera-Mateos et al., bioRxiv, 6 Mar 2025, 27 authors; no journal version linked).
- Local provenance: `data/external/census_pt_segments/SOURCE.txt` ("CELLxGENE Census 2025-11-08").

**Verdict: SUPPORTED.**
- NAR 2025, Introduction: "Census is an Application Programming Interface (API) and data object hosted by CZ CELLxGENE that enables efficient access and custom slices thereof for programmatic access and computational use cases."
- Census docs: "To cite the project please follow the citation guidelines offered by CZ CELLxGENE Discover. To cite individual studies please refer to the tutorial Generating citations for Census slices." So the source study must be cited as well.

**"Six donors" check:** the human pseudobulk has 10 donors. Six of them have ≥ 50 cells in both "kidney proximal convoluted tubule epithelial cell" and "epithelial cell of proximal tubule segment 3" (21_015, HCA_29_ref, HCA_51_ref, HCA_55_ref, PRECISE-34, PRECISE-36), which matches the code's `min_cells=50`.

**Fix:** "…against PT segment pseudobulks (proximal convoluted tubule versus S3) of the six donors with at least 50 cells per segment in the mBDRC renal-cortex dataset[Acera-Mateos 2025], obtained through CZ CELLxGENE Census (release 2025-11-08)[CZI 2025]." Note that the human comparison is PCT (S1+S2) versus S3, not S1 versus S3.

**Writing ideas:**
- Census encourages naming the release and the slice. "Census release 2025-11-08, dataset 09b518f9…" in the Data availability statement makes the profile reproducible.
