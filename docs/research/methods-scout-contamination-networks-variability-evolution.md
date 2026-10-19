# Methods scout: contamination, active modules, differential variability, expression vs sequence evolution (2026-10-18)

A literature review run for the discovery-methods scoreboard (`docs/results/pt-discovery-methods.md`). Four topics
were read from full texts, mostly Europe PMC XML, with each paper's read status tagged. Statements marked *our
inference* are not from any paper. The source texts were kept as local scratch downloads and are not in the
repository.

**How this was read.**
- Four parallel readers covered the four topics, mostly from Europe PMC full-text XML.
- I re-checked these numbers against the downloaded texts:
  - Oliveira's Fig. 2e table (from the figure image);
  - Janssen's 3–35%;
  - smiDE's 178 → 56 → 19;
  - Mitchel's crypt/villus p values;
  - DOMINO's Table 1;
  - Nguyen's 137/150;
  - the DiffVar type I table;
  - memento ρ 0.71;
  - the QRscore AUPRC table;
  - Xuanyuan's 15.3% → 33.2%;
  - Yakubovsky's R = 0.15;
  - Sepp's 23%;
  - the BioMart output;
  - the Ensembl 99 Compara header.
- **Access tags:**
  - [full text] = main text and methods read;
  - [+ supplement] = supplement also read;
  - [main text only];
  - [read via WebFetch] = read through an extraction tool, as quotes rather than line by line;
  - [abstract only] = no claims about methods or results are taken from it.
- Statements marked *our inference* are not from any paper.
- Not repeated: tradeSeq, Lamian, condiments, Droin, Halpern, Ben-Moshe, Dubuis, Tkačik, Squair, Zimmerman, and Crow 2022 (already used in notebook 97).

---

## Short answers

**A. Non-native signal per unit.**
- **Amounts by platform:**
  - Visium HD FFPE colon (Oliveira 2025, Fig. 2e): 0.97–1.73% of marker UMIs land in the adjacent mask. Median displacement of mislocalised UMIs is 3.05–6.01 µm. Density in the adjacent mask is about 0.8–3.6% of the source's.
  - Visium HD outside the tissue (Ren 2025): extra-tissue counts are about 11% of intra-tissue counts, with a tail of several hundred µm.
  - Visium v1: 4.2–14.2% of UMIs are off tissue (10x mouse kidney 7.7%). The chimeric lower bound is 10–15%; the authors estimate 20–30% in total.
  - Droplet mouse kidney (Janssen 2023): 3–35% background per cell, dominated by PT transcripts.
- **Precedent that neighbour composition creates false biology: yes.**
  - Mitchel (Nat Genet 2026): false DE along a positional axis. Goblet cells look "stem/TA-like" in crypts and "enterocyte-like" on villi.
  - Neighbour and region contrasts (smiDE, SPLIT, Proseg).
  - Marker decay that looks like a gradient (SpotClean).
  - Blurred zonation (ResolVI, liver).
  - Ambient-driven batch effects (SoupX).
  - SoupX and DecontX applied to Visium create false spatial patterns (SpotClean supplement).
- **Not found:** any kidney Visium HD contamination estimate.

**B. Fewest non-specific modules.** The answer depends on the null.
- Under gene-score shuffling (EHR), DOMINO is best: 0.84–0.94, against 0.05–0.72 for others.
- Against degree-preserving rewired networks (Lazareva 2021), DOMINO is nearly the only method that beats them.
- Under planted truth (Yang 2024), DOMINO is poor; FDRnet, BioNet, NetMix2 and HotNets keep non-target genes low.
- No benchmark compares module finders with ORA or GSEA under one null.
- Under sample relabeling (Nguyen 2019), GSEA with sample permutation has 0 biased pathways and Fisher ORA has 137 of 150.
- Gene-shuffle nulls are blind to inter-gene correlation; our specimen relabeling is not.

**C. Differential variability (DV).**
- No gene-level DV test is valid at 2 vs 2; published methods need about 8–50 samples per group.
- The workable recipe:
  - a per-specimen NB2 location–scale fit along x;
  - residuals against the mean–variance trend;
  - an explicit unit-size term;
  - specimen-level moderated F and camera;
  - calibration by decoys and a parametric bootstrap.
- No kidney injury study has tested tubule-to-tubule variance.

**D. What predicts conserved expression patterns.**
- **Predictors:**
  - Protein constraint is weak but consistent: r ≈ 0.15–0.38. In Strand's quartiles, regional log-ratio r is about 0.53 in the low-dN/dS quartile against 0.33 in the high one.
  - Duplication is the strongest categorical predictor.
  - Essentiality, disease status and LOEUF also predict conservation.
  - Pleiotropy does not.
  - Expression level confounds all of the above.
- **Gap:** no study relates conservation of a zonation or axial pattern to dN/dS, gene age or paralog status.
- **dN/dS source:**
  - Not in current Ensembl. Release 116 is the final legacy release, and its Compara README marks dN "currently unused".
  - The last release with values is Ensembl 99 (Compara FTP).
  - The Ensembl 80 BioMart archive still serves `mmusculus_homolog_dn` and `_ds`; tested, and it returns values.

---

## A. Contamination in segmentation-aggregated spatial data

**Search note.** I found no Nat Methods 2024 spot-swapping paper. The peer-reviewed paper is Ni et al., Nat Commun 2022 (bioRxiv 2021). A later method, SpaDiff (Cell Rep Methods 2026, PMID 42468536), was skimmed. 10x technical note CG000618 was not found and not read.

### A1. SpotClean / spot swapping

**Citation.** Ni Z, Prasad A, Chen S, …, Newton MA, Kendziorski C. Nat Commun 2022;13:2971. PMID 35624112, doi:10.1038/s41467-022-30587-y. [full text + supplementary information + peer-review file]

**Method.**
- True counts split as Y_gt = S_gt (stays in the spot) + B_gt (bleeds out). Observed X_gj = S_gj + Σ_t B_g,t,j.
- S_gt ~ Pois(μ_gt(1−r_β)).
- B_g,t,j ~ Pois(μ_gt·r_β·[(1−r_γ)·w_tj + r_γ/K]), where:
  - r_β is the global bleed rate;
  - w_tj is a normalised Gaussian kernel exp(−d²/2σ²) in physical distance;
  - r_γ/K is a slide-wide "distal" constant.
- Fitting:
  - r_β, r_γ and σ by least squares on per-spot totals (σ by grid search, the rates by L-BFGS-B);
  - μ by EM;
  - per-spot contamination ĉ_t = expected received counts / expected X_·t.
- Needs at least 25% of spots off tissue.
- The peer-review reply notes that a purely local kernel could not reproduce the off-tissue counts.

**Benchmark.**
- 4.2–14.2% of UMIs fall in background spots across 11 public Visium sets (Supp. Table 1); mouse kidney 7.7%.
- Human skin / mouse gut chimeras:
  - LPSS (the lower bound on cross-species reads in tissue spots) 14.9%, 10.2% and 13.3%;
  - 17.6–26.9% of UMIs in background;
  - about 20–30% swapped in total, per the review reply.
- GFAP and ERBB2 decay with distance from their source at the same rate into adjacent tissue as into empty background (Supp. Figs 8–10). Swapping therefore creates **apparent gradients across boundaries**.
- MSE reduced by more than 20% in most simulations (mouse_kidney SimI 12.1 → 7.8).
- More than 13% of breast "malignant" spots disappear after correction.
- SoupX and DecontX on Visium "impose artificial patterns upon non-SV genes" and inflate SpatialDE FDR (Supp. Figs 13c, 14).

**Limits.** Visium v1 (55 µm spots). Needs off-tissue spots. Isotropic kernel with one global r_β per slide, so no differences by region or condition.

**Adaptation (A-ring).**
- Inputs, per specimen:
  - 2 µm bins outside all segmented units;
  - their distance to the nearest PT unit edge and to the nearest non-PT unit;
  - PT-exclusive genes (Lrp2, Slc34a1, Slc5a2 and human orthologs).
- Model: fit λ(d) = a·exp(−d²/2σ²) + c with an NB GLM, separately for interstitial and off-tissue bins.
- Statistics: σ (bleed length), the near-boundary bleed fraction, and the floor ratio c/in-unit density.
- Decision: bleed is condition-dependent if both IR mice exceed both controls by more than 1.5× in σ or in c/density.
- **Against us:** condition-dependent bleed, which would put a contamination term into whole-PT IR-vs-control contrasts.

**Code.** R/Bioconductor `SpotClean`. SpaDiff is Python.

### A2. Diffusion measured directly in Visium HD

**Citations.**
- (a) Oliveira MF, Romero JP, Chung M, et al. Nat Genet 2025;57:1512–1523. PMID 40473992. [full text; the Fig. 2e table read from the image]
- (b) Ren P, Zhang R, et al. Nat Commun 2025;16:9232. PMID 41107232. [full text; Fig. 2e,f read from the image, so approximate; supplement not read]
- Context:
  - bin2cell, Polański et al., Bioinformatics 2024, PMID 39250728 [full text];
  - ENACT, Kamel et al., Bioinformatics 2025, PMID 40053700 [full text].

**Method.**
- (a) For goblet-cell markers (CLCA1, FCGBP, MUC2), muscularis markers and lymphoid markers, UMIs were counted in hand-drawn source masks against the adjacent mask. Each mislocalised UMI's distance to the source edge was measured. FFPE, CytAssist, Space Ranger 3.0.
- (b) 8 µm bins outside the tissue, by distance to the tissue edge. "Diffusion level" = mean extra-tissue / mean intra-tissue normalised counts.

**Benchmark.**
- (a) Fig. 2e, three samples:
  - source masks hold 98.28 / 98.32 / 99.04% of marker UMIs and adjacent masks 1.73 / 1.68 / 0.97%;
  - median mislocalised distance 3.83 / 3.05 / 6.01 µm;
  - density source vs adjacent 0.19 vs 0.006, 1.10 vs 0.04, and 0.58 vs 0.0045 transcripts/µm². The adjacent density is about 0.8–3.6% of the source's.
- (b) Visium HD FFPE ratio about 0.11 (about 0.22, 0.06 and 0.03 across three tumour types), against about 0.17 for Stereo-seq v1.3. A low tail extends to several hundred µm.
- bin2cell corrects row/column "striping" of 2 µm bins by quantile scaling.
- ENACT: about 25% of bins overlap more than one cell (Xenium-derived synthetic data); area-weighted assignment did best (F1 0.787).
- **Neither bin2cell nor ENACT quantifies diffusion.**

**Limits.** Vendor authors, colon only, 3–4 regions. "Adjacent" is one neighbouring mask, so long tails are missed. Ren measures tissue-to-glass leakage. Nothing on kidney or fresh-frozen tissue.

**Adaptation.** Use these numbers as priors for the planted-contamination test: about 1–4% foreign density within 4–6 µm, plus a slide-wide floor. A-ring checks the same quantities in our specimens.

**Code.** github.com/10XGenomics/HumanColonCancer_VisiumHD. bin2cell and ENACT are Python.

### A3. Segmentation errors confound region and neighbour DE

**Citation.** Mitchel J, Gao T, Petukhov V, Cole E, Kharchenko PV. "Impact and correction of segmentation errors in spatial transcriptomics", Nat Genet 2026;58:434–444, PMID 41559218. [bioRxiv v1 full text, doi 10.1101/2025.01.02.631135; the published version was not read]

**Method.**
- Admixture is shown three ways:
  - scRNA markers appear in other cell types;
  - those markers sit at the cell periphery;
  - stronger membrane stain lies between foreign and native molecules.
- Correction:
  - weighted NMF on per-molecule neighbourhood-composition vectors, within each cell type;
  - admixture factors are labelled by hand;
  - a CRF assigns each molecule to a factor;
  - admixture molecules are removed.

**Benchmark.** CosMx NSCLC, MERFISH ileum, MERFISH hypothalamus, Xenium pancreas.
- Tumour vs stroma fibroblasts: 11 of the top 36 up-genes are malignant markers (GSEA p = 2.8×10⁻⁴).
- **Goblet cells, crypt vs villus:** crypt goblet cells "up-regulate" stem and TA markers (GSEA p = 0.016); villus goblet cells "up-regulate" enterocyte markers (p = 0.003). This is the closest analogue to our positional axis.
- Microglia in the anterior commissure gain oligodendrocyte markers (p = 0.0005).
- Macrophages with a fibroblast neighbour: the DE list is dominated by collagens (p = 0.00029).
- A neighbour classifier beats chance (p = 10⁻⁵), but 4 of its top 7 genes are malignant markers.
- 25% of DIALOGUE's top 32 genes and 30% of scITD's top 23 are admixture.
- After correction, the foreign markers lose significance and native markers stay (Supp. Fig. 2e).

**Limits.** Imaging platforms only. Manual factor labels. Some residual admixture remains.

**Adaptation (A-adj, a negative control).**
- Within each of 8 position bins and each specimen, split PT units by whether a non-PT unit (TAL, DCT, CNT, CD or glomerulus) lies within 10 µm edge to edge.
- Compare the two groups by NB pseudobulk at specimen × bin level, and with an L2 logistic classifier using spatial-block cross-validation.
- Score: the number of foreign-type markers in the top 20 genes, and fgsea on foreign marker sets.
- Null: 1,000 permutations of the adjacency label within specimen × bin.
- Then test whether adjacency probability varies with x or between IR and control.
- **Against us:** a foreign signature, adjacency that varies with x, and overlap with our positional or injury calls.

**Code.** R (pagoda2, scITD); the NMF/CRF repo URL was not located.

### A4. smiDE: overlap-ratio filter and neighbour-expression covariate

**Citation.** Vasconcelos AG, McGuire D, Simon N, Danaher P, Shojaie A. Genome Biol 2026;27:21. PMID 41519776. [full text; supplement not read]

**Method.**
- **ORM_g = mean expression of g in neighbouring other-type cells within r / mean in the analysed type.** Genes with ORM > 1 are excluded.
- Neighbour expression is also added as a covariate: the distance-weighted sum of g in other-type neighbours, rank-inverse-normal transformed, as a log-linear term.
- The rationale is a DAG (Fig. 3b): segmentation error is an unobserved confounder correlated with the contrast.
- Spatial random effects (GP, SPDE, BYM2) can be added.

**Benchmark.**
- CosMx NSCLC: likely false positives (genes "Not Detected in Immune cells" in HPA) were removed for 112/138 (macrophage), 68/72 (CD4), 27/27 (CD8) and 14/16 (Treg).
- MERFISH neurons against microglia exposure, starting from 178 naive hits:
  - 56 (31.5%) removed by the ORM filter;
  - 19 (10.7%) removed by the covariate;
  - 100 (56.2%) removed by spatial modelling;
  - 3 survived.

**Limits.** Cell-level data. The ORM cut-off of 1 is heuristic. The covariate assumes contamination scales with neighbour expression. The random-effect part does not transfer to 2–4 specimens.

**Adaptation (A-cov): the top recommendation in this topic.**
- Per specimen and position bin:
  - ORM_g(x) = mean in non-PT units within 20 µm / mean in PT units;
  - E_gu = Σ_j w_uj·y_gj over non-PT neighbours, with w from the A-ring kernel, rank-inverse-normal within specimen.
- Refit the NB GAM with s(x) + E_gu, and recompute the specimen-level shape LR, moderated F and cameraPR z.
- Flag a gene if ORM > 1 in the bins carrying its effect, or if its shape LR drops by more than 50% when E is added.
- A call is contamination-robust if it keeps ≥ 70% of its effect with flagged genes dropped. This is the existing R4 threshold.
- **Against us:**
  - more than 20% of primary calls, or of the 74 reversed genes, are flagged;
  - or flagged genes are over-represented among calls (Fisher OR > 2).

**Code.** R `smiDE` (Bruker Spatial Biology). The filter and covariate are about 30 lines of Python.

### A5. Spillover-aware correction and segmentation: SPLIT, Proseg, ResolVI

**Citations.**
- SPLIT: Bilous M et al. Nat Methods 2026;23:1152–1162, PMID 42062553. [full text; supplement not read]
- Proseg: Jones DC et al. Nat Methods 2025;22:1331–1342, PMID 40404994. [full text (Nature HTML); supplement not read]
- ResolVI: Ergen C, Yosef N. Nat Methods 2026;23:2003–2014, PMID 42786290. [full text; supplement not read]

**Methods.**
- **SPLIT:**
  - RCTD doublet weights w1 and w2.
  - **Spillover index** = cosine similarity between w2 and the secondary type's share of the cell's neighbourhood (up to 20 cells within 15 µm), with a permutation null.
  - Purification: x = w1·ref1/(w1·ref1 + w2·ref2) ⊙ x_obs.
  - A logistic near-vs-far test counts contaminant markers among the top 20 genes.
- **Proseg:**
  - Cellular-Potts probabilistic segmentation with transcript repositioning.
  - Repositioning prior κ·N(s′, σ_a) + (1−κ)·N(s′, σ_b), with κ = 0.2, σ_a = 4 and σ_b = 0.5 fixed; units are presumably µm.
  - "Spurious" gene pairs are those whose conditional co-expression rises at least 1.5× when nuclei are expanded.
- **ResolVI:**
  - A VAE in which observed counts = own expression (NB) + neighbour-derived expression (NB) + background (Poisson).
  - Neighbour weights have a Dirichlet prior with an RBF distance kernel over 20 neighbours.
  - Has been applied to Visium HD 8 µm bins.

**Benchmark.**
- **SPLIT:** T cells near malignant cells show spurious malignant and epithelial signatures (Fig. 4g). Even nucleus-only segmentation is contaminated. ovrlpy and ResolVI empty many cells, and ResolVI sometimes suppresses real T-cell programs.
- **Proseg:**
  - Image-based defaults leave at least 20% of transcripts unassigned (Cellpose about 26% in CosMx).
  - Tumour-adjacent vs non-adjacent macrophages give tumour genes in DESeq2; Proseg reduces this but does not remove it (Fig. 5d).
  - Small cells are "smothered": median 64 transcripts for T cells vs 180 for macrophages.
- **ResolVI:**
  - Stereo-seq: 14.8% of Slc17a7⁺ cells are Mbp⁺, falling to < 0.01% after correction.
  - Liver: cholangiocyte markers appear pericentrally before correction, a zonation pattern blurred by background (Fig. 4b).

**Limits.** Mostly imaging platforms. SPLIT needs a complete reference. Proseg's prior is not fitted. ResolVI outputs model counts that do not suit our NB specimen-level tests.

**Adaptation (A-erode, after Proseg's expansion test).**
- Rebuild every unit from masks eroded by 1–2 bins (2–4 µm) and dilated by 1–2 bins.
- R_g = (counts per bin in the dilation ring / counts per bin in the eroded core), scaled by the library-wide median. Genes with R_g ≥ 1.5 are boundary-enriched.
- Rerun the positional statistics on the eroded units. A call is robust if it keeps ≥ 70% of its effect. Also test whether R_g varies with x.
- **Against us:** boundary-enriched genes are over-represented among calls, or R_g rises toward S3 or in IR.
- Needs no reference and no model, which suits two species with different panels.

**Code.**
- SPLIT: R/Python, github.com/bdsc-tds.
- Proseg: Rust, github.com/dcjones/proseg.
- ResolVI: Python, scvi-tools.

### A6. Ambient RNA: SoupX, DecontX, CellBender, and a kidney ground truth

**Citations.**
- SoupX: Young MD, Behjati S. GigaScience 2020;9:giaa151, PMID 33367645. [full text]
- DecontX: Yang S…Campbell JD. Genome Biol 2020;21:57, PMID 32138770. [full text]
- CellBender: Fleming SJ et al. Nat Methods 2023;20:1323–1335, PMID 37550580. [main text and methods (HTML)]
- Janssen P et al. Genome Biol 2023;24:140, PMID 37337297. [full text]

**Methods.**
- **SoupX:**
  - b_g = the empty-droplet profile.
  - n_gc = m_gc + N_c·ρ_c·b_g.
  - With genes known to be absent from cell c: ρ_c = Σn_gc / (N_c·Σb_g). There is also an automated mode-of-markers version.
- **DecontX:** LDA-like. Each cell is a θ-weighted mix of its population profile and a contamination profile made from all other populations. Fitted by variational inference.
- **CellBender:**
  - Noise ~ Pois[(1−ρ)·ε·d_drop·χ_ambient + ρ·ε·(y·d_cell + d_drop)·χ̄].
  - True counts are NB with a learned latent prior.
  - The output is tuned by a nominal FPR.
- **Janssen:** pools *M. m. castaneus* and *domesticus* kidney cells, and gets a per-cell background fraction ρ_c from foreign-allele counts with a binomial likelihood.

**Benchmark.**
- SoupX: about 1–2% contamination in species-mixing data; PBMC about 5–6%.
- DecontX: median 1.09% (human cells) and 2.75% (mouse cells), range 0.43–45%, down to 0.25% after correction.
- CellBender: off-target counts per cell, median 225 → 19. At TPR 0.809, its FPR is 0.062 against DecontX's 0.142.
- **Mouse kidney (Janssen):**
  - background is 3–35% of UMIs per cell; frozen-tissue snRNA is highest (35% vs 11%);
  - the PT marker Slc34a1 is detected in about 60% of non-PT cells, falling to 54% (SoupX), 7% (CellBender), 9% (DecontX with an empty-droplet profile) or 2% (DecontX default);
  - marker log2FC shrinks in noisier replicates (Fig. 4).
- SoupX: ambient profiles are specific to each channel and create batch effects (Fig. 4E).

**Limits.** All need empty droplets or known-negative genes. On spatial data SoupX and DecontX create false patterns (A1).

**Adaptation (A-amb, a diagnostic only).**
- b_g from off-tissue bins per specimen.
- ρ_u per PT unit by SoupX's formula over PT-negative genes: Umod, Slc12a1, Aqp2, Nphs2, Ptprc, Hbb, and their human orthologs.
- Test ρ_u ~ s(x) within each specimen, and IR vs control by relabeling.
- If ρ_u varies more than 2-fold, add log ρ_u as a covariate. **Never feed decontaminated counts into the gene models.**

**Code.** SoupX and DecontX (celda) are R; CellBender is Python.

### A synthesis

**How non-native signal is estimated.** Five approaches:
1. Off-tissue or empty profiles.
2. Species or strain mixing (ground truth).
3. Genes known to be absent from a cell type.
4. Neighbour-weighted mixtures: Gaussian kernel plus floor (SpotClean); RBF Dirichlet (ResolVI); w2-vs-neighbourhood cosine (SPLIT); neighbour covariate and ORM (smiDE).
5. Perturbing the segmentation: nuclei vs expanded masks (Proseg).

**Our inference.**
- A within-specimen shape test cancels constant contamination. It does not cancel contamination that changes with position.
- PT neighbours plausibly change along the axis (glomeruli and DCT near S1; TAL, DTL and CD near S3), by condition (IR infiltrates) and by species.
- So A-cov and A-erode matter more than any decontamination.

---

## B. Active modules on interaction networks versus gene-set enrichment

### B1. DOMINO and its false-term benchmark (EMP / EHR)

**Citation.** Levi H, Elkon R, Shamir R. Mol Syst Biol 2021;17:e9593. PMID 33471440, doi:10.15252/msb.20209593. [full text + supplement (Appendix, Expanded View, review file)]

**Method.**
- Input: a yes/no active-gene set (e.g. DE q < 0.05) and a network.
- **Step 0:** Louvain (resolution 0.15) splits the network into "slices" of more than 3 nodes.
- **Step 1:** a slice is relevant if it is hypergeometric-enriched for active genes (BH q < 0.3). It is tested only if it holds ≥ 10% of active genes, or its active fraction is ≥ α = min(0.7, (active/all)(1 + 100/N)).
- **Step 2a:** a prize-collecting Steiner tree (FAST-PCST) reduces each relevant slice to one connected piece.
  - Prizes come from linear-threshold spreading: an inactive node turns on when its active neighbours' summed 1/k influence exceeds 0.5. It earns β^l, where l is the round it turned on and β = max(0, 1 − 3·active/N).
  - Edge cost is 0 if the edge touches an active node, otherwise 1 − 10⁻⁴.
- **Step 2b:** pieces larger than 10 nodes are split by Newman–Girvan, stopping when log(size)/log(N) ≤ modularity.
- **Step 3:** keep modules at hypergeometric Bonferroni q < 0.05.

**How false terms are measured (EMP).**
- Gene activity scores are shuffled uniformly, about 5,000 times, ignoring the network. The method and the GO tests are rerun on each shuffle.
- The null for each term is its best −log10 p over modules; e(t) = 1 − F_t(s).
- A term is "validated" if hypergeometric q ≤ 0.05 and e(t) ≤ 0.05.
- **EHR** = validated terms / all hypergeometric-significant terms. mEHR is the same ratio per module.
- Also scored:
  - richness: non-redundant terms after REVIGO;
  - intra-module homogeneity;
  - robustness: AUPR and F1 of validated terms under 60–90% subsampling, 100 times each.
- GO setup: Biological Process terms of 5–500 genes, with evidence codes IPI, IEP, IGI, IMP, HMP, HGI and HEP excluded. Background = network genes.

**Benchmark.**
- **Table 1, 10 expression datasets on DIP:**

  | Method | EHR | mEHR | F1 | AUPR at 80% | Richness | Homogeneity |
  |---|---|---|---|---|---|---|
  | DOMINO | 0.891 ± 0.129 | 0.868 | 0.486 | 0.776 | 30.5 | 2.376 |
  | NetBox | 0.719 ± 0.425 | | | | | 2.575 (best) |
  | HotNet2 | 0.424 | | | | | |
  | BioNet | 0.398 | | | | | |
  | jAM_SA | 0.236 | | | | | |
  | KPM | 0.149 | | | | | |
  | jAM_greedy | 0.052 | | | | | |

- **GWAS:** DOMINO 0.844, NetBox 0.694, HotNet2 0.081.
- **Larger networks:** DOMINO 0.939 ± 0.046 on STRING > 900 (11,972 nodes), against NetBox 0.18. On HuRI DOMINO is about 0.88–0.94. It beat NetBox on 23 of 24 criteria.
- **Topology-driven terms:** 37% of Biological Process terms (2,450/6,573) are enriched in modules found from topology alone. They were not rejected more often by EMP, so filtering them cannot replace EMP.
- **Cost:** EMP took days on a 44-core server.

**Limits.**
- The null destroys inter-gene correlation. Plain ORA on the same list would score EHR ≈ 1 by construction.
- Binary input only. Thresholds set empirically. Run by DOMINO's own authors. No ORA/GSEA comparator and no planted truth.

**Adaptation (B-D).**
- **Gene statistics:** from notebook 71, not refit. AKI uses its 10,198-gene universe; human vs mouse uses the 1:1 universe of 7,971 genes.
- **Network:** STRING v12.0 (mouse `10090.protein.links.v12.0`, 75 MB), score ≥ 900, restricted to the universe. Background = universe ∩ network.
- **Active set: a fixed top 300 by moderated-F p** (sensitivity 150), the same size in the real contrast and in every relabeling. Notebook 93's ORA looked specific only because the relabelings selected no genes.
- **Sets:** GO Biological Process 5–500 with DOMINO's evidence-code exclusion, plus notebook 71's sets.
- **Nulls:**
  1. EMP with 1,000 shuffles;
  2. 100 degree-preserving rewirings (double-edge swaps, 10×|E|);
  3. the two balanced specimen relabelings, with DOMINO rerun on each relabeling's own top 300.
- **Comparators:**
  - cameraPR;
  - ORA on the same top 300;
  - a degree-only baseline (the 300 highest-degree genes).
- **Keep DOMINO only if all hold:**
  - the swap rule passes (NCDR < 0.25, and each relabeling calls < 5% of tested terms);
  - EHR ≥ 0.8;
  - its validated terms exceed the 95th percentile of the rewired networks;
  - at least one specimen-replicated program is not called by cameraPR or the F-test.
- **Against:**
  - rewired networks give as many terms;
  - or the AKI modules are notebook 93's hub module (Pcna, Ccnd1, Cdkn1a, Myc, Cdk2, Cdk4) and it also appears under relabeling;
  - or its terms are a subset of cameraPR's.

**Code.**
- Python (MIT), github.com/Shamir-Lab/DOMINO: last push 2023-10, Python 3.8, `pip install domino-python` 0.1.1, also on bioconda.
- EMP: github.com/Shamir-Lab/EMP.
- Web server: Levi et al. Bioinformatics 2022, PMID 35139202 [full text]. It converts mouse IDs to human orthologs on its preloaded networks.

### B2. Lazareva 2021: most active-module methods learn degree, not edges

**Citation.** Lazareva O, Baumbach J, List M, Blumenthal DB. Brief Bioinform 2021;22(5):bbab066. PMID 33782690. [read via WebFetch; figures and supplement not seen, so only numbers in quoted sentences are used]

**Method.**
- 10 methods: ClustEx2, COSINE, DIAMOnD, DOMINO, GiGA, GXNA, KeyPathwayMiner, GrandForest, plus the permutation-based Hierarchical HotNet and NetCore.
- 5 networks: BioGRID, APID, STRING, HPRD, IID.
- Random networks, 10 per generator:
  - REWIRED (exact degree kept);
  - EXPECTED_DEGREE;
  - SHUFFLED (gene IDs permuted, so the shape is kept but each gene's degree is not);
  - Barabási–Albert;
  - Erdős–Rényi.
- 5 GEO diseases; seed genes at Mann–Whitney P < 0.001/n.
- Scores: mutual information with the phenotype; −log10 GSEA P of the disease's KEGG pathways; DisGeNET overlap.
- Real vs random compared by one-sided Mann–Whitney.

**Results.**
- Most methods beat SHUFFLED, scale-free and uniform networks, but **not** the degree-preserving generators.
- KEGG and DisGeNET scores rise with a module's mean degree.
- Against the degree-preserving generators, only DOMINO reached significance on KEGG; on DisGeNET, DOMINO and DIAMOnD did.
- Hierarchical HotNet and NetCore failed on all three scores.

**Limits.** The reference sets (KEGG, DisGeNET) are themselves degree-biased. No gene-correlation null. Five diseases.

**Adaptation.**
- The rewired arm of B-D.
- Hub share: the fraction of module genes in the top 5% by degree, real vs rewired.
- **Against:** no difference between real and rewired networks.

**Code.** github.com/dbblumenthal/amim-test-suite.

### B3. Yang 2024: a planted benchmark where DOMINO does poorly

**Citation.** Yang L, Chen R, Goodison S, Sun Y. Brief Bioinform 2024;26(1):bbae692. PMID 39737568. [full text; supplement not read]

**Method.**
- 12 methods on BioGRID, iRefIndex, ReactomeFI and STRING.
- Planted truth: target Reactome or CORUM sets of at least 10 genes. Gene p values are Beta for targets and Uniform otherwise, at signal 0.01–0.11, over 4,400 datasets.
- List-input methods use local fdr < 0.1.
- Metrics: gene F; subnetwork F; FDR = the share of non-target genes in reported subnetworks.

**Results.**
- Best gene F: FDRnet, MuST, BioNet, ROBUST, NetMix2, DIAMOnD, NetCore. Only FDRnet was also good at the subnetwork level.
- FDR controlled by FDRnet, BioNet, NetMix2, Hierarchical HotNet and HotNet2. Hierarchical HotNet was near 0 at signal 0.01–0.08.
- Not controlled: DIAMOnD, DOMINO, RegMOD and ClustEx (the last two near 1).
- "DOMINO did not perform well in our analysis", which the authors blame on its fixed up-front partition.
- Real-cancer estimated FDR: DIAMOnD 0.78, RegMOD 0.91, ClustEx 0.94.
- Methods that correct for degree over-penalise hubs.

**Limits.** Independent gene scores. Mutation-like data.

**Adaptation.** Run DOMINO and FDRnet inside notebook 100's planting machinery (`pseudospace/planting.py`, relabeling swap_IR2A2_Ctrl1A4, 15 sets per replicate with 60% of members planted), with two kinds of plant:
- (a) notebook 71's sets;
- (b) connected STRING subgraphs of 20–60 genes, which is the network's best case.

Score: planted sets recovered at matched relabeled false calls.
- **Keep** if the network method recovers at least 2 more planted sets per replicate than cameraPR in (b), and loses at most 1 in (a).
- **Against:** no gain even in (b).

### B4. Choobdar 2019: the DREAM disease-module challenge

**Citation.** Choobdar S, Ahsen ME, Crawford J, et al. Nat Methods 2019;16:843–852. PMID 31471613. [full text and supplementary figure legends; supplementary tables not read]

**Method.**
- Modules found from topology only on 6 anonymised networks: STRING, InWeb, OmniPath, co-expression, cancer dependency, homology.
- 42 single-network and 33 multi-network methods; modules of 3–100 genes.
- Scored by Pascal: SNP χ² within ±50 kb, LD-corrected, competitive against network genes, BH 5%.
- A method's score = the number of modules associated with at least one of 104 held-out GWAS.

**Results.**
- Top five methods scored 55–60. K1 (diffusion-state distance plus spectral clustering) was the most robust over 1,000 subsamples.
- Only 46% of trait modules were found by more than one method, and 17% across networks.
- Modularity correlated with score at only r = 0.45. Multi-network methods were no better.
- An all-genes background inflated the number of trait-associated random modules (Supp. Fig. 4).
- 33% of STRING consensus trait modules had generic functions.

**Lessons for us.** Use the background = universe ∩ network, and do not merge modules across methods into a consensus.

**Cheap extra.** Partition STRING-mouse once and score the partitions with cameraPR as extra sets under the swap rule. **Against:** relabeled calls grow in proportion to the modules added.

### B5. Network propagation

**Cowen L, Ideker T, Raphael BJ, Sharan R. Nat Rev Genet 2017;18:551–562.** PMID 28607512. [abstract only; paywalled] No claims are taken from it.

**Picart-Armada S, et al. PLoS Comput Biol 2019;15:e1007276.** PMID 31479437. [full text; S1 not read]
- **Method:** 12 rankers (personalised PageRank, diffusion raw/gm/mc/z, knn/wsld, SVMs, random forest, COSNet) against these baselines:
  - EGAD;
  - **pr**: PageRank on topology alone, which ignores the input;
  - randomraw;
  - random.
  22 OpenTargets diseases on STRING and OmniPath, with cross-validation that handles protein complexes.
- **Results:**
  - Keeping complexes in one fold dropped the random forest from about 12 to fewer than 4.5 top-20 hits.
  - Top-20 hits: random 0.1, EGAD about 2, best diffusion about 3.75, random forest about 4.4.
  - pr beat random on every metric.
  - With genetic input: top-20 pr 0.43–0.47 vs best 0.65–0.70; top-100 pr 2.25–2.38 vs best 3.29–3.45. **Propagation adds little over degree.**

**Picart-Armada S, et al. Bioinformatics 2021;37:215–222.** PMID 33070187. [read via WebFetch]
- f = K·y, normalised by a permutation z score or a Monte Carlo p.
- Node bias has a mean term Σ_j K_ij and a variance term, both tied to degree.
- Normalisation helps when true positives are not aligned with degree and hurts when they are.

**Adaptation (B-P, secondary).**
- Random walk with restart (0.5) of z² on STRING-mouse ≥ 900, normalised by 1,000 label permutations within degree bins.
- Use it only as a gene ranking inside notebook 100's planted benchmark. Do not run set tests on propagated scores, because propagation makes neighbouring genes correlated.
- **Against:** pr matches propagation.

**Code.**
- R: diffuStats (b2slab/diffuStats), b2slab/genedise.
- Python: networkx `double_edge_swap`.
- Hierarchical HotNet: raphael-group/hierarchical-hotnet.

### B6. Network and topology methods vs ORA/GSEA on specificity

**Nguyen T-M, Shafi A, Nguyen T, Draghici S. Genome Biol 2019;20:203.** PMID 31597578. [full text; supplement not read]
- **Method:**
  - 13 methods, 8 without pathway topology and 5 with.
  - Positives: 75 human datasets, plus 11 mouse knockouts scored by AUC.
  - **Null:** 2,000 healthy 15-vs-15 relabelings.
  - Bias: skewness of each pathway's null p values (γ₁ > 0.1).
  - Fixed 400 DE genes per dataset.
- **Pathways biased toward small P:**

  | Method | Biased pathways |
  |---|---|
  | Fisher ORA | 137/150 |
  | Wilcoxon rank-sum | 114/150 |
  | CePaGSA | 112/186 |
  | GSEA | 0 |
  | PathNet | 0 (but biased toward large P) |

- Every method except GSEA had at least 66 biased pathways.
- Topology-based methods ranked targets better (p = 8.8×10⁻³). Knockout AUC: ROntoTools 0.799, GSEA 0.763.
- Small pathways were the most biased.

**Buzzao D, Castresana-Aguirre M, Guala D, Sonnhammer ELL. Brief Bioinform 2024;25(2):bbae069.** PMID 38436561. [full text; supplement not read]
- **Method:** 14 methods: overlap (EASE, Fisher), per-gene (CAMERA, fGSEA, GSEA, GSA, GSVA, PADOG, ROAST), topology (CePa-ORA, SPIA) and network crosstalk (ANUBIX, BinoX, NEAT, netPEA). 82 datasets. The negative benchmark resamples gene labels 30 times.
- **Results:**
  - Best G-mean: BinoX and NEAT 0.62, fGSEA and GSEA 0.61.
  - Median FPR < 3%: EASE, Fisher, CAMERA, SPIA.
  - Median FPR > 16%: ROAST, GSVA, BinoX, NEAT.
  - Close to 5%: GSA, fGSEA, GSEA, PADOG, ANUBIX.
  - Crosstalk FPR tracks hub share: ρ 0.66 (BinoX) and 0.71 (NEAT), against 0.06 for ANUBIX.
- **Adaptation:** ANUBIX as one extra row in notebook 71's table (top 300, STRING ≥ 800, swap rule). **Against:** NCDR above cameraPR's, or false calls that track hub density.
- **Code:** R, github.com/MiguelCastresana/anubix; NEAT is the R package `neat`.

**Also read, not central.** Batra 2017, NPJ Syst Biol Appl, PMID 28649433 [full text]: no de novo tool was consistently best, and BioNet and GXNA fail on variance-shift signals.

### B synthesis

- By EHR (gene-score shuffle), DOMINO gives the fewest non-specific terms (0.84–0.94).
- By non-target genes under planted truth, FDRnet, BioNet, NetMix2 and the HotNets do.
- Under sample relabeling (our analogue), GSEA had 0 biased pathways and ORA was worst.
- EHR and gene-label nulls cannot see correlation- or specimen-driven non-specificity. **Every network arm has to be rerun inside our relabelings, with a fixed top-K, and with degree controlled** by rewiring and a degree-only baseline.
- Expected outcome (inference): DOMINO passes EHR, but on AKI it adds little beyond the hub-driven cell-cycle and Myc calls cameraPR already makes.

---

## C. Differential variability along the axis

**The question.** At a fixed PT position, does injury make tubule-to-tubule expression more variable? The test has to separate real heterogeneity from five other sources:
- Poisson noise;
- mean–variance dependence;
- unit and library size;
- coordinate error and the dropped injured S3 units;
- 2-vs-2 replication.

### C1. BASiCS residual overdispersion

**Citations.**
- Eling N, Richard AC, Richardson S, Marioni JC, Vallejos CA. Cell Syst 2018;7:284–294, PMID 30172840. [full text; supplement not read]
- Vallejos CA, Richardson S, Marioni JC. Genome Biol 2016;17:70, PMID 27083558. [full text; additional file 1 not read]

**Method.**
- X_ij ~ Pois(φ_j·ν_j·μ_i·ρ_ij), where ν_j ~ Gamma is technical noise (θ) and ρ_ij ~ Gamma(1/δ_i, 1/δ_i). Marginally this is NB with biological over-dispersion δ_i.
- Eling: log δ_i = f(μ_i) + ε_i, with ε ~ t₅ and f = linear term + 10 Gaussian RBFs in log μ. **ε_i is the residual over-dispersion.**
- The log2 change in δ splits into a mean part f_A(μ_A) − f_B(μ_B) and a residual part ε_A − ε_B.
- A gene is DV if P(|Δε| > ψ0) > α_R, with ψ0 = 0.41 (a 50% change) and α_R set by EFDR 10%. Genes need ≥ 2 detecting cells per group.
- Without spike-ins, technical noise comes from batches, with μ and δ shared across them.

**Benchmark.**
- Vallejos 2016 (74 single cells vs 76 pool-and-split, 9,378 genes):
  - about 2,000 genes more over-dispersed in single cells and about 60 in pool-and-split (EFDR 5%);
  - the mean test called 4,710 genes at τ0 = 0 and 559 at τ0 = 0.4;
  - the authors say FDR is poorly calibrated at τ0 = 0, so a minimum effect threshold is part of the method.
- Eling 2018:
  - the regression prior stabilised δ under subsampling of 939 neurons to 50–500 cells;
  - δ tracked smFISH CV² and ε tracked residual CV²;
  - mixing 5 activated into 93 naive T cells created "DV" genes such as Il2 (Fig. S5): **outliers or contamination create DV**.

**Limits.** Discrete groups. Slow MCMC (more than 11 h per dataset, per Gatlin 2025). Cells treated as exchangeable. With δ shared across batches there is no between-animal dispersion, so 2 vs 2 would be anti-conservative (*inference from Eq. 37–38*).

**Adaptation.**
- In each specimen × bin, the ML log α̂_g minus a robust across-gene trend f_sb(log μ̄) gives ε_gsb.
- Report the trend shift f_IR − f_Ctrl separately as global DV; it is confounded with section quality.
- Use ψ0 = 0.41 as the minimum gene-level effect.

**Code.** R/Bioconductor BASiCS (`Regression=TRUE, WithSpikes=FALSE`); MarioniLab/RegressionBASiCS2017.

### C2. Martinez-Jimenez 2017 (application)

**Citation.** Martinez-Jimenez CP et al. Science 2017;355:1433–1436. PMID 28360329. [main text only (author manuscript); supplementary methods and the 2019 erratum not read]

**Method and results.**
- BASiCS on 1,514 CD4 T cells: young vs old mice, two species, naive vs 3 h stimulation.
- Activation changed 2,063 genes. Among genes with unchanged mean, variability fell (p < 10⁻¹⁰).
- In old mice the core activation program was more variable (p < 0.05) and **expressed in fewer cells**, in both species.
- The variance increase reflects a smaller responding fraction, not noisier responders.

**Adaptation.**
- For injury-program genes, compute the fraction of units with detectable expression per specimen × bin, alongside dispersion.
- Two-state (spared vs injured tubules) predicts a lower responding fraction and higher dispersion. A graded shift predicts neither.
- The AKI universe lacks *Havcr1, Lcn2, Mki67, Top2a*, so use the INJ/TPT/CYC genes that are present.

**Code.** MarioniLab/ImmuneAging2017 (R).

### C3. DiffVar

**Citation.** Phipson B, Oshlack A. Genome Biol 2014;15:465. PMID 25245051. [full text; additional file 1 not read]

**Method.**
- A Levene-type test: each sample's absolute (or squared) deviation from its group mean, scaled by n_k/(n_k+1).
- Any limma design is fitted to the deviations, with empirical-Bayes moderated t and BH.
- Because it is limma, CAMERA, ROAST and TREAT apply directly.

**Benchmark.**
- Type I error on resampled normal kidney samples:

  | Method | Nominal 0.001 | Nominal 0.01 |
  |---|---|---|
  | Bartlett | 0.009 | 0.0257 |
  | F test | 0.0089 | 0.0253 |
  | DiffVar (absolute) | 0.001 | 0.0104 |

- With outliers, Bartlett 0.0138 vs DiffVar 0.00097.
- "Almost perfect FDR control", but no power at n = 5; about 10 per group needed.
- Overlap between the top differentially methylated and top DV CpGs was below chance.

**Limits.** Gaussian scale. Samples are the replicates, so it is pseudoreplication on units.

**Adaptation.**
- The specimen-level engine.
- Response: ε_gsb (C1), or the mean |Pearson residual| from our NB mean model.
- limma with `~ specimen + bin + condition:bin` and weights; moderated F on condition:bin; cameraPR on signed DV z.

**Code.** R `missMethyl::varFit`; easy to reimplement on our moderated-F code.

### C4. MDSeq, Le Priol 2023 and GAMLSS

**Citations.**
- MDSeq: Ran D, Daye ZJ. NAR 2017;45:e127, PMID 28535263. [full text; supplementary methods not read]
- Le Priol C, Azencott C-A, Gidrol X. PLoS Comput Biol 2023;19:e1010342, PMID 36893104. [full text; supporting figures not read]
- de Jong 2019, Physiol Genomics, PMID 30875273. [paper not accessible (403); model read from the authors' code ExpVarQuant.R]

**Methods.**
- **MDSeq:** E[Y] = μ and Var = φμ, with log μ and log φ each linear in their own design. Optional zero inflation, Wald and threshold tests, outlier removal, TMM, Benjamini–Yekutieli.
- **GAMLSS:** `gamlss(y ~ 0 + x + offset(log lib·nf), sigma.fo = ~ 0 + x, family = NBI())`, with an LR test against `sigma.fo = ~ 1`.

**Benchmark (Le Priol: NB simulation, 10,000 genes, 50% DD genes with FC ≥ 1.5, outliers in 10%).**
- GAMLSS and DiffDist were best and needed 20–30 samples per condition to reach AUC 0.8; DiPhiSeq needed 40 and MDSeq about 50.
- A shift in mean "sharply reduces" MDSeq's performance, because under Var = φμ a mean change alone moves φ; NB2 (Var = μ + αμ²) does not have this problem.
- GAMLSS and DiffDist got the sign wrong in 5.7% of calls at n = 50 and 13–14% at n = 20, always for dispersion decreases.

**Limits.** All need about 20–50 independent samples per group.

**Adaptation.**
- **Use NB2.** Per specimen: log μ_i = log L_i + B(x_i)'β and log α_i = U(x_i, m_i)'γ.
- This generalises `nb_gam._ml_alpha` from one α per gene to a spline in x. Significance comes from the specimen-level step (C3).

**Code.** R: MDSeq (zjdaye/MDSeq), gamlss, ExpVarQuant (Vityay/ExpVarQuant).

### C5. memento

**Citation.** Kim MC, …, Ye CJ. Cell 2024;187:6393–6410. PMID 39454576. [full text; supplementary figures not read]

**Method.**
- Hypergeometric capture model with method-of-moments estimators: σ̂² = mean[(Y² − Y(1−q))/(N²q²)] − μ̂².
- **Residual variance** = what is left after a degree-2 polynomial of variance on mean, per group.
- Bootstrap of cells.
- Replicates: estimate per replicate × condition, then regress log residual variance on condition across those points.

**Benchmark.**
- Residual variance against smFISH: ρ 0.71, against 0.56 for the naive estimator and 0.61 for BASiCS (14 genes).
- 2 donors: 761 DV genes; IFN-β reduced the variability of 78% of ISGs.

**Limits.** No continuous covariates (the authors suggest binning). The bootstrap is within replicates, so it is optimistic at 2 replicates (*inference*).

**Adaptation.** Bin x, estimate per specimen, regress across specimens. Use a spatial-block bootstrap within specimen × bin, because Moran's I is 0.06–0.14.

**Code.** Python, yelabucsf/scrna-parameter-estimation.

### C6. Recent and contextual papers

- **QRscore** (Zhou 2025, Cell Rep Methods, PMID 40902592) [full text; supplement not read].
  - NB/ZINB score-weighted rank test.
  - AUPRC at n 113/249: GAMLSS 0.798, QRscore 0.788, DiffVar 0.774, MDSeq 0.701.
  - At 5–15 per group, MDSeq and GAMLSS inflate FDR; at least 8 per group recommended. R/Bioconductor.
- **spline-DV** (Gatlin 2025, NPJ Syst Biol Appl, PMID 40113778) [full text; supplement not read]. Ranks genes only (p < 0.05, no FDR); not recommended.
- **Jiang, Wang & Li 2026** (Cell Syst commentary, PMID 42617588) [full text]. A cell-level Levene test on 70 cells from 11 mice gave 169 DV genes where pseudobulk DE found none: a cautionary case of cells treated as replicates.
- **Bahar Halpern 2015** (Mol Cell, PMID 25728770) [full text; extended procedures not read].
  - Hepatocytes in the same zone vary up to 100-fold because of bursting.
  - Tetraploid cells have lower CV (median −13%) because they average more bursts.
  - Counts were taken in a fixed 500 µm³ volume. **This is the precedent for our unit-size term.**
- **Kidney injury: no DV test found.**
  - Melo Ferreira 2021 (JCI Insight, PMID 34003797) [full text; supplement not read]: Visium, IRI with n = 1. Patchy low-count regions; sham had about 2× the nuclei of "preserved" IRI regions and 4× those of "low" regions. **Injury co-varies with library size.**
  - Xuanyuan 2025 (Nat Commun, PMID 40813851) [full text; supplement not read]: Xenium. Injured PT rose from 15.3% to 33.2% of cells between 4 h and 12 h, spreading from the corticomedullary junction into cortex; healthy-S3 neighbourhoods were lost.

### Design-specific confounds (*our inference*)

1. **Unit size.**
   - A unit summing m sub-units, each NB with dispersion α_c, with a shared tubule multiplier of CV² τ², has α(m) ≈ τ² + α_c/m.
   - The target is τ². Smaller IR units inflate α. Use m = n_spots.
2. **Coordinate error.**
   - Var(e) = σ_x² adds about g′(x)²σ_x² to α, so steep genes gain spurious DV wherever σ_x is larger, possibly in IR.
   - Control: σ_x(x) per condition from notebook 90's decoder, plus fold-protected genes.
3. **Selection.** Deep-S3 units were dropped (IR 251 and 32 vs control 783 and 700; pt-package-issues item 7). This biases IR variance down. Run a sensitivity analysis with the placed c15/c19 units.
4. **Contamination.** A neighbour-composition covariate (links to A).
5. **2 vs 2.**
   - Only three splits exist, so a per-gene permutation p cannot go below 1/3.
   - Calibrate across genes: decoy FDP from the relabelings; a parametric bootstrap with condition-free α(x, m) and real L, m and x; cameraPR at set level.
   - The A2/A4 relabeling carries the A2/A4 shared factor, so it is a conservative null.

### Proposed frozen analysis (C\*)

1. **Mean model per specimen:** NB with offset log(library), 10-df spline in x.
2. **Dispersion model per specimen × bin (8 bins):** α = τ²_gsb + κ_gsb/m by ML, with block-bootstrap standard errors.
3. **Residualise:** ε_gsb = log τ̂² − f_sb(log μ̄), reporting the global shift separately.
4. **Specimen-level test:** limma moderated F on condition:bin, then cameraPR for INJ, TPT, CYC and notebook 71's sets.
5. **Cheap first look:** the MAD of the INJ-score residual and the detection fraction per specimen × bin.

**Decision.** Set-level heterogeneity at bin b requires all four:
- camera FDR ≤ 0.05;
- both IR specimens above both controls;
- |z| above both relabelings;
- ≥ 50% of the effect survives the 1/m, g′²σ_x² and placed-unit checks.

**Mechanistic prediction.**
- A patchy two-state injury predicts excess τ² ∝ p(1−p), largest where the IR mean effect is intermediate.
- A graded shift (the D4 one-Gaussian result) predicts Δτ² ≈ 0.

---

## D. Expression evolution versus sequence evolution

**Our data (notebook 90 D).**
- 7,971 1:1 genes, end contrast E.
- Categories: 189 conserved and 74 reversed (both graded in all four specimens); 229 mouse-only and 70 human-only.
- Replicate r is 0.89 for mice and 0.63 for the human sections. The species-only imbalance is partly human reliability.

### D1. Liao & Zhang 2006

**Citation.** Liao BY, Zhang J. Mol Biol Evol 2006;23:530–540. PMID 16280543. [read via WebFetch; supplement not read]

**Method.**
- GNF atlas: 26 shared tissues, 10,607 ortholog pairs.
- Relative abundance within species.
- Divergence = Euclidean d or 1 − Pearson r. dN and dS from Ensembl (codeml).

**Results.**
- d correlates with dN (r 0.187, P < 1e-80), dN/dS (0.152) and dS (0.119).
- **"We did not observe positive correlations between 1 − r and either dN, dS, or dN/dS"**: the shape measure showed nothing.
- 83.9% of genes diverge less than the neutral expectation.

**Relevance.** Our E is a shape measure. A null result for us would have a precedent.

### D2. Strand 2007: the closest precedent

**Citation.** Strand AD et al. PLoS Genet 2007;3:e59. PMID 17447843. [full text + Fig. 2 image; supplementary tables not read]

**Method.**
- Human and mouse arrays for three matched brain regions.
- limma regional log ratios within each species.
- 8,499 orthologs ranked by Ensembl dN/dS (dS capped at 2× the median) and by % identity.
- Cross-species r of the log ratios computed per quartile.

**Results.**
- Overall r ≈ 0.45.
- **By dN/dS quartile:** about 0.31–0.35 in the highest and 0.52–0.54 in the lowest, in all three contrasts (read from the plot).
- By % identity: from about 0.32 up to 0.55–0.62.
- Top 100 genes: 9 discordant (r 0.86). Top 500: 93 discordant (r 0.75).

**Limits.** No reliability correction per quartile, no covariates, no CIs.

**Adaptation.** Repeat Fig. 2 with E_h against E_m by dN/dS quartile, correcting reliability within each quartile.

### D3. Brawand 2011

**Citation.** Brawand D et al. Nature 2011;478:343–348. PMID 22012392. [main text + online methods; the supplementary note (OU/selection models) not read]

**Method.** 6 organs × 10 amniotes, n ≈ 2. 5,636 amniote 1:1 orthologs. Neighbour-joining trees on 1 − ρ, and a maximum-likelihood test for expression shifts.

**Results.**
- Neural tissue evolves slowest and testis fastest.
- **Kidney** matches liver and heart across amniotes (P > 0.1) and is slower than them in primates (P < 0.05).
- Mouse branches are shorter (purifying selection).
- Divergence saturates.
- About 63% of amniote orthologs have at least one expression shift.
- **No sequence-divergence analysis in the main text.**

**Data.** GEO GSE30352.

### D4. Warnefors & Kaessmann 2013

**Citation.** Warnefors M, Kaessmann H. Genome Biol Evol 2013;5:1324–1335. PMID 23781097. [full text; supplement not read]

**Method.**
- Brawand data: 3,749 amniote and 10,227 primate genes.
- Expression divergence ED = Euclidean distance of log expression over branches (AncML).
- dN from PRANK alignments and codeml.
- Partial Spearman controlling 11 factors: level, τ, neural bias, dS, connectivity, onset, age, family size, length, GC, essentiality.

**Results.**
- ρ(ED, dN) = 0.27 (amniotes) and 0.38 (primates).
- Among highly expressed genes it persists after all 11 controls (partial ρ 0.16, P < 1e-15).
- **Among lowly expressed genes, controlling for level alone abolishes it.**
- Essentiality lowers ρ from 0.28 to 0.24.

### D5. Paralogs and the ortholog conjecture

**Citations.**
- Kryuchkova-Mostacci N, Robinson-Rechavi M. PLoS Comput Biol 2016;12:e1005274, PMID 28030541. [full text; S1 not read]
- Huminiecki L, Wolfe KH. Genome Res 2004;14:1870–1879, PMID 15466287. [read via WebFetch]

**Kryuchkova-Mostacci.**
- Tissue-specificity τ correlates between orthologs at 0.74–0.89 in tetrapods, declining linearly with divergence time.
- Between paralogs it is significantly lower, even for primate-specific paralogs, and does not depend on paralog age.
- After small-scale duplication, the lower-expressed copy loses specificity.
- Code: Figshare 10.6084/m9.figshare.3493010.v2 (R).

**Huminiecki (similar profile = R > 0.8).**

| Ortholog class | Similar profile |
|---|---|
| Strict 1:1 | 13.7% (181/1,325) |
| One species-specific duplication | 7.4% (11/148) |
| Several duplications | 4.9% (5/102) |
| Independent duplications in both species | 1.9% (1/52) |

Mean Ka/Ks is 0.17 for orthologs against 0.71–0.88 for paralogs.

### D6. Cardoso-Moreira 2019 and 2020

**Citations.**
- Cardoso-Moreira M et al. Nature 2019;571:505–509, PMID 31243369.
- Cardoso-Moreira M et al. Cell Rep 2020;33:108308, PMID 33113372.
- [full text for both; supplementary tables not read]

**2019.**
- Early-development genes are less tolerant of loss-of-function variants (P < 1e-10).
- Pleiotropy correlates with pLI (r 0.29).
- **Genes with new trajectories are "as pleiotropic as" genes with conserved ones** (ED Fig. 10d).
- 93–96% of changes are organ-specific; rodents have more changes.
- Kidney is not in the change analysis, and dN/dS was not tested.

**2020.**
- 51% of genes (6,816/13,471) differ between human and mouse in at least one organ.
- Disease genes are depleted for differences (brain P 0.002, liver 9e-5).
- Stark breadth differences are depleted for essential genes (4% vs 11%, P 8e-6) and disease genes (16% vs 26%, P 0.02).
- Organ specificity correlates at r 0.86.

**Data.** evodevoapp.kaessmannlab.org.

### D7. Supporting papers

- **Chen J et al. Genome Res 2019;29:53–63** (PMID 30552105) [full text; supplement not read].
  - OU model, dX = σdB + α(θ − X)dt, with evolutionary variance σ²/2α, over 17 mammals.
  - 83% of genes fit stabilising selection better.
  - Variance against phyloP: r −0.25.
  - Kidney high-variance genes are enriched for oxidoreductases.
  - Sequence-conserved but expression-divergent genes are enriched for transcription regulators (e.g. NOTCH1).
  - Lineage-shift calls had an FDR of 54–78%.
  - The EVEE portal returned HTTP 502 on 2026-10-18; data at GEO GSE106077.
- **Sepp M et al. Nature 2024;625:788–796** (PMID 38029793) [full text, relevant sections].
  - Cerebellum: 23% of genes preserved, 17% diverged, 60% intermediate.
  - Presence/absence differences in 15.3% of genes.
  - **Gains are less constrained (LOEUF) and more cell-type-specific.**
  - **A random forest on exon length and GC could not predict the calls** (ED Fig. 11d): the artifact control to copy.
- **Yakubovsky O et al. Nature 2026;653:1148–1157** (PMID 41986723) [full text].
  - Human liver Visium against mouse.
  - Portal-bias correlation between species is **Spearman R = 0.15** (P = 5.2e-25).
  - VNN1, THBD and INMT are pericentral in mouse but periportal in human. GLUL, ADH4, ALB and GLS2 are conserved.
  - 12,410 1:1 orthologs from BioMart. **No dN/dS, gene-age or paralog analysis.**
  - *Vnn1* is also on our reversed PT list.
- **Harnik Y et al. Nature 2024;632:1101–1109** [abstract and extended data legends only]. Human vs mouse enterocyte zonation. No claims taken.

**Not found.** No study relates conservation of a zonation or axial pattern (liver lobule, villus, nephron) to dN/dS, gene age or paralog status.

### What predicts conserved vs changed expression

| Predictor | Evidence | Size |
|---|---|---|
| Protein constraint (dN, dN/dS, % identity) | weak, consistent | r 0.15–0.19 (arrays); ρ 0.27–0.38 (RNA-seq), 0.16 after 11 covariates; Strand quartiles about 0.53 vs 0.33; a shape-only measure showed nothing (Liao & Zhang) |
| Duplication | strongest categorical predictor | similar profiles 13.7% → 7.4% → 1.9% |
| Constraint, essentiality, disease | predicts conservation | LOEUF (Sepp); 4% vs 11% essential; disease depletion |
| Expression level | confounder | abolishes the low-expression correlation |
| Pleiotropy / τ | not reliable | no difference in Cardoso-Moreira 2019 |

### Where to get human–mouse dN/dS (tested 2026-10-18)

- **Ensembl status:**
  - Ensembl **116** (June 2026) is the **final release on the legacy platform**.
  - `www.ensembl.org/biomart/martservice` returns HTTP 308 to `jun2026.archive.ensembl.org`.
  - REST `/info/software` returns release 116.
- **dN/dS is gone from current Ensembl.**
  - `hsapiens_gene_ensembl` has 14 `mmusculus_homolog_*` attributes and no `_dn` or `_ds`. The same holds in releases 114, 110 and 105.
  - The Compara 116 README says "dn : non-synonymous mutation rate (currently unused)".
- **Last release with values: Ensembl 99** (Nov 2019).
  - Range-read headers: releases 96 and 99 have dn/ds filled; 100 and 104 have NULL.
  - Full file: `https://ftp.ensembl.org/pub/release-99/tsv/ensembl-compara/homologies/homo_sapiens/Compara.99.protein_default.homologies.tsv.gz` (170 MB compressed; not downloaded).
  - Columns: `gene_stable_id, protein_stable_id, species, identity, homology_type, homology_gene_stable_id, homology_protein_stable_id, homology_species, homology_identity, dn, ds, goc_score, wga_coverage, is_high_confidence, homology_id`.
  - Filter to `homology_species == mus_musculus` and `homology_type == ortholog_one2one`.
- **BioMart that still serves dN/dS:** the Ensembl 80 archive, `https://may2015.archive.ensembl.org/biomart/martservice`. It is marked "will be retired soon". The Ensembl 75 archive returned 503.
- **The query I ran** (POST `query=` XML):
  - dataset `hsapiens_gene_ensembl`;
  - filter `hgnc_symbol` = SLC5A2, APOC3, VNN1, SLC34A1, LRP2, NOTCH1, PDK4, CUBN;
  - attributes `ensembl_gene_id, external_gene_name, mmusculus_homolog_ensembl_gene, mmusculus_homolog_orthology_type, mmusculus_homolog_orthology_confidence, mmusculus_homolog_perc_id, mmusculus_homolog_dn, mmusculus_homolog_ds`.
  - Response: HTTP 200, 805 bytes, 0.8 s, all 8 `ortholog_one2one`.

  | Gene | dN | dS | dN/dS |
  |---|---|---|---|
  | APOC3 | 0.2672 | 0.6025 | 0.44 |
  | CUBN | 0.1955 | 0.6557 | 0.30 |
  | VNN1 | 0.1607 | 0.6764 | 0.24 |
  | LRP2 | 0.1401 | 0.6720 | 0.21 |
  | PDK4 | 0.0402 | 0.6264 | 0.064 |
  | SLC34A1 | 0.0469 | 0.8047 | 0.058 |
  | SLC5A2 | 0.0503 | 0.9034 | 0.056 |
  | NOTCH1 | 0.0395 | 1.0616 | 0.037 |

  - The Ensembl 80 confidence flags disagree with 116 (four of five 0s are 1 in 116), so do not filter on the Ensembl 80 flag.
  - Anecdote, not a test: reversed APOC3 and VNN1 have high dN/dS, while reversed PDK4 and NOTCH1 are very low. That fits Chen's sequence-conserved, expression-divergent regulators.
- **Ensembl 116 attributes for ortholog QC:**
  - `mmusculus_homolog_associated_gene_name, _orthology_type, _subtype` (the speciation node, not gene age);
  - `_orthology_confidence, _perc_id, _perc_id_r1, _goc_score, _wga_coverage`.
- **Paralogs:** `hsapiens_paralog_associated_gene_name, hsapiens_paralog_subtype, hsapiens_paralog_orthology_type, hsapiens_paralog_perc_id`, and `mmusculus_paralog_*` on the mouse dataset. Example: VNN1 → VNN2 (Mammalia, 62.6%).
- **Constraint:** gnomAD v4.1 `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv` (95.5 MB). Columns include `lof.pLI` and `lof.oe_ci.upper` (LOEUF).
- **Gene age:** GenTree (gentree.ioz.ac.cn responds; paper not read) or Ensembl paralog subtype.
- **Recommendation:**
  - Ensembl 99 FTP for dN and dS, with dS ≤ 2× the median;
  - Ensembl 116 for % identity (sensitivity) and for GOC/WGA/confidence QC;
  - Ensembl 80 BioMart only for quick lookups.

### Testable predictions

**Primary contrast:** reversed (74) against conserved (189). Both are graded in all four specimens, so detection and reliability are matched.

**Covariates:**
- log mean count and detection per species;
- |E|;
- CDS length and GC;
- probes per gene in each panel;
- GOC/WGA;
- capped dS.

**Inference:** logistic regression, plus a permutation of the predictor within expression × detection deciles (2,000 permutations).

1. **P1 · Constraint (primary).**
   - Test A: log dN/dS is higher in reversed genes. Expect AUC ≈ 0.6, which gives z ≈ 2.5 for 74 vs 189 (z ≈ 2.0 at AUC 0.58). Sensitivity: % identity.
   - Test B, Strand-style: reliability-corrected r(E_h, E_m) per dN/dS quartile among genes graded in at least one species. Corrected Q1 − Q4 ≥ 0.1, beyond the permutation band.
   - **Against:** adjusted AUC < 0.55 with a CI that includes 0.5, **and** Q1 − Q4 inside the null band.
2. **P2 · Duplicability.** A young within-species paralog (Mammalia or younger) predicts reversal.
   - Test: CMH OR across expression deciles.
   - **Against:** OR < 1.2 or a CI that includes 1. Report the count first, because power may be low.
   - Distinct from notebook 88, which covered non-1:1 families.
3. **P3 · LOEUF.** Higher in reversed genes. **Against:** AUC < 0.55.
4. **P4 · Evolvability.** Kidney OU variance (Chen) is higher in reversed genes. **Against:** AUC < 0.55. Blocked until the data are reachable.
5. **P5 · Cross-organ lability.** Genes with reversed liver zonation (Yakubovsky, Zenodo) are over-represented among our reversed PT genes. Seeds: Vnn1, and Dpp4 (human-only here, human-pericentral in liver). **Against:** the Fisher OR CI includes 1.
6. **P6 · Pleiotropy (exploratory).** The prior is no difference. Add a positional τ to ask whether reversed genes keep their localisation.

**Artifact control (after Sepp).** A random forest using only probe count, CDS length, GC and detection must fail to separate reversed from conserved (AUC ≤ 0.6). If it succeeds, panel artifacts are live and P1–P5 cannot be read biologically.

---

## The three most promising new analyses, ranked

### 1. Neighbour-contamination audit of the positional, injury and reversed-gradient calls (A-ring → A-cov/ORM, with A-erode and a planted arm)

**Why first.**
- It is the only one of the three that could overturn claims we already make.
- There is direct precedent for spillover creating false DE along a positional axis (Mitchel's crypt/villus goblet cells), and for marker decays that look like gradients (SpotClean).
- Within-specimen shape tests do not cancel contamination that varies with position (*our inference*). PT neighbours change along S1→S3 and with IR.
- It reuses notebook 71's gene models and notebook 100's planting code.

**Prediction.**
- The fitted bleed is small and condition-independent: σ about 3–6 µm, adjacent density about 1–4% of source (Oliveira's numbers).
- ≥ 70% of each robust pathway's effect, and of the 74 reversed-gradient genes, survives adding neighbour expression E_gu and dropping ORM > 1 genes.
- Flagged genes are not over-represented among calls.
- Contamination planted at 1–4% with the fitted kernel, varying with x or with condition and carrying no true signal, keeps the relabeling-calibrated FDP ≤ 0.1.

**Result that counts against it.** Any of:
- more than 20% of primary calls or reversed genes flagged;
- Fisher OR > 2 for flagged genes among calls;
- both IR mice above both controls by more than 1.5× in bleed σ or floor;
- the eroded-vs-dilated ratio R_g rising toward S3 or in IR;
- planted FDP > 0.1 at 2% contamination.

### 2. Sequence constraint and duplicability against human–mouse gradient reversal (D-P1, D-P2, plus the artifact control)

**Why.**
- It is cheap: about 170 MB from Ensembl 99 and 95 MB from gnomAD, applied to the existing notebook 90 categories.
- No precedent exists for axial or zonation patterns.
- It gives the species result an evolutionary reading.
- Liver (Yakubovsky R = 0.15) shows zonation is also poorly conserved there, with Vnn1 reversed in both organs.

**Prediction.**
- Reversed genes have higher log dN/dS than conserved ones (adjusted AUC ≈ 0.6).
- The Strand-style reliability-corrected r(E_h, E_m) rises by ≥ 0.1 from the highest to the lowest dN/dS quartile.
- Genes with a young paralog are over-represented among reversed genes (CMH OR ≥ 1.2).

**Result that counts against it.**
- Adjusted AUC < 0.55 with a CI that includes 0.5, **and** Q1 − Q4 inside the permutation band. This closes the hypothesis. Liao & Zhang's null for shape divergence on arrays makes this a real possibility.
- Separately, if the probe/length/GC/detection random forest separates reversed from conserved (AUC > 0.6), the categories may be panel artifacts and no biological reading is allowed.
- Power is modest (z ≈ 2.5 at AUC 0.6). Test B carries most of the power.

### 3. Position-resolved differential variability under injury (C\*)

**Why.**
- It asks a question our mean-based pipeline cannot: are injured tubules unequal, rather than uniformly shifted?
- It is untested in kidney.
- It builds on `nb_gam` (α as a function of x and m), limma-style moderation and cameraPR.

**Prediction.**
- IR raises residual tubule-to-tubule dispersion τ² of the INJ, CYC and TPT programs, in both IR mice, at positions where the IR mean effect is intermediate (excess τ² ∝ p(1−p)).
- The responding fraction falls where dispersion rises (Martinez-Jimenez's pattern).

**Result that counts against it.** Any of:
- set ε contrasts inside the relabeling spread;
- a global shift affecting all genes equally (technical);
- the effect absorbed by the 1/m unit-size or g′²σ_x² coordinate-error terms;
- present in only one IR mouse;
- reversed by adding the dropped deep-S3 units.

A null after residualisation would support the graded-shift reading from the existing D4 result.

**Cheap first look:** per specimen × bin, the MAD of the INJ-score residual and the injury-gene detection fraction.

**Not ranked: B (DOMINO, propagation).**
- Network calls look specific only under gene-shuffle nulls.
- Most active-module methods learn degree.
- The planted benchmark disfavours DOMINO.
- Propagation adds little over a degree-only baseline.
- Our cameraPR with relabeling already handles correlation, which those nulls cannot.
- If tried at all: DOMINO with a fixed top 300 under our relabelings and a rewired network, or ANUBIX as one extra row in notebook 71's table.

---
