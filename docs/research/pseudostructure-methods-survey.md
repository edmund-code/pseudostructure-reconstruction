# Literature-guided design for repeated PT reconstruction

Updated 2026-10-04. Bounded, targeted survey of trajectory inference, reference
alignment and optimal transport, not an exhaustive review or a novelty audit.
Searches used public method names and generic scientific questions only. Parallel
CLI/authentication was unavailable, so retrieval used the available web search and
primary-source pages instead. No private profiles were sent to an external service.

| Primary source | Relevant method or evidence | Implication for our experiment |
| --- | --- | --- |
| [Haghverdi et al., DPT (2016)](https://www.mdc-berlin.de/research/publications/diffusion-pseudotime-robustly-reconstructs-lineage-branching), DOI 10.1038/nmeth.3971 | Diffusion-based distances order profiles using a chosen starting state. | Required graph baseline. Use a training-only S1 root; do not claim a graph distance is physical nephron length. |
| [Faure et al., scFates (2023)](https://pubmed.ncbi.nlm.nih.gov/36394263/), DOI 10.1093/bioinformatics/btac746; [official API](https://scfates.readthedocs.io/en/latest/api.html) | Principal graphs, pseudotime and downstream bifurcation analysis; the installed curve routine uses ElPiGraph. | Required nonbranching curve baseline. A curve is a legitimate backbone; novelty need not reside in inventing a new curve optimizer. |
| [Albergante et al., ElPiGraph (2020)](https://www.mdpi.com/1099-4300/22/3/296), DOI 10.3390/e22030296; [author API](https://elpigraph-python.readthedocs.io/en/stable/elpigraph.computeElasticPrincipalCurve.html) | Elastic principal graphs; API exposes initial node positions/edges and point weights. | Test an ordered anatomical initial polyline as a shared physical-object prior. Initialization is not enforced anatomical order, and generic curve flexibility does not establish longitudinal anatomy. Reuse the existing backend rather than inventing another optimizer. |
| [Street et al., Slingshot (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6007078/), DOI 10.1186/s12864-018-4772-0 | Cluster-level lineage structure initializes simultaneous principal curves; single-lineage initialization follows cluster centers instead of PC1. | Our anatomical initialization ablation has methodological precedent. Better initialization alone is not a distinct reconstruction contribution. |
| [Saelens et al., trajectory benchmark (2019)](https://www.nature.com/articles/s41587-019-0071-9), DOI 10.1038/s41587-019-0071-9 | Compared 45 methods on 110 real and 229 synthetic datasets; suitability depends on dimensions and topology. | Fix the nonbranching PT scope and compare matched inputs, prediction and stability; avoid claiming a universal trajectory winner. |
| [Campbell and Yau, PhenoPath (2018)](https://research.manchester.ac.uk/en/publications/uncovering-pseudotemporal-trajectories-with-covariates-from-singl/), DOI 10.1038/s41467-018-04696-6 | Latent pseudotime with covariate effects/interactions. | Separating a coordinate from specimen/domain effects has precedent; the physical repeated-object setting and validation contract must carry our contribution. |
| [Campbell and Yau, Order Under Uncertainty (2016)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005212), DOI 10.1371/journal.pcbi.1005212 | Probabilistic pseudotime and propagation of ordering uncertainty into downstream differential analysis. | Preserving assignment ambiguity has precedent. For repeated cuts, distinguish prediction at a coordinate mean from averaging predictions across possible positions; neither yields calibrated anatomical uncertainty automatically. |
| [Velten et al., MEFISTO (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8828471/), DOI 10.1038/s41592-021-01343-9 | Smooth/non-smooth factor variation and group alignment using supplied continuous covariates. | Useful rationale for separating reproducible spatial programs from amplitude variation; our latent position is not an observed MEFISTO covariate. |
| [Welch, Hartemink and Prins, MATCHER (2017)](https://link.springer.com/article/10.1186/s13059-017-1269-0), DOI 10.1186/s13059-017-1269-0 | Gaussian-process latent trajectories and quantile alignment to a uniform master coordinate integrate different molecular assays. | Shared latent alignment has precedent. Quantile matching can confound differences in sampling coverage with position in our cuts; do not import uniform occupancy as anatomical evidence. |
| [Schiebinger et al., Waddington-OT (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6402800/), DOI 10.1016/j.cell.2019.01.006 | Distributional couplings between known collection times, with growth/death modeling. | Do not replace collection time by unknown spatial position and inherit identifiability. An OT component must earn value beyond a mean reference. |
| [Nitzan et al., novoSpaRc (2019)](https://www.nature.com/articles/s41586-019-1773-3), DOI 10.1038/s41586-019-1773-3; [author documentation](https://novosparc.readthedocs.io/) | Spatial gene-expression reconstruction by transport under structural correspondence assumptions. | Soft assignment to a canonical anatomical object has precedent. Our winding PT excludes using tissue x/y distances as longitudinal-nephron geometry. |
| [Liu et al., PASTE2 (2023)](https://genome.cshlp.org/content/genome/early/2023/08/08/gr.277670.123.full.pdf), DOI 10.1101/gr.277670.123 | Partial fused Gromov–Wasserstein alignment allows incompletely overlapping tissue slices. | Partial matching is relevant to uneven PT coverage; slice alignment is not our task. Avoid transporting every population into a fixed full-axis occupancy. |
| [Klein et al., moscot (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11864987/), DOI 10.1038/s41586-024-08453-2 | Scalable, multimodal transport for temporal, spatial and spatiotemporal mappings. | Reuse an established OT implementation if transport helps; scale or multimodality alone does not identify an unknown anatomical coordinate. |
| [Alpert et al., cellAlign (2018)](https://www.nature.com/articles/nmeth.4628), DOI 10.1038/nmeth.4628; [author implementation](https://github.com/shenorrLabTRDF/cellAlign) | Dynamic-time-warping alignment of interpolated expression trajectories; local and global alignment differ. | Canonical reconstruction can use registration, but complete global alignment may force missing states to match. |
| [Sumanaweera et al., Genes2Genes (online 2024; issue 2025)](https://www.nature.com/articles/s41592-024-02378-4), DOI 10.1038/s41592-024-02378-4 | Distribution-aware gene alignments include matches, warps and mismatches. | Preserve molecular divergence rather than converting every disease/species change into positional warping; infer structure before downstream gene-specific comparisons. |
| [Chakraborty and Panaretos, functional registration (accepted author version)](https://arxiv.org/abs/1702.03556) | Identifiability of amplitude versus coordinate warping depends on structural conditions; roughness penalties alone can fail. | A smooth fitted atlas is not evidence that specimen shifts and position have been identified. Restrict nuisance variation and test ordering separately from numerical gauge. |
| [Panaretos and Zemel, point-process registration (2016)](https://arxiv.org/abs/1603.08691), DOI 10.1214/15-AOS1387 | Links repeated point-process registration, amplitude/phase variation and Wasserstein geometry, including over-registration. | Relevant precedent for repeated-object registration and a reason to inspect forced matching. Their observed point-process coordinate is not our unknown longitudinal coordinate, so the theory is not an identification guarantee here. |
| [Lause, Berens and Kobak, analytic Pearson residuals (2021)](https://link.springer.com/article/10.1186/s13059-021-02451-7), DOI 10.1186/s13059-021-02451-7 | Analytic count residuals use an exposure-dependent NB variance; the authors compare residual and variance-stabilizing preprocessing. | Our existing transform follows this precedent. For repeated cuts, test whether using one training-reference scale preserves placement under read perturbation; it changes noise weighting and is not the paper's Pearson transform. |
| [Hafemeister and Satija, regularized NB normalization (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6927181/), DOI 10.1186/s13059-019-1874-1 | Models sequencing depth and emphasizes preserving biological heterogeneity while controlling observation effects. | Exposure handling matters before latent-space fitting. Their single-cell results do not establish behavior for mixed-cell tubule aggregates; inspect actual count perturbations and gene transfer. |

Publisher/author abstracts, indexed primary-source excerpts and official method
APIs were checked across these entries. Accessible full text was read for
Slingshot and Genes2Genes; several publisher/PMC pages blocked full-text retrieval.
The table supports bounded design choices, not claims of complete full-text review.
PASTE2 is linked to the published article, rather than treating its earlier PMC
preprint as the final version. Each entry should be revisited for a manuscript's
precise implementation or historical claims.

## What is specific to our problem

Our prior is a shared ordered physical PT object, observed through unordered cuts
of many realizations with unknown nephron membership. A method suitable here can
combine existing tools: a canonical backbone, specimen-balanced molecular
reference distributions, restricted nuisance variation, partial-coverage-aware
projection and visible assignment ambiguity. These ingredients are proposed
design choices, not established novelty. The meaningful test is whether they
improve reproducible organization over DPT, scFates and coarse segments in the
actual healthy sections while preserving molecular deviations as deviations.

## Hypotheses and a real-data decision sequence

1. **Anatomical initialization:** endpoint or soft ordinal starts may improve the
   repeated atlas's held-out prediction/stability over PC1. Hold its likelihood
   and penalties fixed. Improved training convergence alone is insufficient.
2. **Reference reconstruction:** the canonical atlas may transfer withheld gene
   programs better than DPT/scFates coordinates learned from identical modeling
   genes. Use both control-transfer directions and rotated gene folds.
3. **Projection:** a distribution-based frozen projector may add value beyond
   nearest-neighbor coordinate interpolation. First compare all axes through an
   identical neighbor projector to distinguish axis learning from projection.
4. **Restricted variation:** if simple geometry transfers better than isotropic
   atlas likelihood, test training-only residual covariance or restricted specimen
   calibration on actual data, with unchanged query/gene benchmarks.
5. **Optional transport:** only afterward test whether partial/unbalanced
   distributional matching improves common programs or assignment stability.
   Do not impose matched abundance, full support, real time or biological drift.

Query labels, withheld response genes, x/y and glomerular depth never select the
molecular coordinate. Existing coarse labels retain upstream marker information;
withheld markers are not completely independent of that annotation. With two
healthy mice all comparisons remain descriptive. Continued use of the same
cohort for method development must be stated; it is not a new external test set.

## Follow-up after actual-data reference tests

The functional-registration abstracts above were retrieved from the authors'
arXiv records; their formal conditions and proofs have not been audited for this
application. They motivate a caution, not a theorem about our estimator. Large
numerical shifts across gene panels can coexist with similar global ordering;
within-segment ordering must also be checked. Subsequent specimen calibration
should freeze its learned parameters before evaluation, explicitly restrict which
molecular directions it can remove, and test uneven coverage. Smoothness or
transport alone must not turn molecular differences into apparent position.

The [PhenoPath primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6015076/)
is also available in full text. Its introduction and model overview confirm the
shared latent axis with covariate effects/interactions; that is relevant precedent
for the umbrella formulation, not evidence that our repeated PT coordinate is
anatomically identifiable.

## Observation scale follow-up

Primary text for Lause et al. was read through its analytic-residual derivation,
and Hafemeister/Satija through the introduction and normalization diagnostics;
this is bounded methods reading, not a full review of either paper. PMC retrieval
was intermittent; the Lause publisher full text was accessible. The proposed
fixed-reference rate representation is our algebraic ablation, not a claimed new
normalization method or an implementation of scTransform. Its conditional mean
is invariant to common count/exposure scaling before clipping, but sampling noise
remains heteroscedastic. Both representation and frozen atlas must earn useful
actual-data behavior, including within-segment gene programs and read stability.

This observation-scale experiment was prepared but not run. It is deferred after
the project goal checkpoint: a preprocessing improvement alone would not establish
repeated-structure reconstruction. The next specified test concerns groupwise
canonical inference with matched pooled trajectory baselines.
