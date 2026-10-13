# Pathways a segment-level analysis puts in one segment: where they sit inside it (notebook 84)

Frozen on 2026-10-12, before any fit or within-segment position of notebook 84 was computed. Before freezing, I had
seen only:
- notebook 83's draft of 16 published programs;
- the coordinate records of notebooks 53, 72 and 74.

No pathway curve inside a segment had been computed.

## Question

A segment-resolution analysis calls pathways "enriched in PTS3" (or PTS1, PTS2). At that resolution these pathways
are one class. Along a continuous coordinate, do they sit at different places inside the segment (its first,
middle or last third)? And do they differ more than random groups of the same segment's genes?

## Why this design

- **The discrete call comes from data that are not ours.** Chen et al. 2021 microdissected mouse PTS1, PTS2 and PTS3
  with replicates (JASN; PMID 33769951; GEO GSE150338). Enrichment in their segment contrast is what a segment-level
  study reports.
- **Positions are read on LRS-Z, not the package coordinate.**
  - By the frozen C14 rule, within-segment results on the package coordinate are descriptive: notebook 74 found 10
    of 12 cells noise-dominated, and only mouse S3 was not.
  - LRS-Z (notebook 53's zone-anchored landmark ratio) has 0 of 12 noise-dominated cells on the same units, and
    gene-fold agreement of 0.93–0.97 (notebook 72; issue 3).
  - LRS-Z has fold versions that drop one gene fold's landmarks. So each gene is placed on a coordinate that did not
    use it (fold protection).
  - The package curves (notebook 71) enter as a concordance check.
- **The null keeps every gene's own position.** Gene identities are permuted within the segment's gene pool,
  stratified by expression. The question is whether pathway membership groups genes by position more than chance
  does.

## Labels

- Notebook: `analysis/notebooks/84_pt_within_segment_pathways.ipynb`
- Logic: `84.within_segment_pathways.1`
- Module: `pseudospace/within_segment.py`, tested in `tests/test_within_segment.py`
- Results: `results/pt_within_segment_pathways/` (`tables/`, `figures/`, `run_record.json`)

## Inputs (read only)

| Source | What 84 reads |
|---|---|
| Chen 2021 | `data/external/segment_microdissection/GSE150338_Mouse_TPM_Replicates.txt.gz` (sha256 `daf8b6f4…5488`). PT samples: PTS1 5, PTS2 8, PTS3 6 (`m{mouse}_{segment}_{a,b,c}`). Non-PT structures: DTL1, DTL2, DTL3, ATL, MTAL, CTAL, DCT, CNT, CCD, OMCD, IMCD, GL |
| 70, 71 | `package_inputs.hvm_inputs` (mt-free, as 71 and 81): counts, library, labels, 11,055 universe genes, 1,513 gene sets. 71's stage cache: NB specimen curves (the package arm) |
| 72 | `transitions.csv`: the package segment units |
| 53 | `pt_reconstruction_v3/coordinates/lrs_z.csv`: `position` and `fold_0`–`fold_2`, matched to units by `old_structure_id` (set M, as 81) |
| 74 | the base integration features (pass 1 ∪ pass 2 `highly_variable_for_harmony`) and the 15 panel and anchor genes: package-arm fold protection (74's rule) |

**Guards (not cached).**
- G1: the Chen file's sha256.
- G2: `hvm_inputs` reproduces 71's universe, sets and units.
- G3: set M is the same size as in 81 (11,045 units).
- G4: the LRS-Z fold assignment `sha256_folds(gene, 3)` equals 53's rule; fold k omits only fold-k landmarks.

## D · The segment-level call (Chen 2021)

- **D1 · Segment contrast.**
  - Genes: Chen symbols in 71's universe; duplicated symbols keep the row with the highest mean PT TPM. Genes need
    mean TPM ≥ 1 in at least one of PTS1, PTS2 and PTS3.
  - Scale: log2(TPM + 1). Technical replicates (a, b, c) of one mouse and segment are averaged, giving 3, 5 and 3
    samples.
  - For segment k: the effect is mean_k − (mean_j + mean_l)/2. The residual variance comes from the three-group
    one-way model (8 df).
  - The moderated t is `program_localization.moderated_t`: limma squeezeVar with a trend on the mean log2 TPM. The
    effect is divided by √(1/n_k + (1/n_j + 1/n_l)/4).
  - Samples are treated as independent: no mouse blocking, because mouse IDs are not known to be shared across
    segments.
- **D2 · Segment pathways.**
  - cameraPR (`pathway_pipelines.camera_pr`, inter-gene correlation 0.01) of 71's 1,513 sets on the D1 t of
    segment k.
  - A "segment-k pathway" is direction Up with BH FDR ≤ 0.05.
- **D3 · Core genes.** A segment-k pathway's core is its members that meet all of:
  - t ≥ 2 and effect ≥ 0.5 log2 in segment k;
  - PT-specific: mean TPM in PTS-k above the mean TPM of every non-PT structure (the bound against neighbour signal,
    issue 24);
  - usable in both mice: detection ≥ 5% on set M units.

  A pathway enters the analysis with at least 5 core genes. The **pool** of segment k is every universe gene that
  meets these criteria, whatever its pathway.
- **D4 · Representatives.**
  - Pathways are taken in cameraPR p order (ties by name).
  - A pathway becomes a representative unless its core overlaps an earlier representative's core by ≥ 0.6 (overlap
    = |A ∩ B| / min(|A|, |B|), the project's redundancy measure).
  - Each representative's theme lists the pathways it absorbed.
  - Tests use representatives; every pathway is described.

## C · Where each pathway sits inside the segment (ours)

- **C1 · Gene curves on LRS-Z, fold-protected.** Per mouse, unit score = log2((count + 0.5) / library × the median
  library). Each score is then:
  - centred per gene;
  - reduced by the unit's median usable gene (81's `UNIT_TYPICAL` rule).

  Gene g sits on LRS-Z fold `sha256_folds(g)`. Per specimen, its curve is a least-squares cubic B-spline with 12
  basis df (9 interior knots at the fold coordinate's quantiles), evaluated on that fold's common grid.
  - **Segment units per fold:** 72's pooled logistic transitions on that fold's positions.
  - Curves are interpolated to one grid of 151 points over 0–3 segment units. Points outside the fold grid's range
    are left empty.
- **C2 · Pathway curve and within-segment profile.**
  - Pathway curve: the median over core genes (log2).
  - In segment k (grid points with k ≤ u < k + 1), the profile v is the means of the first, middle and last thirds,
    minus the segment mean.
  - Also recorded:
    - range r, the maximum minus the minimum of the curve inside the segment;
    - the peak position, the argmax inside the segment;
    - the peak third, the argmax of v.
  - "Flat inside the segment" means r < 0.25 log2.
- **C3 · Confirmatory test, one per segment (mouse S1, S2, S3; Holm over the three).**
  - Statistic: T_k = Σ over representatives of (v_p^A − v̄^A) · (v_p^B − v̄^B), where A and B are the two mice and
    v̄ is the mean over representatives. This is the between-pathway spread of the within-segment profile, counted
    only where the mice agree.
  - Null: 2,000 permutations (seed 84) of gene identities within the segment-k pool, stratified by tertile of log
    mean count over mouse units. The same permutation applies in both mice, and every core maps to its permuted
    genes.
  - p = (1 + #{T_null ≥ T}) / 2,001.
  - At least 3 representatives are needed; otherwise "not testable".
- **C4 · Reportable position (descriptive), per representative and segment.** All of:
  - not flat in either mouse on LRS-Z;
  - the same peak third in both mice on LRS-Z;
  - the package arm agrees where it can be evaluated.

  The package arm uses 71's curves under 81's rule (`reference_curve` centring, `program_curve` median), on core
  genes outside the base integration features and the 15 panel and anchor genes, with ≥ 5 such genes. "Agrees" means
  the same peak third in both mice; otherwise the position is "LRS-Z only".
- **C5 · Human (descriptive).** Applies C1–C3 to the two human sections for the same mouse-defined representatives.
  It reports T and p, but no claim: one donor, and the pathways come from mouse data.

## What would support the claim, and what would not

- **Claim, for segment k:** "pathways that a segment-level analysis calls segment-k-enriched sit at different places
  inside segment k". It needs:
  - Holm-adjusted p ≤ 0.05 for T_k;
  - at least two representatives with reportable positions in different thirds.
- **Adds nothing in segment k:** T_k is not significant. The figure is then drawn and captioned "no within-segment
  separation".
- **Not testable:** fewer than 3 representatives. Too few PTS-k pathways at the segment level is itself the reading.

## Sensitivities (descriptive)

- Core without the PT-specific filter.
- Core threshold t ≥ 3.
- Spline basis df 6 (81's basis).
- Unstratified permutation.
- Every pathway in place of representatives, for the descriptive table only.

## Figures (no words inside the plot area)

- **84a–c, one per segment.**
  - Left: Chen's segment means of each representative's core (median over core genes of centred log2 TPM, PTS1,
    PTS2, PTS3). Here they all look alike.
  - Right: the LRS-Z curves (mean of the two mice) as a heatmap over 0–3 segment units, with rows ordered by peak
    position inside segment k and a mark at each row's peak.
- **84d.** Line plots of representatives with reportable positions in different thirds: the earliest, the middle
  and the latest. Each specimen gets its own line.

## Outputs

- `chen_segment_t.csv`
- `segment_pathways.csv` (cameraPR, cores, themes)
- `positions.csv` (per representative × specimen × arm)
- `tests.csv`
- `sensitivities.csv`
- `how_the_protocol_was_read.csv`
- `run_record.json`
