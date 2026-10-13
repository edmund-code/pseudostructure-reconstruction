# Pathways a segment-level analysis puts in one segment: where they sit inside it (notebook 84)

Frozen on 2026-10-12, before any fit or within-segment position of notebook 84 was computed. Before freezing, I had
seen only:
- notebook 83's draft of 16 published programs;
- the coordinate records of notebooks 72 and 74.

No pathway curve inside a segment had been computed.

**Amendment, before any fit (user instruction, 2026-10-12).** The first frozen version (commit 0281964) read
positions on LRS-Z. The user ruled: "only use the elastic curve as the primary coordinate, not anything else". The
elastic curve, the package's stage-04 pseudospace, is now the only coordinate.

The elastic curve's within-segment weakness is handled with the project's own tools:
- 74's fold-protection rule;
- an expression- and integration-stratified null;
- the C14 reading.

No other coordinate enters. Nothing had been computed under the first version.

## Question

A segment-resolution analysis calls pathways "enriched in PTS3" (or PTS1, PTS2). At that resolution these pathways
are one class. Along the elastic curve, do they sit at different places inside the segment (its first, middle or
last third)? And do they differ more than random groups of the same segment's genes?

## Why this design

- **The discrete call comes from data that are not ours.** Chen et al. 2021 microdissected mouse PTS1, PTS2 and PTS3
  with replicates (JASN; PMID 33769951; GEO GSE150338). Enrichment in their segment contrast is what a segment-level
  study reports.
- **Positions come from the elastic curve:** notebook 71's NB specimen curves, in 72's segment units.
- **The null keeps every gene's own position.** Gene identities are permuted within the segment's gene pool,
  stratified by expression tertile and by whether the gene built the curve (an integration feature). The question is
  whether pathway membership groups genes by position more than chance does.
- **C14 frames the reading.** Notebook 74 found the elastic curve's fold noise larger than the within-segment spread
  in 10 of 12 cells. Only mouse S3 is not noise-dominated: species-mean fold noise ÷ spread 0.98, against 2.0 in
  mouse S2 and 4.2 in mouse S1. So, by the frozen rule:
  - every within-segment statement here is descriptive;
  - the test measures whether pathways separate inside a segment beyond random gene groups;
  - mouse S3 is the cell where the curve's position noise is below the within-segment spread.

## Labels

- Notebook: `analysis/notebooks/84_pt_within_segment_pathways.ipynb`
- Logic: `84.within_segment_pathways.1`
- Module: `pseudospace/within_segment.py`, tested in `tests/test_within_segment.py`
- Results: `results/pt_within_segment_pathways/` (`tables/`, `figures/`, `run_record.json`)

## Inputs (read only)

| Source | What 84 reads |
|---|---|
| Chen 2021 | `data/external/segment_microdissection/GSE150338_Mouse_TPM_Replicates.txt.gz` (sha256 `daf8b6f4…5488`). PT samples: PTS1 5, PTS2 8, PTS3 6 (`m{mouse}_{segment}_{a,b,c}`). Non-PT structures: DTL1, DTL2, DTL3, ATL, MTAL, CTAL, DCT, CNT, CCD, OMCD, IMCD, GL |
| 70, 71 | `package_inputs.hvm_inputs` (mt-free, as 71 and 81): counts, labels, 11,055 universe genes, 1,513 gene sets. 71's stage cache: the condition model's NB specimen curves (101-point grid), convergence and separation flags |
| 72 | `transitions.csv`: the pooled S1→S2 and S2→S3 cuts and segment units |
| 70 run | the base integration features (pass 1 ∪ pass 2 `highly_variable_for_harmony`; 74's rule) and the 15 panel and anchor genes |

**Guards (not cached).**
- G1: the Chen file's sha256.
- G2: `hvm_inputs` reproduces 71's universe, sets and units, and 71's cached condition model reproduces its
  T_spatial (as 81).
- G3: the logistic transitions reproduce 72's (as 81).

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
  - usable in both mice: 71's curve converged and not separated, detection ≥ 5% (81's rule).

  A pathway enters the analysis with at least 5 core genes. The **pool** of segment k is every universe gene that
  meets these criteria, whatever its pathway.
- **D4 · Representatives.**
  - Pathways are taken in cameraPR p order (ties by name).
  - A pathway becomes a representative unless its core overlaps an earlier representative's core by ≥ 0.6 (overlap
    = |A ∩ B| / min(|A|, |B|), the project's redundancy measure).
  - Each representative's theme lists the pathways it absorbed.
  - Tests use representatives; every pathway is described.

## C · Where each pathway sits inside the segment (elastic curve)

- **C1 · Pathway curve.** Per mouse:
  - Gene curves are 71's NB specimen curves, centred by their mean and by the specimen's median usable gene at each
    grid point (`pf.reference_curve`).
  - The pathway curve is the median over core genes (`pf.program_curve`, log2), on 72's segment units (grid
    0.11–2.87).
- **C2 · Within-segment profile.**
  - In segment k (grid points with k ≤ u < k + 1), the profile v is the means of the first, middle and last thirds,
    each relative to their mean. S1's first third starts at 0.11 and S3's last third ends at 2.87.
  - Also recorded:
    - range r, the maximum minus the minimum of the curve inside the segment;
    - the peak position, the argmax inside the segment;
    - the peak third, the argmax of v.
  - "Flat inside the segment" means r < 0.25 log2.
- **C3 · Test, one per segment (mouse S1, S2, S3; Holm over the three).**
  - Statistic: T_k = Σ over representatives of (v_p^A − v̄^A) · (v_p^B − v̄^B), where A and B are the two mice and
    v̄ is the mean over representatives. This is the between-pathway spread of the within-segment profile, counted
    only where the mice agree.
  - Null: 2,000 permutations (seed 84) of gene identities within the segment-k pool. They are stratified by tertile
    of log mean count over mouse units × integration feature (yes/no). The same permutation applies in both mice,
    and every core maps to its permuted genes.
  - p = (1 + #{T_null ≥ T}) / 2,001.
  - At least 3 representatives are needed; otherwise "not testable".
- **C4 · Reportable position (descriptive), per representative and segment.** Both of:
  - not flat in either mouse;
  - the same peak third in both mice.

  A further column, "fold-protected", records whether the same peak third holds in both mice on core genes outside
  the base integration features and the 15 panel and anchor genes (74's rule; ≥ 5 genes needed).
- **C5 · Human (descriptive).** Applies C1–C3 to the two human sections for the same mouse-defined representatives.
  It reports T and p, but no claim: one donor, the pathways come from mouse data, and human S2 is the weak cell.

## What would support the statement, and what would not

- **Descriptive statement, for segment k:** "pathways that a segment-level analysis calls segment-k-enriched sit at
  different places inside segment k along the elastic curve". It needs:
  - Holm-adjusted p ≤ 0.05 for T_k;
  - at least two representatives with reportable positions in different thirds.

  Under C14 this stays descriptive. It is strongest in mouse S3, the one cell that is not noise-dominated.
- **Adds nothing in segment k:** T_k is not significant. The figure is then drawn and captioned "no within-segment
  separation".
- **Not testable:** fewer than 3 representatives. Too few PTS-k pathways at the segment level is itself the reading.

## Sensitivities (descriptive)

- Core without the PT-specific filter.
- Core threshold t ≥ 3.
- Cores outside the integration features (74's fold protection), for the test as well as the positions.
- Unstratified permutation.
- Every pathway in place of representatives, for the descriptive table only.

## Figures (no words inside the plot area)

- **84a–c, one per segment.**
  - Left: Chen's segment means of each representative's core (median over core genes of centred log2 TPM, PTS1,
    PTS2, PTS3). Here they all look alike.
  - Right: the elastic-curve pathway curves (mean of the two mice) as a heatmap over the PT in segment units, with
    rows ordered by peak position inside segment k and a mark at each row's peak.
- **84d.** Line plots of representatives with reportable positions in different thirds: the earliest, the middle
  and the latest. Each specimen gets its own line.

## Outputs

- `chen_segment_t.csv`
- `segment_pathways.csv` (cameraPR, cores, themes)
- `positions.csv` (per representative × specimen)
- `tests.csv`
- `sensitivities.csv`
- `how_the_protocol_was_read.csv`
- `run_record.json`
