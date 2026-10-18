# Do the mouse AKI calls hold in human AKI? A KPMP benchmark (notebook 98)

**The question.** Notebook 94 showed which of our mouse AKI calls replicate in another lab's mouse IRI. The paper's
AKI section has a translational question as well: which positional injury changes of the mouse PT appear in human
acute kidney injury?

The KPMP/HuBMAP atlas (Lake et al. 2023, Nature, PMID 37468583; already local) has single-nucleus PT profiles. They
come from AKI biopsies and reference tissue, with PT-S1, S2 and S3 labels and injured states (aPT, dPT, cycPT). This
notebook scores the same call lists on human AKI, plus notebook 95's signalling finding.

Protocol frozen on 2026-10-18, before any expression comparison of KPMP donors. Only the metadata (donors, disease
category, region, sex, PT subclass counts) had been tabulated.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/98_pt_human_aki_benchmark.ipynb`; logic `98.human_aki.1`; seed 98.
- Results `results/pt_human_aki_benchmark/`.
- No new module: it reuses `external_replication`, `species_level.limma_trend_species`, `regulator_activity` and
  `pathway_pipelines.camera_pr`.
- **Data:** `data/external/kpmp_lake2023/lake2023_snRNA_v1.5.h5ad` (CELLxGENE, raw counts in `raw/X`, human symbols in
  `var.feature_name`).
- **Our calls and programs:** exactly notebook 94's lists and claim-check programs, read the same way. The one change
  is in S-P80: notebook 80's L23 programs (Lake 2023's aPT lists) are excluded, because they were derived from these
  data. Kirita's K-* programs are kept, since they are independent of KPMP.
- **PROGENy:** native human weights (`data/external/regulator_networks/progeny_human.tsv`, 500 genes per pathway), the
  same file as notebook 95.

## 2 · Donors, nuclei and statistics (fixed now)

- **Nuclei:**
  - `subclass.l1 == 'PT'`, excluding `dPT/DTL` (ambiguous with thin limb);
  - region C (cortex), because AKI biopsies are cortical;
  - donors with ≥ 100 such nuclei.
- **Groups:**
  - AKI: `disease_category == 'AKI'`;
  - reference: Healthy_reference_tissue, Healthy_stone_donor and Tumor_nephrectomy.
  - Sensitivity: reference = Healthy_reference_tissue only.
- **Genes:**
  - human symbols map to mouse through the one-to-one ortholog map (`ortholog_map_used.csv`);
  - E_h = (notebook 94's candidate universes) ∩ mapped genes, with whole-PT CPM ≥ 1 in ≥ 3 donors;
  - mitochondrial genes are out.
- **Level:**
  - whole-PT pseudobulk per donor;
  - log2 CPM, voom-style;
  - limma-trend with robust eBayes, AKI against reference, centred by the median (`limma_trend_species`).
- **Location:**
  - per donor, d = log2 CPM(PT-S3) − log2 CPM(PT-S1);
  - only donors with ≥ 20 nuclei in each;
  - the same limma fit on d.
- **Sets:** camera as in notebook 94 on the level t, the location t and |location t|.
- Sex and age are not modelled. This is a known limit of a two-group fit, and the groups' sex and age are reported.

## 3 · Scoring and rules (notebook 94's, with one contrast instead of time points)

- Lists, axes, internal signs, matched nulls (1,000), rule 1 (replicates) and rule 2 (comparison margin 0.10, n ≥ 20)
  are exactly notebook 94's.
- Rule 3 (stability across time points) is replaced. A comparison is **stable** if it also holds under the
  healthy-reference-only sensitivity contrast.
- **The F-test verdict:** "replicates", "partly" or "does not replicate", as in notebook 94's rule 4.
- **Claim checks:** K1–K9 as in notebook 94, on the human contrast, plus one:
  - **K10:** PROGENy EGFR, MAPK and PI3K activity up in human AKI. Activity is ulm on the level t with 1,000
    permutations within expression deciles (notebook 95's estimator). It holds when q ≤ 0.05 and z > 0 for each.
- **Cross-species summary:** for each list, its precision in human against its precision in Kirita at 12 h (94):
  which mouse calls generalise to human AKI.

## 4 · Outputs and figures

- **Tables:**
  - `donors.csv`, `universe.csv`, `external_genes.csv.gz`, `external_sets.csv.gz`;
  - `list_scores.csv`, `comparisons.csv`, `claim_checks.csv`;
  - `cross_species.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 98a: precision against the null band, human, every list;
  - 98b: human against mouse (Kirita 12 h) precision per list;
  - 98c: claim checks, human and mouse side by side.

## 5 · Caveats

- Human AKI here is a mix of causes, biopsy times and severities. The mouse model is a single, timed ischemic insult.
- Biopsies are small cortical cores. S3 sits partly in the medulla, so human S3 is under-sampled.
- The reference group mixes healthy reference, stone donors and tumour nephrectomy.
- KPMP profiles nuclei; we profile Visium HD tubule units.
- Agreement is a lower bound on shared biology.
