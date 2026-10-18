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

## Notebook 98 results

`analysis/notebooks/98_pt_human_aki_benchmark.ipynb` (logic `98.human_aki.1`) writes to
`results/pt_human_aki_benchmark/`. It followed the protocol frozen in 1727af4, and its runtime is 35 s. Notebook 94's
read-only call lists passed their size guard (G4).

**How the protocol was read.**
- Notebook 94's scoring code is used unchanged, except:
  - two contrasts replace the time points;
  - there are no assigned-nuclei axes;
  - S-P80 drops Lake's L23 lists.
- Duplicated human symbols are summed.
- K10's strata are 91's expression deciles.
- One post hoc section (9.3).

**Donors.**
- AKI: 14 donors, 17,353 cortical PT nuclei.
- Reference: 7 healthy reference, 4 stone donor and 3 tumour nephrectomy donors.
- Location needs ≥ 20 S1 and S3 nuclei, which leaves 10 AKI against 13 reference donors (primary) and 10 against 7
  (healthy only).
- E_h has 10,470 genes and 1,484 library sets.
- Externally at q ≤ 0.05: 2,961 genes change level, but only 42 change location (9 in the healthy-only contrast).
  Location in human biopsies is close to unmeasurable.

**PT states (post hoc 9.3; share of cortical PT nuclei):**

| Group | S1 | S2 | S3 | aPT | dPT | cycPT |
|---|---|---|---|---|---|---|
| AKI | 0.11 | 0.06 | 0.09 | **0.58** | 0.15 | 0.013 |
| Healthy reference | 0.14 | 0.03 | 0.23 | 0.06 | **0.54** | 0.001 |
| Stone donor | 0.42 | 0.09 | 0.12 | 0.21 | 0.15 | 0.003 |
| Tumour nephrectomy | 0.16 | 0.07 | 0.06 | 0.37 | 0.33 | 0.004 |

The reference group is itself dominated by injured-state nuclei, so the frozen whole-PT contrast partly compares adaptive
with degenerative PT.

**Agreement with our L** (genes, ρ):
- 0.229 (primary) and 0.209 (healthy only), against 0.524 for Kirita at 12 h;
- 0.17 on identity-retaining nuclei (post hoc);
- location ρ is about 0;
- set level z agree at ρ 0.48.

**Lists on their primary axis** (human; enrichment over the matched null):

| List | Human | Mouse IRI 12 h (94) |
|---|---|---|
| G-L91 · level | 743/2928 (1.38, R) | 1.67, R |
| G-C71 · level | 1084/5256 (1.13, R) | 1.36, R |
| S-L91 · level | 30/538 (1.55, R) | 1.74, R |
| S-C71 · level | 17/459 (1.06) | 1.50, R |
| S-P80 · level | **13/22 (19.9, R)** | 2.05, R |
| G-M91, S-M91, S-N93 · location | 1/327, 0/236, 0/133 | 7.9, 5.3, 5.0, all R |
| G-F71 · unsigned location | 27/4296 (1.42, R) | 1.77, R |
| S-F71, S-F76 · unsigned location | 1/102, 1/105 (not R) | 6.8, 8.0, R |

**Verdicts.**
- **Positional F-test in human AKI: "partly".** Its genes replicate weakly (1.4×); its sets do not (1 of 102), because
  human location is not measurable here.
- **Best on level: the published programs (S-P80).** 0.59 against ≤ 0.25, stable in the healthy-only contrast (0.64).
  Among our own lists, 91's level screen does best (sets 1.55×, genes 1.38×), and clustering's sets do not replicate.

**Claim checks** (human z against mouse z at 12 h):

| Check | Human | Mouse |
|---|---|---|
| K1 injury markers up | 3.1 ✓ | 9.9 ✓ |
| K2 cell-cycle program up | 3.3 ✓ | 2.7 ✓ |
| K3 transport (T_PT) down | −3.6 ✓ | −11.1 ✓ |
| K4 T_PT loss larger in S3 | −0.4 | −3.5 ✓ |
| K7 male-biased PT down | −1.9 | −5.9 ✓ |
| K8 fasting-down down / fasting-up up | −0.2 / −4.8 | −7.8 ✓ / 1.0 |
| K9 induced up / lost down | 4.9 ✓ / −0.1 | 6.8 ✓ / −7.1 ✓ |
| K5, K6 (S3 location) | none hold | 3 hold |

- On identity-retaining nuclei only (post hoc), K1 and K9-induced still hold. K2 and K3 do not (1.7, −0.9): in human
  AKI, transport loss and cell-cycle gain come with the shift into injured states.
- **K10 fails.**
  - PROGENy EGFR (3.9) and PI3K (2.9) rise in human AKI, at q ≤ 0.01.
  - MAPK does not (0.5).
  - Human AKI's strongest footprints are JAK-STAT (6.9), TNFα (4.0) and NF-κB (3.0), with hypoxia and androgen falling.

**What this answers.**
1. **The core injury programs generalise to human AKI:** injury markers, the induced program, cell-cycle gain and
   transport loss. In human they come with the shift into adaptive and degenerative PT states, rather than with change
   inside identity-retaining cells.
2. **Positional (where-along-the-PT) calls cannot be tested in human biopsies.**
   - S3 is under-sampled in cortical cores.
   - Location signal is about zero, so neither the positional F-test's sets nor 91's location calls replicate.
   - This is a limit of the human data, not evidence against the mouse calls, which replicate in mouse (94).
3. **Mouse-specific readings stay mouse-specific:** the male-biased program loss, the fasting signature and the lost
   program do not appear in human AKI. The fasting-up program even falls.
4. **Method ranking for AKI level:**
   - published programs first;
   - then 91's level screen;
   - clustering's level sets do not generalise across species (1.06×), though they replicate in mouse.
5. **Signalling:** EGFR and PI3K activation is shared by mouse IRI and human AKI. MAPK is mouse-only here.
