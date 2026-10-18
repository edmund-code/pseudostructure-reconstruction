# Regulators and signalling footprints along the PT (notebook 95)

**The question.** Every pathway method so far scores a pathway by its *members*. Footprint methods score a regulator
or a signalling pathway by its *targets*: the genes it is known to switch on or off (PROGENy, Schubert 2018, PMID
29295995; decoupleR and CollecTRI, Badia-i-Mompel 2022, PMID 36699385; both read in full,
`docs/research/continuous-analysis/b-continuous-pathways.md` row 8). Two questions:
1. Does target-based activity give PT biology that is specific, replicated between animals and confirmed in
   independent data? The biology covers healthy zonation in each species, the species comparison and injury.
2. Does it beat member-based sets scored the same way?

We ask it on the same scorecard: our specimen controls, and the external benchmarks of notebook 94 (Kirita 2020 IRI)
and of the cell atlases.

Protocol frozen on 2026-10-18, before any regulator or footprint activity was computed on our data or on external
data. The network files had been downloaded and their structure tabulated: rows, sources, complexes, sign columns and
PROGENy pathways. Notebook 94's protocol and code existed, but it had produced no result yet.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/95_pt_regulator_activity.ipynb`; logic `95.regulators.1`; seed 95.
- Results `results/pt_regulator_activity/`.
- Module `pseudospace/regulator_activity.py`, tests `tests/test_regulator_activity.py`.
- **Networks** (OmniPath web service, 2026-10-18, academic licence; outside Git under
  `data/external/regulator_networks/`). Their sha256 is recorded by the notebook's guard G1:
  - `collectri_human.tsv`: `interactions?datasets=collectri&genesymbols=yes&organisms=9606&loops=yes`, 64,978 edges
    and 1,209 sources;
  - `progeny_human.tsv`: `annotations?resources=PROGENy`, 14 pathways with a weight and p for each gene.
- **Orthologs:** `results/minimal_pt_scfates/ortholog_map_used.csv` (one-to-one, as everywhere).
- **Our statistics** (read only):
  - notebook 90's end contrasts per specimen (`d_end_contrasts.csv.gz`);
  - notebook 91's gene table (L, the four pair L, L R1, L R2, Δ late, control expression);
  - notebook 71's AKI library sets.
- **External statistics:**
  - notebook 94's `external_genes.csv.gz` (Kirita level and location per time point, its t*) and its
    `pseudobulk_counts.csv.gz` (sham segment pseudobulks per mouse);
  - Census human PT segment pseudobulks (`data/external/census_pt_segments/`, donors, "proximal convoluted tubule"
    against "proximal tubule segment 3").

## 2 · Families of sources (each scored the same way)

| Family | Sources | Weights |
|---|---|---|
| REG | CollecTRI regulons | ±1. NF-κB and AP-1 dimers collapse to NFKB and AP1 (decoupleR's rule). An edge is −1 when it is inhibition-only and +1 otherwise. Other complexes are dropped. |
| FOOT | PROGENy, 14 pathways | the published weights of each pathway's 500 genes with the smallest p (decoupleR's default) |
| HALL | MSigDB Hallmark sets (71's library) | 1 for members. These are response signatures, so they sit between members and footprints. |
| MEM | KEGG and Reactome sets (71's library) | 1 for members: pathway membership |
| mRNA | the regulator's own gene (REG sources that map to one mouse gene) | the gene's own statistic, read as its z |

- Sources need ≥ 10 targets or members in the universe of the layer.
- Human symbols map to mouse through the one-to-one map. Unmapped targets are dropped.

## 3 · The estimator (one for every family)

- **ulm** (decoupleR's univariate linear model): for a statistic vector y over the layer's genes, the t value of the
  slope in y ~ 1 + w. Non-targets have w = 0.
- **Null:** 1,000 permutations of y within the deciles of control expression (91's log2 CP10k; deciles over the layer's
  genes). z = (t − mean) / SD of the permuted t, with a normal two-sided p and BH within family × vector.
- The permutation treats genes as independent, so z is a ranking device. Specificity comes from relabeled contrasts,
  replication and external data, as in notebook 91.

## 4 · Layers, statistic vectors and calls

**Z · Healthy zonation, per species** (the 7,971 genes of notebook 90's end contrasts; E = mean of the last 30% −
the first 30% of a specimen's own curve).
- Vectors: E for Ctrl1A2, Ctrl1A4, HUK1_COR1 and HUK1_MED1.
- **Positional source in a species:** q ≤ 0.05 in both of that species' specimens, with the same sign. Human
  replication is across two sections of one donor.
- **Species categories:**
  - conserved: positional in both species, with the same sign;
  - reversed: positional in both, with opposite signs;
  - mouse only / human only.
- **External:**
  - mouse: Kirita's sham mice, y = the mean over the three mice of log2 CPM(PTS3) − log2 CPM(PTS1) (94's
    pseudobulks, on 94's E ∩ the layer's genes);
  - human: Census human donors, y = the mean over donors of log2 CPM(segment 3) − log2 CPM(convoluted), for donors
    with both cell types and ≥ 50 cells in each.
  - Replicated means external q ≤ 0.05 with the same sign.

**I · Injury level** (91's 10,361 shared genes).
- Vectors: L (condition), the four IR × control pair L, L R1 and L R2.
- **Called:** q ≤ 0.05 on L; all four pairs have t of the same sign; both relabelings have |t_R| ≤ 0.5 |t|.
- **Specificity count:** the same q rule with R1 or R2 as the contrast.
- **External:** Kirita level t at 94's t* (and at every t). Replicated means external q ≤ 0.05 with the same sign.

**IL · Injury location (descriptive; 91 has no signed relabeled location).**
- Vector: Δ late − L.
- External: Kirita location t at t*.
- Reported as agreement only.

## 5 · Scorecard and verdicts

Per family and layer:

| Criterion | Layer | Better |
|---|---|---|
| between-mouse agreement: Spearman ρ of z (Ctrl1A2 against Ctrl1A4) over the family's sources | Z | higher |
| external precision enrichment, mouse: (replicated share of mouse-positional sources) / (replicated share of all the family's sources) | Z | higher |
| external precision enrichment, human (Census), for human-positional sources | Z | higher |
| specificity: mean relabeled count / condition count | I | lower |
| pair agreement: share of q ≤ 0.05 sources (condition) whose four pair t share their sign | I | higher |
| external precision enrichment at t*: (replicated share of called sources) / (replicated share of all) | I | higher |

- **Enrichment p:** hypergeometric over the family's sources.
- **Dominance** (notebook 93's `dominance`, tolerances 0.05 for ρ, shares and ratios, 0.25 for enrichments) gives
  "better", "worse", "equal" or "mixed" per layer.
- **Comparisons:** REG against MEM, REG against HALL, FOOT against HALL, FOOT against MEM, HALL against MEM.
- A family with fewer than 5 calls on a layer is "undefined" there. This applies to FOOT in most layers, since it has
  only 14 sources.
- FOOT is also compared pathway by pathway with its named member and response sets:

  | PROGENy | Hallmark | KEGG |
  |---|---|---|
  | TNFa | TNF-alpha Signaling via NF-kB | TNF signaling pathway |
  | NFkB | — | NF-kappa B signaling pathway |
  | p53 | p53 Pathway | p53 signaling pathway |
  | Hypoxia | Hypoxia | HIF-1 signaling pathway |
  | TGFb | TGF-beta Signaling | TGF-beta signaling pathway |
  | JAK-STAT | IL-6/JAK/STAT3 Signaling | JAK-STAT signaling pathway |
  | EGFR | — | ErbB signaling pathway |
  | MAPK | — | MAPK signaling pathway |
  | PI3K | PI3K/AKT/mTOR Signaling | PI3K-Akt signaling pathway |
  | VEGF | — | VEGF signaling pathway |
  | WNT | Wnt-beta Catenin Signaling | Wnt signaling pathway |
  | Estrogen | Estrogen Response Early | Estrogen signaling pathway |
  | Androgen | Androgen Response | — |
  | Trail | — | Apoptosis |

  For each row: the injury-level z of each, the pairs' signs, the relabeled z, and the external z at t*.
  Descriptive.
- **mRNA:** for REG calls, the share whose TF mRNA L has the activity's sign. This is descriptive, since TF activity
  is often post-transcriptional.

## 6 · Recovery checks (named now, from Kirita 2020's main text, read in full; a miss is reported)

| Check | Source | Expectation in our IR against control |
|---|---|---|
| RC1 | Hnf4a, Hnf1b and Pbx1 regulons drive differentiation genes in the successful-repair PT cluster | REG HNF4A, HNF1B, PBX1: injury level down |
| RC2 | NFkB and Relb regulons are specific to failed-repair PT | REG NFKB (and RELB if a source): level up |
| RC3 | Tcf7l1 regulon is specific to failed-repair PT | REG TCF7L1: level up |
| RC4 | Myc is expressed by the acute injured states (mRNA claim; the regulon is our extension) | REG MYC: level up |

- RC2 and RC3 concern a state that is rare before 2 days, so a miss may reflect our unknown time point.
- Each check is also read in Kirita's own data at every time point. Kirita's claims came from SCENIC on these nuclei,
  so the external reading of RC1–RC4 is a consistency check, not independent support.

## 7 · Outputs and figures

- **Tables:**
  - `sources.csv`: family, targets per layer;
  - `activity_zonation.csv`, `activity_injury.csv`: t, z, q per vector;
  - `calls.csv`;
  - `scorecard.csv`, `verdicts.csv`;
  - `foot_vs_sets.csv`, `recovery_checks.csv`, `mrna_agreement.csv`;
  - `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 95a: zonation, mouse z against human z for REG, with conserved and reversed sources marked;
  - 95b: the scorecard by family;
  - 95c: injury z against Kirita z at t*, per family;
  - 95d: PROGENy against its named sets, injury z with pairs and relabelings.

## 8 · Caveats

- The networks are mostly human, and PROGENy comes from acute perturbations of cell lines. Transfer to the adult
  mouse PT in situ is untested, and this notebook is that test.
- Regulons overlap pathways (HNF4A targets include many PT transporters). A regulator call can restate a pathway call
  under a regulator's name. Activity is not proof of a regulator's action.
- One human donor: "human positional" is replication across two sections.
- Injury location has no relabeled control, so it stays descriptive.
