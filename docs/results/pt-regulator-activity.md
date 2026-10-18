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

## Notebook 95 results

`analysis/notebooks/95_pt_regulator_activity.ipynb` (logic `95.regulators.1`) writes to `results/pt_regulator_activity/`.
It followed the protocol frozen in d64bfd1. Guards G1–G2 passed and t* is 12 h. Runtime is 25 s. Module
`pseudospace/regulator_activity.py` with `tests/test_regulator_activity.py` (3 tests).

**How the protocol was read.**
- **Correction:** G1's edge count. The frozen 64,978 was `wc -l` minus the header, but the file's last line has no
  newline, so it holds 64,979 rows. The 1,209 sources match.
- **Post hoc, not judged:** section 7.1, threshold-free agreement per family.
- Network sha256: CollecTRI `fb2dec0ffda1…`, PROGENy `a1a0f938211d…`.

**Sources.**

| Family | Z (7,971 genes) | I (10,361 genes) |
|---|---|---|
| REG (CollecTRI regulons) | 322 | 390 |
| FOOT (PROGENy) | 14 | 14 |
| HALL (Hallmark) | 50 | 50 |
| MEM (KEGG, Reactome) | 1,292 | 1,399 |

**Scorecard** (protocol section 5):

| | REG | FOOT | HALL | MEM |
|---|---|---|---|---|
| Z · between-mouse ρ | 0.929 | 0.657 | 0.925 | 0.909 |
| Z · mouse-positional sources (scorable), replicated in Kirita's sham | 27, 1 | 0, 0 | 3, 2 | 22, 8 |
| Z · mouse external enrichment | 11.9 (p 0.08) | – | – (< 5) | 58.1 (p < 0.001) |
| Z · human-positional sources | 1 | 2 | 3 | 3 |
| I · q ≤ 0.05 (condition; mean of relabelings) | 135; 83 | 8; 5 | 33; 22.5 | 514; 173.5 |
| I · specificity (relabeled / condition) | 0.615 | 0.625 | 0.682 | **0.338** |
| I · pair agreement | 1.00 | 1.00 | 1.00 | 1.00 |
| I · called; replicated in Kirita at 12 h | 35; 8 | 3; 3 | 8; 7 | 205; 133 |
| I · external enrichment | 1.17 (p 0.37) | – | 1.46 (p 0.09) | **1.87 (p < 0.001)** |

**Verdicts** (dominance):

| Comparison | Zonation | Injury |
|---|---|---|
| REG vs MEM | worse | worse |
| REG vs HALL | equal | mixed |
| FOOT vs HALL | worse | better* |
| FOOT vs MEM | worse | worse |
| HALL vs MEM | equal | worse |

\*FOOT against HALL on injury rests on specificity alone (0.625 against 0.682). FOOT's external enrichment is undefined
with 3 calls.

**Post hoc: threshold-free agreement** (Spearman ρ over each family's sources, not judged):

| | REG | FOOT | HALL | MEM |
|---|---|---|---|---|
| Z · Ctrl1A2 vs Ctrl1A4 | 0.929 | 0.657 | 0.925 | 0.909 |
| Z · COR1 vs MED1 (one donor) | 0.738 | 0.899 | 0.760 | 0.479 |
| Z · mouse vs Kirita sham | 0.638 | 0.868 | 0.586 | 0.576 |
| Z · human vs Census | 0.400 | 0.640 | 0.415 | 0.359 |
| Z · mouse vs human | 0.174 | −0.055 | 0.133 | 0.231 |
| I · L vs Kirita 12 h | 0.621 | 0.719 | 0.721 | 0.688 |
| I · R1 vs Kirita 12 h / 14 d | 0.380 / 0.491 | 0.367 / 0.662 | 0.227 / 0.651 | 0.228 / 0.592 |
| I · R2 vs Kirita 12 h / 14 d | −0.038 / −0.208 | 0.398 / −0.262 | 0.113 / −0.362 | 0.222 / −0.224 |

**What the families find.**
- **Zonation, mouse.** Few sources pass q ≤ 0.05 in both mice, although their rankings agree at ρ 0.91–0.93. The
  per-specimen end contrasts are heavy-tailed, so a mean-based ulm t is small for most sources.
  - MEM recovers notebook 90's positioned sets. S1-high: amino-acid, inorganic-ion and SLC transport, basigin, mineral
    absorption. S3-high: PPAR signalling, peroxisome, peroxisomal lipid metabolism, butanoate.
  - REG adds regulator names. S1-high: HNF1A, KLF15, VDR, NR1I2 and NR1I3 (PXR, CAR), SMAD3, MITF, TBX2. S3-high: PPARA,
    BHLHE40/41, MECP2, REST.
  - In human only NR1I3 (S1-high in both species) is positional, and in Kirita's sham mice only KLF15 replicates.
- **Injury level.**
  - MEM calls 205 sets, 133 of which replicate in Kirita.
  - REG calls 35 regulons, 28 of them up. Only 8 replicate (DMTF1, HMGA2, HR, POU4F1, TBP, TCF4, TCF7L2, ZBTB17),
    which is no better than all regulons.
  - HALL calls 8, of which 7 replicate: androgen response, G2-M, IL-6/JAK/STAT3, mitotic spindle, PI3K/AKT/mTOR, ROS
    and UV response up.
  - FOOT calls 3 and all 3 replicate: **EGFR, MAPK, PI3K**.
- **PROGENy against its named sets (injury z; Kirita z at 12 h).**
  - **Footprints see what member sets miss** for receptor-tyrosine-kinase signalling:
    - EGFR 10.6 (Kirita 7.8), against KEGG ErbB 1.6;
    - MAPK 8.4 (6.8), against KEGG MAPK 1.9.

    This is Schubert's rationale: these pathways act post-translationally, so their members barely change while their
    targets do.
  - **Footprints miss what response sets see:**
    - p53 −0.2, against Hallmark 3.9 and KEGG 3.9;
    - hypoxia −0.8, against Hallmark 3.2;
    - JAK-STAT −0.2, against Hallmark IL-6/JAK/STAT3 6.1;
    - androgen −2.0, against Hallmark 3.0;
    - Trail −0.1, against KEGG apoptosis 5.5.
  - TNFα, NF-κB and TGFβ are strong in all three.
- **Recovery checks (Kirita's regulons).**
  - RC1 holds for HNF1B (down, z −3.2). It fails for HNF4A (+1.2) and runs opposite for PBX1 (+4.1).
  - RC2 holds in direction for NFKB and RELB (z 8.9; Kirita z about 9 at every time). They fail the call rule,
    because relabeling 1 is large.
  - RC3: TCF7L1 has fewer than 10 targets, so it is not a source.
  - RC4 holds: MYC is up (z 7.0; Kirita z 5–6 at 4–12 h).
- **mRNA.** Of 24 REG calls with a mapped TF gene, 67% move in the activity's direction.

**What this answers.**
1. **Target-based activity is not a better way to call PT biology here.**
   - Regulons are less specific than member sets on injury (0.62 against 0.34 relabeled per condition hit).
   - They replicate externally no better than chance (1.17×, against 1.87× for member sets).
   - Their zonation calls rarely carry to human or to Kirita.
   - Under the frozen scorecard both footprint families are "worse" than KEGG and Reactome membership.
2. **Footprints add one specific, replicated piece of biology:** EGFR, MAPK and PI3K signalling rise in the injured
   PT. Each is called (pairs agree, relabelings small) and replicated in Kirita (z 5.3–7.8), and member sets miss
   EGFR and MAPK. Use PROGENy beside the sets for signalling, not instead of them.
3. **Regulator names are hypotheses, not calls.** The positional regulons are HNF1A, KLF15, VDR and PXR/CAR toward
   S1, and PPARA toward S3. They sit on the same transporter and peroxisome genes as the member sets, so they restate
   the sets' biology under a regulator's name.
4. **Post hoc: relabeling 1 resembles injury itself.** IR2A2 + Ctrl1A2 against IR2A4 + Ctrl1A4 correlates with
   Kirita's IRI response at ρ 0.23–0.38 at 12 h and 0.49–0.66 at 14 d, in every family. Relabeling 2 correlates near 0
   or negatively.
   - The A2 sections share an injury- or repair-like program in both conditions, as the cohort memory noted for
     stress genes.
   - So the AKI calibration failures (0.120, 0.144) partly count real, shared biology as false calls. This fits notebook
     94: calls that fail internal calibration replicate externally.
