# Is human PT axially flatter than mouse PT? Literature, first principles, and which stories survive

Workstream 3, phase 3. Written 2026-10-04.

## Protocol for the small quantitative checks (written before running them)

These use **public data** (Lake et al. 2023 snRNA v1.5; CELLxGENE mouse snRNA pseudobulk) and our
**already saved** contrasts (`results/pt_literature_deep/gene_level_lake_contrasts.csv`). No model is
fitted to our private data. Script: `SCRATCH/ws3/zonation_check.py`; outputs:
`results/pt_literature_deep/zonation_check/`.

**Q1 · Does tissue state flatten human PT zonation?** Lake nuclei from cortex (`region C`), all
disease groups (normal, AKI, CKD). Donors with ≥ 30 nuclei in cortical PT-S1 and in cortical PT-S3.
- Gene set G (fixed by mouse data only): genes in our ortholog universe with |male mouse snRNA
  S3 − S1| ≥ 1 log2 and pooled Lake cortical PT CPM ≥ 10.
- Donor amplitude k_d: least-squares slope through the origin of the donor's human (S3 − S1) on the
  mouse (S3 − S1) across G. k = 1 means mouse-like amplitude; k = 0 means no zonation of mouse-zonated genes.
- Injury burden: fraction of the donor's cortical PT nuclei labelled aPT, dPT or cycPT among
  PT-S1, PT-S2, PT-S3, aPT, dPT and cycPT (`dPT/DTL` excluded from both).
- Decision: **supports tissue-state flattening** if Spearman(k_d, injury burden) ≤ −0.3 with a donor
  bootstrap 95% interval below 0 **and** the AKI or CKD median k_d is below the normal median;
  otherwise **no support**. Caveat fixed in advance: injured nuclei leave the S1/S3 classes for
  aPT/dPT, so the remaining S1/S3 nuclei are the healthier ones and this test is conservative.

**Q2 · Is our donor's flattening unusual for human?** Our human-on-mouse slope is put on the Lake
scale by calibrating each of our species against the mouse snRNA reference on G:
k_ours = slope(our human S3 − S1 on mouse snRNA) ÷ slope(our mouse S3 − S1 on mouse snRNA).
Compare with the distribution of k_d over Lake **normal** donors. Descriptive only.

**Q3 · Which story genes go beyond flattening?** Flattening model: human contrast = k̄ × mouse
contrast, with k̄ the median k_d of the 7 primary normal Lake donors (cortex, ≥ 50 nuclei per
segment; notebook 38). For each R5 gene, observed = mean Lake human S3 − S1 over those 7 donors,
predicted = k̄ × mean of male and female mouse snRNA S3 − S1, SE = donor SD/√7.
- *Reversal beyond flattening*: observed has the opposite sign to the mouse contrast,
  |observed| ≥ 0.25 and at least 6 of 7 donors share its sign.
- *Human-specific gradient*: |mouse contrast| < 0.25 or the gene is absent in mouse, and
  |observed| ≥ 0.5 with at least 6 of 7 donors sharing the sign.
- *Exceeds flattening*: same sign as the mouse contrast, |observed| > |predicted| and
  |observed − predicted| > 2 SE (human more zoned than flattening predicts).
- *Flatter than predicted*: |observed| < |predicted| and |observed − predicted| > 2 SE, without
  meeting the reversal rule.
- *Consistent with flattening*: |observed − predicted| ≤ 2 SE and no reversal.
Rules are applied in the order listed above.
The same classification is repeated on our own segment means with k_ours, as a check.

---
## Bottom line

1. **Human PT really is less axially differentiated than mouse PT, but the deficit is smaller,
   and more specific to particular gene programs, than the mouse-selected slope implies.** In
   independent snRNA data (Lake 2023 cortex versus the mouse atlas), human S3 − S1 contrasts are
   about 0.12 of mouse contrasts on genes the mouse zonates. That number is lowered by a selection
   asymmetry: on genes the *human* zonates, mouse contrasts are 0.62 of human (0.75 after the
   within-species selection baseline). In distribution terms human has 3× fewer genes with
   |S3 − S1| ≥ 1 log2 (326 vs 989) and half the amplitude (100 most zonated genes: median 2.3 vs
   5.0 log2). Canonical transport-segment markers keep rodent-like zonation in human (SGLT2,
   SGLT1, AQP1, OATs, SLC5A8, SLC7A13). The mouse-only zonation sits in metabolic, xenobiotic and
   sterol programs, several of them androgen-linked in male mice.
2. **Tissue state modulates human zonation but cannot produce it.** Across 33 Lake donors, the
   amplitude falls with the injured-PT fraction (Spearman −0.42, 95% interval −0.66 to −0.10;
   pre-specified rule met). The effect is small: extrapolated to zero injury the amplitude is 0.147
   of mouse, against 0.12 in normal donors. Region is not the explanation either: mouse cortical
   medullary-ray PST (microdissected PTS2) already differs from S1 as much as outer-stripe PTS3
   does (1,273 and 1,235 genes ≥ 1 log2). Our donor (0.26 on a platform-calibrated scale) is not
   flatter than the Lake donors (0.08–0.16).
3. **Stories beyond flattening:** Dcxr, Ugt3a1 and Acox2 (reversals in Lake and in our data);
   Rbp4, Aox1 and Ugt1a9 (human-specific gradients in both). Gatm, Acaa2 and Gclm reverse in Lake,
   but **in our own data they match the flattening prediction**, so those claims rest on the atlas.
   **Consistent with flattening:** Acadm, Gclc, Slc6a19, the sterol genes, Crot/Nudt19, Lpl, Inmt
   and Cyp24a1. These should be presented as instances of the global difference.

## 1 · Is weaker axial differentiation of human PT reported?

**Direct statement: NOT FOUND.** Eight searches found no study stating that human PT segments are
less differentiated along the axis than rodent segments. The classic human and primate
ultrastructure papers exist: *Human renal ultrastructure. I. Proximal tubule of healthy
individuals* (1966, [PMID 5920363](https://pubmed.ncbi.nlm.nih.gov/5920363/)) and *Ultrastructure
of the proximal tubule of the rhesus monkey kidney* (1969,
[PMID 4980756](https://pubmed.ncbi.nlm.nih.gov/4980756/)). Neither has an abstract in Europe PMC
and neither is open access, so I cannot cite their content. A six-species study (including man,
rat and mouse) of where Bowman's capsule epithelium meets PT epithelium found that position varies
by species; human and rat are alike
([PMID 8376192](https://pubmed.ncbi.nlm.nih.gov/8376192/), abstract). It says nothing about
zonation along the PT.

**Counter-evidence: canonical transporters are zoned in human as in rodents.**
- SGLT2 is in S1/S2 and SGLT1 in S3 brush border; the renal locations "largely resembled those in
  rats and mice" ([PMID 25304002](https://pubmed.ncbi.nlm.nih.gov/25304002/), abstract).
- AQP1 protein and mRNA are highest in the human proximal straight tubule, weaker in the
  convoluted part ([PMID 9013443](https://pubmed.ncbi.nlm.nih.gov/9013443/), abstract).
- Human OAT1–3 stain more in S1/S2 than S3 (PMID 27053689, round-2 check).
- Lake cortex reproduces these: SLC5A2 −2.6, SLC5A1 +4.1, SLC5A8 +4.8, SLC7A13 +6.0, SLC22A7 +2.7
  log2 (S3 − S1), close to mouse.

The flat human genes are mostly not transporters. Examples: SLC5A12 (−0.2 vs mouse −8.8),
SLC34A1 (−0.6 vs −3.7), Cyp7b1 (−0.3 vs +7.3), Gatm, Acaa2, the sterol and FAO programs.

**Single-cell atlases.** Lake 2023 could not resolve PT-S1 from PT-S2 in its single-cell data
and kept them merged ([PMC10356613](https://pubmed.ncbi.nlm.nih.gov/37468583/), OA full text). In
its snRNA, PT-S2 is the smallest normal cortical class (1,744 nuclei against 5,902 PT-S1 and 3,826
PT-S3; from the dataset's obs). Healthy mouse snRNA resolves S1, S2 and S3 clusters
([PMID 32571916](https://pubmed.ncbi.nlm.nih.gov/32571916/), OA full text). I found no published
side-by-side measure of how distinct the subclasses are in each species. Our symmetric check above
is the first quantitative version we have.

**Architecture and hormones (first principles).**
- Mouse PT sex differences are concentrated in S2/S3
  ([PMID 36758122](https://pubmed.ncbi.nlm.nih.gov/36758122/), abstract). Healthy human kidney has
  few sex-biased genes, and only 9 are shared with mouse (PMID 19277126, round-2 check). So part
  of the mouse "late program" is androgen-induced late-PT expression that a male human PT does not
  mount. Female mice still show most of these programs (notebook 38), so androgen is part of the
  story, not all of it.
- Human nephron models still take axial transporter densities from rat. The 2019 human-nephron
  model uses rat values except for a few fitted adjustments
  ([PMC6405173](https://pubmed.ncbi.nlm.nih.gov/30802242/), OA). Human axial data are scarce enough
  that the field borrows rodent zonation; a positional human map fills that gap.
- I found no verified source comparing PT length or segment proportions between human and mouse.
  Any argument from "longer human PT, more gradual transitions" stays uncited.

## 2 · Could tissue state explain it?

**Literature: yes, partly, and the markers are known.**
- The three kinds of reference tissue differ in injury markers, cell-state proportions and age
  signatures ([PMID 42160485](https://pubmed.ncbi.nlm.nih.gov/42160485/), OA full text): tumour
  nephrectomy, living-donor pre-transplant biopsy and healthy-volunteer biopsy.
  - Early stress genes (ATF3, DUSP1, FKBP5, FOS, JUNB, JUN) are highest in tumour-nephrectomy and
    living-donor tissue.
  - HAVCR1 (KIM-1) is higher in living donors and IL18 in tumour nephrectomy; LCN2 is absent in all
    reference groups.
  - Adaptive/maladaptive PT (aPT) is marked by VCAM1, and age-associated PT genes are enriched in
    aPT.
- Reviews document procurement artefacts in human reference kidney tissue
  ([PMID 36811638](https://pubmed.ncbi.nlm.nih.gov/36811638/);
  [PMID 40062478](https://pubmed.ncbi.nlm.nih.gov/40062478/); abstracts).
- In mouse injury, failed-repair PT down-regulates terminal-differentiation markers (Slc5a12,
  Slc22a30, Slc7a13; the first is S1-restricted and the last rises toward S3 in mouse
  microdissection). Acute injured states form
  "injured S1/2" and "injured S3" clusters (PMID 32571916, OA). This is a mechanism by which injury
  could compress zonation.
- Cold ischaemia changes the kidney transcriptome by compartment in mouse Visium, mainly OXPHOS in
  the inner medulla ([PMID 42323650](https://pubmed.ncbi.nlm.nih.gov/42323650/), abstract). This is
  evidence of procurement effects, not of loss of zonation.
- No study measuring PT zonation as a function of ischaemia time, procurement route or donor age
  was found.

**Data (Q1, public Lake data, rule pre-specified).** Across 33 donors (13 normal, 14 CKD, 6 AKI),
zonation amplitude declines with the fraction of PT nuclei in injured states (aPT, dPT, cycPT):
Spearman −0.42, 95% interval −0.66 to −0.10. Medians: AKI 0.106, normal 0.120, CKD 0.125. The
decision rule says **supports**. The effect is modest: a linear fit gives 0.147 at zero injury
against 0.105 at 90% injury, while mouse is 1 by construction. Even "normal" Lake donors carry
24–87% of PT nuclei in injured states, so the rule describes reference tissue as it is collected.
**Conclusion: tissue state explains at most a small part of the flattening.**

**Markers to check in our human sections** (for workstream 2):
- early stress: ATF3, FOS, JUN, DUSP1;
- injury: HAVCR1, IL18, LCN2;
- failed repair: VCAM1;
- loss of terminal differentiation: SLC5A12, SLC7A13 (from the mouse work).

Compare with Lake normal donors. Our donor's amplitude (0.26, platform-calibrated) is already not
flatter than theirs.

## 3 · Which stories survive "beyond flattening", and how to frame it

Flattening model: human contrast = k̄ × mouse contrast, k̄ = 0.119 (median of the 7 primary Lake
donors). In our data k_ours = 0.256 on the snRNA-calibrated scale. Table:
`results/pt_literature_deep/zonation_check/story_genes_vs_flattening.csv`.

| Story | Gene | Lake class | Our data | Reading |
|---|---|---|---|---|
| Reversal | **Dcxr** | reversal (7/7 donors) | reversal | Beyond flattening in both: the best example |
| Reversal | Ugt3a1, Acox2 | reversal (6/7) | reversal | Beyond flattening in both |
| Lost in mouse | **Rbp4**, Aox1 | human-specific gradient (7/7) | human-specific | Beyond flattening: mouse does not express them |
| Detox | Ugt1a9 | listed as reversal* | human-specific | Human-specific (mouse microdissection ≤ 0.8 TPM). *The mouse snRNA signal is Ugt1a-cluster multimapping. |
| Creatine | **Gatm** | reversal (+0.51; 7/7) | consistent (−0.78 vs −0.81 predicted) | Beyond flattening **only in the atlas** |
| FAO | **Acaa2** | reversal (+0.58; 7/7) | flatter than predicted | Beyond flattening in the atlas; compatible in ours |
| FAO | Acadm | consistent | consistent | **Explained by flattening** |
| Detox | Gclm / Gss | reversal / exceeds | flatter | Atlas only |
| Detox | Gclc | consistent | consistent | **Explained by flattening** |
| Detox | Cyp2e1 | listed as reversal* | — | Presence difference: human CYP2E1 is about 3 CPM, so its gradient is noise |
| Early program | Slc6a19 | consistent | exceeds | Mostly flattening; Slc9a3 reverses in the atlas |
| Sterol, Lpl, Inmt, Crot/Nudt19, Cyp24a1 | | consistent or flatter | flatter | Explained by flattening (or presence differences) |

(Our-data classes use no standard error, so only the sign-based classes, reversal and
human-specific, are meaningful there.)

**How the Discussion should frame flattening, if workstream 2 confirms it.**
1. State it as a result. Human PT carries fewer and shallower S1→S3 gradients than mouse PT. This
   is replicated in an independent multi-donor human atlas and is not explained by cortex-only
   sampling or tissue state. The human program is largely a subset of the mouse one, and the
   canonical transport markers keep rodent-like zonation.
2. Change the null for species-by-position claims. "Human flat, mouse zoned" is the expected
   pattern, not a finding. A gene-level positional claim needs a deviation from human = k × mouse:
   - a reversal (Dcxr, Ugt3a1, Acox2);
   - a human-specific gradient (RBP4, AOX1, UGT1A9);
   - or human zonation where mouse has none.
3. Keep Gatm and Acaa2 as atlas-supported reversals, with the caveat that our own data cannot
   separate them from flattening. Move Acadm, Gclc, Slc6a19 and the sterol program into the
   flattening paragraph rather than presenting them as independent biology.
4. The R3 early/late landscape may be mostly flattening, as workstream 2 suspects. Mouse
   late-rising metabolic programs (FAO, sterol, glutathione, androgen-linked genes) are exactly the
   ones human does not zone. The Discussion should say this rather than read the landscape as
   positional biology.
5. Translational point: rodent models over-represent axial specialization of PT metabolism. Toxic
   and metabolic effects confined to a mouse segment may be spread along the whole human PT.
   Computational human nephron models currently borrow rodent zonation (PMC6405173).

## Limits of this check

- Both species' segment labels are snRNA cluster annotations. If the human clusters partition a
  continuum less cleanly, amplitudes shrink. That is part of what "flatter" means here and cannot
  be separated without human anatomical sampling.
- The mouse snRNA reference pools cortex and outer stripe. The region-matched microdissection
  check argues this does not matter.
- Lake donors are older and include tumour-nephrectomy and biopsy tissue. Our donor's source and
  age are unknown (`docs/paper/open_items.md`).
- The symmetric amplitude comparison and the region-matched microdissection check were added after
  the pre-specified Q1–Q3. They are exploratory.

## Files

- `SCRATCH/ws3/zonation_check.py` (script)
- `results/pt_literature_deep/zonation_check/`: `summary.json`, `lake_donor_amplitude.csv`,
  `story_genes_vs_flattening.csv` and the cached Lake pseudobulk `lake_cortex_s1_s3_all_disease.npz`
