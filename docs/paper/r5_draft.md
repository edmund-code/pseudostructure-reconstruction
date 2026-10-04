# R5 Results and Discussion — draft (workstream 3, phase 2)

Source of every number: `analysis/notebooks/38_pt_biology_external_depth.ipynb` (logic version
`38.literature_deep.2`), outputs in `results/pt_literature_deep/`. Gene-level evidence:
`story_gene_evidence.csv`; story verdicts: `story_verdicts.csv`; genome-wide check:
`genome_wide_lake_check.csv`; Figures 11–12 in `results/pt_literature_deep/figures/`.
Citations marked [V] were verified in phase 2 (abstract or open-access full text; details,
including which papers are not open access, in `SCRATCH/ws3/literature_phase2.md`); [R2] were
verified in the round-2 literature check (`docs/results/pt-pathway-literature-novelty.md`).

**Conditional on workstream 1.** Nothing below states or depends on the R3 early/late landscape
(human-high pathways peaking early, mouse-high late). That pattern could be produced by
attenuated human positions and is under test in workstream 1. The R5 claims are tested on
**reviewed segment labels**, not on the coordinate, and against an independent human atlas whose
labels we did not make; the segment-label attenuation estimate is reported in R5.0. Where a
sentence would link an R5 story to the landscape it is marked ⟦if WS1 clears attenuation⟧.

---

## Results

### R5.0 An independent human reference with true S1/S2/S3 labels in cortex

No human study has microdissected PT segments, and notebook 34's human reference compared
convoluted PT with S3. Here the human comparator is the KPMP/HuBMAP single-nucleus atlas
(Lake et al. 2023 [V], CELLxGENE dataset `a12ccb9b`): nuclei from normal tissue with author labels
PT-S1, PT-S2 or PT-S3, restricted to cortex (`region C`), and the 7 donors with at least 50 nuclei
in each segment (4 female, 3 male); 12 donors at 20 nuclei served as a sensitivity set. The mouse
comparators are the snRNA atlas by sex (12 male and 12 female donors) and microdissected PT
segments from male (Chen et al. 2021 [V]) and, new here, female mice (Chen et al. 2023 [V], GEO
GSE212213). Every claim was written down, with the rules for keeping or dropping it, before these
data were opened (notebook 38, section 0).

*Genome-wide*, the cortex-only human reference reproduces our human zonation (Spearman 0.43 over
6,787 genes; 93.1% direction agreement for 335 strongly zonated genes) and our species-by-position
interaction: for the 10% most position-dependent genes, Spearman 0.70 with male mouse snRNA as the
mouse arm and 96.7% direction agreement for 276 genes strong in both (female mice: 0.67 and 93.7%;
female microdissection: 0.64 and 91.5%).

*Attenuation.* If our human segment means were blurred, every human contrast would shrink by a
common factor relative to an external reference. Among 70 control genes chosen only from the
external data (strongly zonated in the same direction in both species), our human contrasts scale
to the human atlas as our mouse contrasts scale to the mouse atlas: relative human attenuation
*a* = 0.86 (gene-bootstrap 95% interval 0.64–1.35). Every "human flat" claim below also survives
division of our human contrast by *a*. This is a segment-label check; whether the continuous
coordinate attenuates human positions is workstream 1's question.

*Probe panels.* Our two species were measured with different Visium HD probe panels. Of 11,071
tested genes, 82.1% carry the same number of probes in both panels; 610 have one probe in one
species and two or more in the other (197 with a single mouse probe). A gene with unequal probes
cannot be a headline finding in this section (table `probe_counts.csv`).

Of 39 claimed genes tested (three probe-panel examples are reported separately), 30 met every
pre-specified test (*headline*), 5 were *supporting* and 4 did *not replicate* (Figure 12). Two of
the 30 (OAT1, NaDC3) pass only against outer-stripe mouse S3 and reverse against the cortical
mouse segment (R5.6). The four lead stories are below.

### R5.1 Creatine synthesis: the first step is confined to early PT in mouse and spans the whole human PT

GATM (AGAT) makes guanidinoacetate, the committed step of creatine synthesis, and the kidney is the
main site of that step [R2: PMID 17928413]. GATM is a human PT enzyme whose mutation causes renal
Fanconi syndrome [R2: PMID 29654216]. In our data mouse Gatm falls steeply from S1 to S3 while human
Gatm stays high along the whole PT (Figure 11). Every reference agrees: the human atlas shows no
fall in cortex (S3 − S1 = +0.51 log2; +0.57 in men, +0.47 in women), whereas mouse falls by
6.8 (male snRNA), 6.3 (female snRNA), 12.1 (male microdissection) and 11.0 (female microdissection)
log2 units; rat falls by 3.9. Gatm is not a region artefact: even the cortical medullary-ray
segment of mouse (microdissected PTS2) has 23.5 TPM against 4,270 TPM in PTS1. Probes are balanced
(3/3), mouse Gatm is if anything female-biased so male mice understate it, and GATM protein is
high in human tubules (Human Protein Atlas, "enhanced" reliability).

The rodent half is not new: transamidinase activity in microdissected rat nephron is confined to
S1 and S2, with S1 higher [V: PMID 1378964], and our mouse transcript profile replicates that in a
second rodent. What is new is the human half and therefore the contrast: no study has resolved
GATM along the human PT, and human kidney does release guanidinoacetate in vivo, as rat kidney
does [V: PMID 17928413], so the species differ in where along the tubule the step runs, not in
whether the kidney runs it.

The second step, GAMT, is classically hepatic, with renal guanidinoacetate exported to liver
[V: PMID 17928413]. Within the mouse PT, however, Gamt rises from S1 to later segments in
microdissection (male 31 → 85 → 82 TPM; female 26 → 96 → 61), and in rat, so the two steps are
spatially separated along the mouse tubule; no segment-resolved GAMT study exists in any species
[V: 4 targeted searches]. **This part rests on external data only:** our Visium data do not
resolve Gamt (flat in both species; 30% mouse detection; unequal probes 3/4). Human GAMT is flat
in the atlas.

### R5.2 Mitochondrial β-oxidation: two steps of one pathway are zoned in opposite directions in mouse and flat in human

In mouse, Acadm (MCAD) rises from S1 to S3 while Acaa2 (mitochondrial 3-ketoacyl-CoA thiolase) is
confined to S1; in human both are nearly flat (Figure 11). Both genes passed every test: human
atlas S3 − S1 = +0.46 (ACADM, about a fifth of the mouse rise and in the same direction) and +0.58
(ACAA2, opposite to the mouse fall), against mouse +2.1 / +2.3 (Acadm, male/female snRNA) and
−4.0 / −4.2 (Acaa2); microdissection agrees in both sexes (Acadm 318 → 3,387 → 3,569 TPM; Acaa2
1,371 → 66 → 2 TPM in males). Both hold against the cortical mouse segment (PTS2), neither is
sex-biased in mouse PT, probes are balanced, and ACADM and ACAA2 proteins are high in human
tubules (HPA), so "flat" in human is not "absent". Peroxisomal and fatty-acid uptake members of
the program behave the same way (Nudt19, Slc27a2 headline; Crot supporting because human
detection in our data is 8%, although the atlas detects CROT and HPA shows high PT protein). Two
qualifications: Acsm3's late rise is a **male-mouse** program (female microdissection 0.0 / 1.3 /
0.3 TPM), consistent with the known sex dimorphism of mouse PT fatty-acid metabolism [R2: PMID
34651140; PMID 36758122]; and Hadh did not replicate (both species fall alike in the external data).
A 2026 intravital study shows that late mouse PT uses lipase-rich lysosomes to mobilize lipid
droplets while early PT catabolizes filtered protein [V: PMID 41803103], an independent,
mechanistic reason to expect late mouse PT to be specialized for lipid handling.

Classical rat microchemistry measured β-oxidation *capacity* along the nephron and found
3-hydroxyacyl-CoA dehydrogenase activity roughly uniform across cortical proximal and distal
segments [V: PMID 7200500; PMID 2260682], with the β-oxidation/citrate-synthase ratio highest in
mid-proximal segments [V: PMID 6588129]. Our mouse Hadh is likewise nearly flat (microdissection
490 / 601 / 375 TPM), so an activity assay of that step would not see the gene-specific pattern
reported here: two other steps of the same pathway, Acadm and Acaa2, are zoned in opposite
directions in mouse. No segment-resolved data for ACADM, ACAA2 or HADH, and no human
segment-resolved fatty-acid-oxidation data, were found [V]. The difference from rat enzyme
activity (different species, enzyme activity versus mRNA) is a stated limitation.

### R5.3 Where along the tubule drug-handling capacity sits differs between species

Human PT expresses the glucuronidation enzyme UGT1A9 most strongly in early PT, declining toward
S3, and mouse PT does not express Ugt1a9 at all (microdissection ≤ 0.8 TPM in every segment, both
sexes, with three probes in each panel): human atlas S3 − S1 = −0.67. UGT1A9 is among the most
abundant renal UGTs in human kidney [V: PMID 25650382]; the only human kidney UGT
immunohistochemistry used a pan-UGT1A antibody and described a uniform tubular distribution
without separating convoluted from straight PT [V: PMID 17698974], so an isoform-specific,
early-high UGT1A9 profile is unreported. Mouse PT instead places glutathione synthesis late: Gclc,
Gclm and Gss rise from S1 to S3 in both sexes and in microdissection (female Gclc 37 → 199 → 490
TPM), while human expression is flat (GCLC +0.61 against mouse +3.3/+2.9). The axial glutathione
literature is species-dependent: in isolated perfused rabbit PT, GSH *synthesis* is highest in S1
while *cellular concentration* is highest in S3 [V: PMID 9612330], and in rat glutathione is
highest in convoluted and early straight PT [V: PMID 998800]; no segment-resolved GCLC or GCLM data
exist in any species [V]. One mouse injury-model scRNA-seq study reports Gclc as "exclusively
expressed in the distal nephron" [V: PMID 38313210]; microdissected healthy mouse PT contradicts
that for the PT (male 41 / 871 / 445 TPM). Mouse Cyp2e1 is restricted to S1–S2 in our data and in
microdissection; mouse renal CYP2E1 is a known cortical, proximal-tubular, male and
testosterone-regulated enzyme [V: PMID 2260985; PMID 7839370], and it is lower but present in female
mice (2.6 / 93.7 / 0.2 TPM). In human, CYP2E1 mRNA is essentially absent from our donor and from
the atlas (peak 3 CPM); the literature conflicts (no p-nitrophenol oxidase activity in six human
kidneys [V: PMID 7839370] versus CYP2E1 protein reported in human PT cells [V: PMID 38042273]),
so we state low or absent mRNA, not absent protein. EPHX1, a known human proximal and distal
tubular enzyme [V: PMID 7769232], did **not** replicate positionally (the two human references
disagree in sign and HPA detects no PT protein) and is dropped. Human AOX1 rises toward late PT
(atlas +1.13); it is a known human tubular protein [V: PMID 11510964], and humans carry one AOX
gene where mouse and rat express four tissue-restricted isoforms [V: PMID 26920149]. No aldehyde
oxidase gene (Aox1, Aox3, Aox4) is expressed in mouse or rat PT microdissection; that absence is
our reanalysis, not a literature statement.

*Verdict note.* Under the pre-specified rule this story is **supporting, not lead**: its lead gene
GCLC carries 6 human probes against 3 mouse probes. Every biological test on GCLC passed, the
probe-balanced members GCLM and GSS are headline, and an excess of human probes would raise human
sensitivity rather than create a mouse-specific rise; the rule is applied as written and the
decision to promote it is left to the coordinator.

The translational reading is about *position*: a rodent model reproduces neither where human PT
conjugates (early) nor where mouse PT synthesizes glutathione and activates acetaminophen (S1–S2
and late), so S3-targeted injury in rodents need not predict the human segment at risk.

### R5.4 Genes whose species difference changes sign along the PT

A single whole-PT fold change cannot represent genes whose species difference reverses along the
axis (notebooks 31–32: 70 such genes, 24 with a null whole-PT comparison). The clearest is Dcxr
(dicarbonyl/L-xylulose reductase): it falls from S1 to S3 in mouse and rises in human. All five
external references agree (human atlas +1.53, +1.59 in men and +1.48 in women; mouse −0.46/−0.54
snRNA, −1.94/−2.03 microdissection; rat −1.45), detection is high in both species, probes are
balanced, the mouse fall is not to zero, and the reversal holds against the cortical mouse segment.
DCXR protein is known in murine proximal tubule [V: PMID 11882650] and DCXR loss tracks human CKD
[V: PMID 31217356]; its axial distribution is unreported in either species.
All eight reversal genes met the pre-specified rules, with two qualifications: Cyp24a1's human
profile differs between sexes in the atlas (+1.44 in men, −1.03 in women), and IGFBP4's human rise
is an mRNA result with no detectable tubular protein (HPA), so neither is named in the text.

### R5.5 Supporting programs

- **Mouse has lost a late-PT program that human and rat keep (RBP4).** Human RBP4 rises steeply
  toward late PT (atlas +3.46), rat RBP4 rises toward S3 (+5.4 log2; rat in situ hybridization
  placed RBP mRNA in outer-stripe S3 [R2: PMID 2469758]), and mouse Rbp4 is 0.0 TPM in every
  microdissected segment of both sexes with three probes per panel. The rat study itself warns
  that kidney RBP mRNA (S3) and immunoreactive RBP (convoluted PT, reabsorbed ligand) sit in
  different places and that perinephric fat expresses RBP mRNA [V: PMID 2469758]; our claim is
  therefore an mRNA claim, and ambient signal from perirenal fat in a cortex section is untested. AOX1 belongs
  here too but is supporting (2 human versus 3 mouse probes).
- **Human keeps early-PT transporters in late PT** (Slc6a19, Slc9a3, Slc4a4). All three met every
  test, including against the cortical mouse segment (SLC6A19: human atlas −0.30 against mouse −5.3;
  SLC9A3: +0.72 against −4.8). Phase 1 treated this as one hypothesis rather than three findings,
  and it remains one: human cortical S3 retains an S1/S2-type transport program that mouse loses.
  Mouse B0AT1's S1–S2 restriction is known [R2: PMID 35979966]; the human profile is not reported.
- **Mouse late-PT sterol synthesis** (Hmgcs1, Hmgcr headline; Cyp51, Sqle supporting) holds in
  female mice and in the cortical PTS2 segment, which weakens the phase-1 region objection, but
  rat runs the opposite way, so it is a mouse, not a rodent, program; it stays supporting.
- **Mouse-only late lipid genes** (Lpl, Inmt; Pank1, Coasy headline). Me1 did not replicate as
  mouse-only: the human atlas detects ME1 (20 CPM), falling toward S3; the mouse late rise has a rat
  precedent (malic enzyme highest in proximal straight segments [V: PMID 6588129]). Inmt is
  sex-biased in mouse PT [V: PMID 40259083].

### R5.6 What did not hold, and why it matters for cross-species maps

- **Serine synthesis is not human-specific.** Our mouse data barely detect Phgdh (1% of
  structures) and Psat1 (6%), but each has a single mouse probe against three human probes, and
  mouse microdissection shows both S1-restricted in both sexes (Phgdh 89 / 0.4 / 0.7 TPM in males).
  Human early-high zonation replicates (atlas −1.93 and −2.19). The correct statement is
  "S1-high in both species"; the apparent species difference is a probe-panel artefact. Gstp1 (one
  mouse probe; 234–425 TPM in mouse microdissection) is dropped for the same reason.
- **The organic-anion module flips with anatomy.** OAT1 and NaDC3 pass against outer-stripe mouse
  S3 but reverse sign against the cortical medullary-ray segment, where mouse OAT1 peaks (970 TPM);
  human OAT1–3 are already known to be S1/S2 > S3 [R2: PMID 27053689]. It belongs in R2/R3 as the
  clearest case of a whole-PT null hiding opposite regional differences, not in R5 as biology.
- Gamt (in our data), Hadh, Ephx1 and Me1 did not replicate (above); they are the four
  *not replicated* rows of Figure 12, shown alongside the passes.
- A pre-specified check of human outer-medulla PT could not be run: no atlas donor had at least 30
  nuclei in cortical S1, cortical S3 and medullary S3.

---

## Discussion (draft sections)

**What is new.** To our knowledge no study has compared human and mouse PT gene programs as a
function of position along the tubule. Cross-species kidney work to date compares cell types,
whole-PT expression or disease states: a unified mouse–human atlas found little overlap of
single-cell differential genes but agreement at pathway level [V: PMID 37639336], and a 2025
cross-species atlas defines conserved cell states across human and rodent models [V: PMID
40775269, abstract]. Our contribution is the positional comparison and a short list of programs
in which even the gene-level contrast is unreported: early-confined mouse versus PT-wide human
creatine synthesis, oppositely zoned mitochondrial β-oxidation steps in mouse, the placement of
conjugation, glutathione synthesis and CYP2E1, and genes such as DCXR whose species difference
changes sign. Each was checked against a cortex-only multi-donor human atlas, both mouse sexes
and direct microdissection, with the rules fixed in advance, and the list includes what failed.
The human atlas compares with mouse only through cell-state correlation heatmaps [V: Lake 2023,
OA full text]; the 2025 mouse sex atlas compares whole-PT sex differences with human snRNA, not
position [V: PMID 40259083, OA full text]; regional human kidney profiling has no mouse arm
[V: PMID 38513647, abstract]. The 2025 cross-species atlas is not open access and its PT
annotation depth could not be checked; this is the one residual novelty risk. Classic rat
microdissection already placed transamidinase in S1–S2 [V: PMID 1378964] and β-oxidation capacity
along the nephron [V: PMID 7200500], so the novelty is the human arm and the gene-level contrast,
not rodent zonation as such.

**Placement, not presence.** The recurring pattern is that human and mouse PT express the same
enzymes and transporters but place them differently: several mouse genes are sharply restricted
to S1 (Gatm, Acaa2, Slc6a19) or to S2–S3 (Acadm, glutathione and sterol synthesis, Cyp2e1), while
the human counterparts are expressed along the whole PT. A whole-PT comparison sees these as
modest magnitude differences or misses them entirely. ⟦if WS1 clears attenuation: this is
consistent with the R3 observation that species differences concentrate at the ends of the
axis.⟧ Whether human PT is genuinely less sharply partitioned, or whether human S1/S2/S3 are
less distinct as cell states, is open; the attenuation check argues against a simple labelling
artefact in our data, but the human atlas labels are themselves cluster annotations, and the same
atlas could not separate PT-S1 from PT-S2 in its single-cell data [V: Lake 2023]. We found no
published account of how human S1/S2/S3 map onto the rodent convoluted / medullary-ray / outer-
stripe scheme [V: 4 searches], which is why the human arm is best described as a position along a
continuous coordinate rather than as transferred rodent segment labels. Within-segment structure
is not unique to our map: the 2025 mouse spatial atlas describes a molecular subdivision of S3
[V: PMID 40259083].

**Relation to existing atlases.** Human atlases label PT-S1/S2/S3 as discrete classes of
dissociated nuclei [V: Lake 2023], and mouse microdissection gives three anatomical points with no
human counterpart [V: Chen 2021; Chen 2023]. Our map adds a continuous axis measured on intact
tubules in tissue in both species. We used the newest independent atlases as replication rather
than competition, and their agreement with our interaction (96.7% for the most position-dependent
genes) is the main evidence that the map is not specific to our specimens.

**Precedents for species differences in PT handling.** Species differ in renal transporter and
drug-metabolizing enzyme composition at the protein level [R2: PMID 31123036; PMID 38711199], and
human and mouse rely on different phosphate transporters [V: PMID 42494142]. What has not been
shown is that such differences are organized along the tubule. For drug development, the
position matters because PT injury is segment-specific: a rodent's S1–S2 Cyp2e1 and late
glutathione synthesis do not mirror human early conjugation.

**Limitations.** Two male mice and two cortex sections from one male human donor; specimens, not
structures, are the replication units, and no statement here is a population estimate for
humans. Human sections contain no outer stripe, so mouse outer-stripe S3 has no human
counterpart; we therefore also compared against the mouse cortical medullary-ray segment, which
changes the conclusion for the OAT module. The two species were measured with different probe
panels; genes with unequal probes are flagged and cannot be headline results, and one apparent
species difference (serine synthesis) turned out to be a panel artefact. Mouse PT sex
differences concentrate in S2/S3; female mice confirm most late mouse programs but not Acsm3.
External human positional data are snRNA cluster annotations, not anatomy. Everything is
mRNA-level and descriptive; protein evidence is presence or absence (HPA), not zonation.

**What would move this forward.** Female mice and female human donors measured on the same
platform, matched cortical sampling in both species, and in situ validation of the
position-defining genes (GATM, ACAA2, ACADM, DCXR, RBP4, AOX1, UGT1A9) on human sections where S1
and S3 can be distinguished.

---

## Figure legends (draft)

**Figure 11 | R5 story genes against direct segment measurements.** Left, fitted human (orange)
and mouse (blue) curves along the PT coordinate (notebook 12 model) with each specimen's binned
means; shaded bands mark the reviewed S1/S2/S3 ranges. Right, each external dataset's S1, S2 and
S3 values centred on its own S1 (log2): human cortex snRNA (Lake et al. 2023, 7 donors), mouse
snRNA by sex, mouse microdissection by sex (Chen et al. 2021, 2023) and rat microdissection (Lee
et al. 2015). Panel titles give the pre-specified verdict and probes per gene (M/H).

**Figure 12 | Replication matrix.** One row per claimed gene, one column per pre-specified test
(filled ✓ pass, open × fail, grey not testable, blank not claimed), with verdict, probe counts and
post-hoc annotations (absent in female mouse; human sexes or human references disagree; sign flips
against the cortical mouse segment).

---

## Appendix · Proposed edits to committed files (not made; for the coordinator)

`docs/results/pt-pathway-literature-novelty.md`
1. **Misquote, Tier 1 glutathione row:** it says "rabbit PT GSH synthesis is highest in S3 (PMID
   9612330)". The abstract reports synthesis highest in **S1** and cellular GSH concentration highest
   in S3; rat glutathione is highest in convoluted and early straight PT (PMID 998800). Replace.
2. **Creatine row:** add that rodent early-PT restriction is KNOWN (rat transamidinase activity only
   in S1 and S2, S1 > S2; PMID 1378964); the species contrast and the human profile remain NOT FOUND.
3. **Tier 2 serine synthesis:** demote (single mouse probe for Phgdh and Psat1; mouse
   microdissection S1-high in both sexes; human early-high replicates). Remove Gstp1 from the human
   early-high row (1 vs 3 probes; mouse expresses it at 234–425 TPM).
4. **OAT module:** move from Tier 1 to the method illustration; note that the sign reverses against
   the cortical medullary-ray mouse segment (PTS2).
5. Add a probe-ratio column to the Tier 1 and Tier 2 tables (`results/pt_literature_deep/probe_counts.csv`),
   and cite PMID 38313210 (Gclc "distal-only" in injury scRNA) with the microdissection counterpoint.
6. **CYP2E1:** state the human conflict (PMID 7839370 vs PMID 38042273) and word our result as
   low/absent human mRNA.

`docs/paper/outline.md` R5 row: replace "serine synthesis" with the drug-handling placement story
and note that the detox story is *supporting* under the pre-specified rule (GCLC probe imbalance)
unless the coordinator decides otherwise.
