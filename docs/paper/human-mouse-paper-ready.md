# The human-vs-mouse axis: what is paper-ready after the package switch

Written 2026-10-20. The species-axis companion to `paper-ready-results.md`, which ranks both axes together. This
page goes deeper on one axis because that is where the work is now.

Sources: `docs/results/pt-package-issues.md` (what the switch broke), `pt-package-species-calls.md` (77),
`pt-package-continuous-vs-discrete.md` (74), `pt-coexpression-conservation.md` (97), `pt-non-one-to-one.md` (86),
`pt-positional-decomposition.md` part II (91), `pt-gene-atlas.md` and `pt-evidence-model.md` (112, 113),
`pt-published-findings.md` (81, 82), `pt-disease-genetics.md` (108), `pt-transporter-handoffs.md` (101).

## 0 · The one structural conclusion

**The species axis is a gene-level result. The package switch is what made that true, and the paper should follow
it rather than resist it.** Three independent findings force it:

- **Atlas replication of the 83 species calls is not set-specific.** They replicate against random sets and against
  noise calls, but not against random genes of the same positional strength (p 0.35 and 0.64; on the retired
  coordinate these were 0.0013 and 0.0057). Replication supports **genes**, not gene sets (77).
- **There is no pathway-level species program.** Of 381 sets, none is reversed, conserved or species-dominant at
  q ≤ 0.05; the set cosines simply spread around the genome-wide 0.23 (z mean −0.20, SD 0.93) (91 S2).
- **On this axis the positional F-test is the most informative screen we have** — coefficient 0.88 of the nine
  species lines, above zonation class 0.64, atlas level 0.54 and handoffs 0.40. That is the exact opposite of its
  AKI ranking (0.05, interval crossing zero) (113).

So: genes lead, the F-test is the primary screen, and the pathway list is description rather than evidence.

## 1 · The strongest species results, ranked

### S1 · The reversed gradients, corroborated three independent ways — the headline
About a quarter of the gradient is shared between species, and **74 genes reverse direction** (90 D, 91 S). Three
separate methods, none of which can see the others' statistic, agree on them:

1. **Co-expression rewiring, which needs no levels at all.** Median rewiring D is 0.761 for reversed genes against
   0.213 for conserved (Mann–Whitney p < 10⁻¹⁰, on the 31 reversed genes inside 97's 3,679-gene network set), and
   the rewiring itself **replicates in independent atlases** (ρ 0.30; 0.34 among genes whose neighbourhoods
   reproduce in both species) (97).
2. **Convergence across the nine species screens.** The k ≥ 5 species core is 23 genes after the contamination
   filters, and it *is* a reversal list: Slc16a9, Lap3, Mme, Hpgd, Slc17a3, Slc5a12, Slc27a2, Vnn1, Acot12.
   Convergence on this axis is sharp — enrichment over the expression-stratified null 1.31× at k ≥ 3, 2.65× at 4,
   **7.36× at 5, 58× at 6** — and validated by annotation the screens never saw (anchored share 0.154 background
   → 0.870) (112).
3. **The two lists overlap gene by gene.** Of 97's 25 "rewired in both species" genes, Mme, Lap3, Hpgd, Slc17a3,
   Slc16a9, Acot12 and Ugt3a1 are on the independent reversal list.

And it is **guarded by a negative**: reversed and conserved gradients are equally constrained in sequence
(AUC 0.49–0.51, CMH OR 0.78) and an artefact forest cannot separate them (0.54), so the list is not an artefact of
relaxed constraint on fast-evolving genes (106).

Three methods, two of them level-free, one externally replicated, plus a sequence control. Nothing else on this
axis has that much behind it.

### S2 · Conservation rises with gradient strength — **still on the retired coordinate**
The old central claim (index ≈ 0.55 in the top 1% of genes, 0.1–0.3 among moderately zonated genes that still
reproduce within species, near zero below; direction agreement 0.5 → 0.85–0.89 with strength) comes from notebook
49, which reads `results/minimal_pt_scfates/`. **It has not been re-run on the package units.** The human S2
relabelling (the package calls 47–60% of the old human S2 units S1, and S2−S1 classified genes fall from 9,013 to
6,269) acts directly on the gradient this index is built from, so the number can move.

This is the single blocking item on the species axis. Everything in S1 and S3–S6 is already on the package.

### S3 · The conserved core is transport, and it sits where disease sits
- **16 genes conserved in both species' co-expression and in gradient**: Slc4a4, Slc6a19, Slc7a8, Dpep1, Gpx3,
  Apoe, Iyd, Lgmn, Tpp1, Cotl1, Sptan1, Itpr2, Fam151a, Tbx2, Atp6v0c, Gnptg — core PT transport and handling (97).
- Transporters are over-represented among conserved genes (33 of 43 zonated transporters; a hand count, still open
  item 18).
- **11 conserved transporter handoffs**, confirmed in microdissection and in the atlases, overlapping and inside
  S2: amino acids, organic anions, SGLT2→SGLT5. Handoffs are *not* a family-wide principle — that question fails —
  which is what makes the 11 worth naming (101).
- **Mendelian proximal-tubulopathy genes sit toward S1 in both species** (mean E −1.47 mouse / −0.56 human, beyond
  matched nulls), 25 of 27 literature positions are reproduced, and OCT2 rises toward S3 in mouse but is flat in
  human (108).

### S4 · The only external cross-species positional arbiter
Cheval et al. 2012 (PMID 23056504) microdissected S1 and S3 in **both** human and mouse. Our species × position
interaction replicates in their data (82, B2b). This is the one external check the species axis has, it is
program-level, and it is what licenses any set-level reading at all. It also qualifies the novelty claim: they
report global divergence but no per-gene or per-program profile comparison.

### S5 · Non-1:1 orthology is where the mouse's zonation hides — a self-contained second section
A 1:1 analysis structurally cannot see any of this (86):
- **Strong-zonation share by orthology class:** 1:1 3.8%, 1:many 9.9%, many:many 14.2%, but many:1 only **1.7%**.
  The excess sits in mouse genes of **mouse-expanded** families, not in mouse single-copy genes whose human side
  expanded. MH OR 1.84–4.05 under every sensitivity. 41 strong probe-clean non-1:1 mouse genes, 25 of them outside
  the shared space (Cyp2d12, Nat8f5, Slc22a19, Scd1, Ugt2b38, Serpina1a, the Cyp4a cluster, Slco1a1/1a6 …).
- **Sex-biased mouse PT genes are enriched in non-1:1 families**: any bias OR 1.50 (p 0.001), **male bias 1.90
  (p 10⁻⁴)**, female 1.16 (p 0.40). This is an open question in the gene-age literature, never tested in kidney,
  and the sex bias comes from external female-vs-male data, so our all-male cohort cannot confound it.
- **Tandem paralogs share positions, and it is specific to paralogs**: median r 0.31 for tandem pairs under 1 Mb
  against 0.07 for dispersed paralogs and 0.04 for non-paralog neighbours within 1 Mb (human 0.48 / 0.07 / 0.04);
  similarity falls with distance; and tandem status predicts shared position while sequence identity does not
  (t 3.7 against 0.5). Co-regulation of tandem duplicates beats a plain neighbourhood effect.
- Counterweights in the same notebook: the human counterparts of strongly zonated mouse non-1:1 genes are zonated
  about as often as 1:1 partners (20% against 22%); family sums fail calibration (0.089 against the 0.075 limit),
  so no family-level call is claimable; paralogs are measured as reliably as 1:1 genes (65% against 70% agreement
  with microdissection).

### S6 · Where the species differ, the old answer was wrong
36 human-high and 47 mouse-high calls; **most differ most in S3, in both directions**, and human-high calls split
between S1 and S3. The old "human-high in S1, mouse-high in late S2" pattern does not survive the switch. 89% of
resolved calls still summarise mouse zonation the human sections do not share (77).

## 2 · What the switch broke, and must be stated in the paper

1. **"Continuous beats discrete" is weak at set level.** The segment F-test passes the swap rule with 76 species
   calls, 71 of them shared with the positional F-test; the spline adds 12 redundancy clusters (10 near misses) and
   the steps add 5. **The supported message is narrower and should replace the old one:** testing *shape against
   specimen shape noise* finds species differences that level-based clustering and conventional DE either miss or
   confuse with probe-panel differences (71, 74, 75).
2. **Within-segment claims are forbidden by the frozen noise gate.** Human gene-fold agreement 0.49, equal depth
   0.73 (human S1 0.55), species registration gap 0.37, **10 of 12 cells noise-dominated**, and one fold refit
   reverses the mouse PT. The transporter micro-gradients (mouse 115, human 62 genes) and the within-segment
   gradient panels are therefore descriptive, hatched, and not claims (72, 74).
3. **Our landmark coordinates beat the package within segments** (gene-fold 0.93–0.97, 0 of 12 noise-dominated),
   but the elastic curve is frozen as the only primary coordinate. Their role — comparator, dropped, or a
   within-segment add-on — is still undecided and is a user decision.
4. **Calibration is narrow:** the species F-test passes at 0.071 against the 0.075 limit (71).
5. **The human labels moved.** Human PT-S2 shrinks to 632 units from 1,813; HUK1_MED1's segment transitions cross
   on the curve; human units are effectively undefined for within-species statements.
6. **Human zonation is the weakest measurement in the study.** Against microdissected human S1/S3, our
   cortex-section human S3 − S1 agrees at **ρ 0.14**, against 0.46 for mouse. Only two mouse gradients are
   confirmed at program level and no human position is (81, 82).

## 3 · Three confounds specific to this axis

1. **The probe panels differ, so levels are never readable.** Probe exposure of the calls: F-test 0.17, clustering
   0.69 (71, 104). Genes with no included probe are structural zeros — mouse Gsta1, Nat8f6/7, Apol9b and human
   UGT2B15, AKR1C1 — so a "human-only alpha-GST" result would be a panel artefact (86).
2. **Neighbour signal is asymmetric between species, which is worse than neighbour signal being present.** Mouse S1
   units carry DCT, podocyte and Ren1 signal at −2.0 to −5.0 centred, where our human S1 units carry at most −1.2.
   A mouse-only S1 bias can therefore be neighbour RNA rather than biology; S100g is the worked example, with a
   mouse S1 signal of −5.0 in our units and −0.07 in Census mouse PT. Species calls are 1.4–1.9× enriched for
   edge-enriched genes, the F-test least at 1.44× (81 issue 24, 104). The human S3 end has its own flagged set:
   SLC12A1, UMOD, EGF, KCNJ1, AQP2 (109).
3. **The human tissue is one donor, two cortex sections.** There is no human counterpart to mouse outer-stripe S3
   at all, and the two "replicates" are sections of one block, not donors.

## 4 · Recommended Results structure for the species paper

| # | Section | Rests on |
|---|---|---|
| R1 | Cohort, segmented structures, the reconstructed axis and its measured limits | 70, 72 (state §2.2 here, not in the discussion) |
| R2 | A gene-level species comparison: what shape-against-noise testing finds that level screens and conventional DE do not | 71, 74, 75 (the narrowed claim of §2.1) |
| R3 | **Reversed gradients, corroborated by a level-free method, by convergence across screens, and in external atlases; not explained by sequence constraint** | 90 D, 91 S, 97, 112, 106 |
| R4 | The conserved core: transport, the 11 handoffs, and disease genes toward S1 in both species | 97, 101, 108 |
| R5 | Non-1:1 orthology: where mouse zonation hides, and the sex program inside it | 86 |
| R6 | Conservation rises with gradient strength | **49, after the re-run** |
| R7 | Negatives: no pathway-level species program; no sequence prior; co-expression not useful genome-wide; paralog substitution fails all three tests | 91 S2, 106, 97, 88 |

R3 should lead the biology because it is the only claim with three independent methods and an external
replication behind it. R6 cannot be written until notebook 49 is re-run.

## 5 · What to package, for the species use case

The priority differs from the AKI one in `paper-ready-results.md`:

1. **The positional F-test with its relabeling calibration** — still first, and here it is the *best* screen rather
   than the highest-volume one (113: 0.88). But for cross-species use it must ship **with the probe-exposure and
   neighbour-marker filters wired in as a required pre-step**, not an option, because that is where this axis
   actually fails. `pseudospace/nb_gam.py`, `pathway_calibration.py`, `contamination.py`,
   `pseudospace/species_calls.py`.
2. **The level-free comparison layer:** gradient divergence plus co-expression rewiring
   (`pseudospace/coexpression_conservation.py`). Across platforms levels are unusable, and this is the only
   cross-species measure that needs none. Notebook 97's "not shown useful" verdict is about its genome-wide rule
   (median neighbourhood reproducibility 0.67 against a 0.7 threshold); as a **corroboration layer on an existing
   gene list** it worked and replicated externally. Ship it with that scope stated.
3. **Orthology handling:** the 1:1 space, family sums, and the non-1:1 class comparisons
   (`pseudospace/non_one_to_one.py`, `ortholog_pairing.py`), with the family-sum calibration failure wired in as a
   refusal rather than a footnote.

**Do not ship** set-level species calls, or any cross-species level comparison.

## 6 · The blocking item

Re-run notebook 49 (conservation by gradient strength) on the package units and labels. It is the last SCF13-era
number in the species story, it is the paper's old central claim, and the human relabelling acts directly on it.
