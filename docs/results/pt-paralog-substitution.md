# Species differences carried by other paralogs (notebook 88)

Frozen on 2026-10-13, before any statistic of notebook 88 was computed. The user asked, after notebook 86: "what are
questions that can only ask [with] non 1:1 orthologs. ex. are there differences in human mouse but are accounted for
by multiple genes in the other species? these type of questions. explore a bit more".

**What I had seen before freezing.**
- Notebook 86's results (`pt-non-one-to-one.md`), including H7b's table (71's calls against family-unit calls), H8's four
  families and their r values, and the W curves' zonation calls.
- Population counts for this notebook (no statistic): of 698 export pairs inside multi-member orthogroups, 307 have
  both genes PT-expressed on U_PT and 106 of those are 71 species calls; 33 of the 106 have another PT-expressed
  probe-clean family member in at least one species (25 in mouse, 18 in human); 20 zonated genes have ≥ 2 PT-expressed
  probe-clean family members in the other species; 1 family has ≥ 2 zonated members in both species.

## Questions only genes outside 1:1 orthology can answer

The 1:1 analysis pairs one human gene with one mouse gene. When a family has several members, a human–mouse difference
for the kept pair can have three readings that only the rest of the family can separate:
1. **Substitution (P1):** another member of the family in the other species carries the profile that the kept partner
   lacks. The function moved to a paralog, not away.
2. **Compensation (P2):** no single member matches, but the family sums agree: the difference is in how the family
   divides the work, not in what the family does.
3. **Division of labour (P3):** a single gene in one species spans what several genes split in the other.
Two further questions are descriptive here: whether position or sequence pairs the members of many:many families
(P4), and which species differences only the family units can see (P5). P6 lists, for each strongly zonated gene
without an ortholog, the other species' genes that occupy the same place.

## Labels

- Notebook: `analysis/notebooks/88_pt_paralog_substitution.ipynb`
- Logic: `88.paralog_substitution.1`
- Results: `results/pt_paralog_substitution/`
- Reuses `pseudospace/non_one_to_one.py` (no new module code is planned).

## Inputs (read only)

- Notebook 86: its W runs on U_PT for both species (stage caches `w_mouse_pt_shared`, `w_human_pt_shared`: centred
  curves, typical shapes, zonation), `w_calls_{mouse,human}.csv`, `h7_family_units.csv` and the family definitions
  (rebuilt with 86's rule from 62's inventories: orthogroups at support ≥ 3 that are not 1:1, family-clean panel
  members).
- Notebook 71: its gene statistics (condition q) and the export's kept pairs (`human_symbol`).
- Ensembl 98 human–mouse homologies (P4).
- Stage-01 native counts, the export's library and 71's knots and grid (for the null sums, through 86's `w_run` rule
  with 86's fixed shapes).

**Guards.** 86's guards G1–G4 are rerun; 86's cached W payloads must list exactly the genes of its `w_calls` tables;
refitting 86's mouse U_PT run on 200 genes must reproduce its centred curves (atol 1e-8).

## Definitions

- **Curve:** the two specimens' mean centred curve (86's W4, log2), on the grid points that all four specimens cover.
- **r(a, b):** Pearson correlation of two curves; **D(a, b):** the RMS difference of the two mean-centred curves (log2).
- **Strata** (P1, P3 nulls): within a species, quintiles of log2 CP10k and of Z among PT-expressed probe-clean genes
  (86's H4 strata, extended to all PT-expressed genes). **Expression deciles** (P2 null): deciles of log2(CP10k + 0.01)
  over all panel genes of the species on U_PT.
- **Random replacement:** a gene of the same species and stratum (or decile), outside the orthogroup, probe-clean, not
  mt, drawn independently for every replaced member (seed 88).

## P1 · Substitution

- **Pairs:** 71's species calls (condition q ≤ 0.05) among export pairs (m, h) in multi-member orthogroups (either gene's
  class 1:many, many:1 or many:many), both genes PT-expressed and probe-clean, with ≥ 1 other PT-expressed probe-clean
  member of the orthogroup in at least one species.
- **Gain** in a direction with alternatives: mouse side, r_alt = max over other mouse members m′ of r(h, m′) and
  gain = r_alt − r(h, m); human side likewise with r(m, h′). A pair's gain is its larger direction.
- **Substituted:** gain ≥ 0.5 and r_alt ≥ 0.5.
- **Null:** every alternative replaced at random (1,000 draws), same counts per direction. p = (1 + #draws with at least
  the observed number of substituted pairs) / 1,001, one-sided. Per pair, the empirical p of its gain is reported.

## P2 · Compensation

- **Pairs:** P1's 71-called pairs (both PT-expressed and probe-clean), whose family holds more than the kept gene in at
  least one species (86's family-clean members).
- **Sums:** per species, the summed counts of the family's family-clean panel members, fitted with 86's W model and 86's
  typical shapes (so a sum is read against the genes' typical shape).
- **Statistic:** L = log2(D(sum_h, sum_m) / D(h, m)). Accounted for: L ≤ −1 (the difference at least halves).
- **Null:** pseudo-families: the kept gene plus, for each other member, a random gene of that member's expression decile,
  summed and fitted the same way (200 draws per family).
- **Test:** each family's percentile of L among its null; the percentiles against one half, Wilcoxon, one-sided
  (lower). Beside it: the observed number accounted for against the null mean.

## P3 · Division of labour

- **Anchors:** zonated probe-clean genes S in a non-1:1 orthogroup whose other species holds ≥ 2 PT-expressed
  family-clean members E.
- **Statistic:** δ = r(S, sum_E) − max r(S, Eᵢ); the sum is fitted as in P2. "Sum beats best": δ > 0.
- **Null:** the E members replaced at random (strata), summed and fitted (200 draws). p = (1 + #draws with at least the
  observed number of anchors where the sum beats the best member) / 201, one-sided.

## P4–P6 (descriptive)

- **P4:** families with ≥ 2 PT-expressed probe-clean members in both species and ≥ 1 zonated member: the member pairing
  that maximises the summed r (Hungarian) against the pairing that maximises summed Ensembl identity; the share of
  positional pairs that are sequence pairs, against the share expected by chance (1/k).
- **P5:** 86's family units called (q joint ≤ 0.05) with relabeled p > 0.05 in both relabelings, whose family has no
  71-called kept pair; with each species' largest member Z and the family-sum curves. Family calls are not claimable
  one by one (issue 31).
- **P6:** for each strongly zonated probe-clean gene without an ortholog (no_ortholog or no_consensus), the five
  PT-expressed genes of the other species with the highest r.

## Multiplicity and reading

P1, P2 and P3 are the primary tests, reported raw and after BH over the three. A positive P1 or P2 means that some
human–mouse differences of the 1:1 analysis are differences in which paralog does the work. A negative result means
that, on this curve, the 1:1 differences stand when the rest of the family is counted.

## Bounds

As 86: two male mice, one male human donor, human positions the weakest measurement, probe panels per species. 71's
species calls are themselves calibrated only at the set level. The tests run across genes and families.

## Figures (plots only)

- 88a: P1: observed gains against the null; substituted pairs (human gene, kept mouse partner, substituting member).
- 88b: P2: observed L against the null percentiles; examples where family sums agree and the kept pair does not.
- 88c: P3: δ against the null; examples.

## Notebook 88 results

`analysis/notebooks/88_pt_paralog_substitution.ipynb` (logic `88.paralog_substitution.1`) writes to
`results/pt_paralog_substitution/`. It followed the protocol frozen in 09184cc. Guards G1–G5 passed: 86's cached W runs
list its PT-expressed genes, and 200 refitted mouse genes reproduce their curves. 87 grid points (curve position
0.05–0.82) are covered by all four specimens.

**Populations.** 93 of 71's species calls are kept pairs in multi-member orthogroups with both genes PT-expressed and
probe-clean. The pre-freeze count of 106 had not applied the probe rule. 32 of the 93 have an alternative member in
at least one species (24 in mouse, 18 in human).

### Primary tests: none of the three holds

| Test | Observed | Matched null | p | BH |
|---|---|---|---|---|
| P1 · substitution | 6 of 32 pairs | mean 7.5 (95th percentile 11) | 0.85 | 0.89 |
| P2 · compensation | median percentile 0.57; 2 of 57 families with the difference halved | null mean 8.6 halved | 0.89 | 0.89 |
| P3 · division of labour | sum beats best member for 2 of 22 anchors | mean 2.6 | 0.72 | 0.89 |

**Reading.** On this curve, the human–mouse differences of the 1:1 analysis stand when the rest of the family is
counted. Another paralog matches the partner's profile no more often than a random gene of the same expression and
zonation does. Family sums close the difference no more than random additions do. A single-copy gene is not matched
better by its co-orthologs' sum than by their best member.

### Details

- **P1, the six "substituted" pairs** (descriptive; per-pair empirical p in brackets):
  - human ABCB1 is tracked by mouse Abcb1b (r 0.67), not the kept Abcb1a (0.01) (p 0.015);
  - human AKR1C3 by Akr1c21 (0.73), not Akr1c18 (−0.58) (p 0.048);
  - human ALDH3B1 by Aldh3b3 (0.95), not Aldh3b1 (−0.92) (p 0.055);
  - mouse Zfp30 by human ZFP82 (p 0.08), Dhfr by DHFR2 (p 0.20), Tubb2a by TUBB2B (p 0.32).

  Abcb1b and Akr1c21 come closest to a substitution, but the count as a whole is what random genes give.
- **P2 · why sums rarely close the difference.** Adding random genes to the kept gene pulls a family sum toward the
  typical gene's shape and shrinks the human–mouse difference (null median L −0.03 overall, down to −2.5 for the 7-member
  Akr1c family). Real paralogs shrink it no more (median L −0.02). The 2 families whose difference halves (Akr1c18/AKR1C3
  L −2.27; Clec2h/CLEC2D −1.32) halve it less than their pseudo-families do (−2.54, −1.59). This fits 86's finding that
  paralogs co-position: summing a family adds more of the same profile.
- **P3:** the two anchors whose sum beats its best member are Ugt2b34 (δ 0.02) and Ifitm2 (0.005). The human UGT2B7 is
  matched by mouse Ugt2a3 (r 0.84) better than by the Ugt2b sum (0.44). The MT family (mouse 2, human 6 PT-expressed)
  has negative r in every comparison: the human MT genes do not share the mouse Mt1/Mt2 profile, summed or not.
- **P4 · position does not follow sequence** (descriptive): in 6 families with ≥ 2 members on both sides, the
  positional pairing equals the sequence pairing in 1 of 12 pairs, against 4.4 by chance. Examples: Gnai3 pairs
  positionally with GNAI1 (0.74), and Ugt2b34 with UGT2A3 (0.86).
- **P5 · differences only the family units see:** 19 family units are called (relabeled p > 0.05 in both) with no
  71-called kept pair. In most, the mouse family is zonated and the human one is not: Serpina1 (mouse Z 5.3, human
  0.3), Ceacam (4.4 against 0.6), Isoc2, Ifitm. In the Oat-like Slc22 family both sides are zonated (Slc22a19 Z 5.0;
  SLC22A24 2.6). In the other direction, human CYP2B6 is zonated (Z 3.4) while no mouse Cyp2b member is PT-expressed.
  Family calls are not claimable one by one (86, issue 31).
- **P6 is uninformative.** The curves are smooth (6-df splines), so on 87 grid points random human–mouse gene pairs
  already reach |r| 0.83 (5th–95th percentile −0.84 to 0.83). The top match among about 10,000 genes always reaches
  r ≈ 0.99, so "the gene in the other species at the same place" is not resolvable by profile correlation on this
  curve.

### Recorded readings

1. The P1/P2 population applies the frozen probe rule (93 pairs, not the pre-count's 106).
2. P6's lack of resolution was found after the run, from the spread of random-pair r; it is reported, not tested.
3. Figure 88a's titles are on two lines.

### Tables and figures

- Tables: `p1_substitution.csv`, `p1_null_counts.csv`, `p2_compensation.csv`, `p3_division_of_labour.csv`,
  `p4_positional_vs_sequence.csv`, `p5_family_only_differences.csv`, `p6_niches.csv`, `primary_tests.csv`.
- Figures: 88a substitution, 88b compensation, 88c division of labour.

### Relation to SAMap (added 2026-10-14, after the run)

P1 is the "paralog substitution" analysis of SAMap (Tarashansky et al. 2021, *eLife*, PMID 33944782; full text read).
SAMap calls a substitution when a paralog's expression correlation beats the ortholog's by more than 0.3. It has no null
and does not treat specimens as replicates. P1 moves the question onto the PT curve and adds the matched random-gene
null. Without that null, P6's smooth-curve correlations would have produced substitutions by chance. Cite SAMap if P1
is reported. `non-one-to-one/samap_assessment.md` says why running SAMap itself would add no testable result here.
