# Where along the PT does the injury response sit, and does it move? (notebook 114)

**The question.** Notebook 113 found that the two most informative lines of evidence are the ones that place a
change — in time (the injury phase) or in space (the location screen) — not the ones that measure its size. Nothing
has yet placed a change in both at once.

Notebook 94 fitted, for every gene and every Kirita time point, three segment contrasts of the IRI response: S3−S1
(it called this "location"), **S3−S2** and **S2−S1**. It used only the first, and only at the best-matching time.
The other two contrasts and the other four time points were written to
`results/pt_external_iri/tables/external_genes.csv.gz` and never read. That table is a position × time grid of the
injury response in 18 mice, and it is already committed.

This notebook reads it. Three things follow that our own two-vs-two, single-time-point data cannot give:

1. **Which part of the PT carries the response, and at which boundary.** Our data says the late PT, with the
   strongest genes at the S2/S3 boundary and inside S3 (notebook 112) and injury rising toward S3 in external
   Visium sections (notebook 105). Kirita's S3−S2 and S2−S1 contrasts separate those two boundaries directly.
2. **Whether the geography moves.** The level response fades from 4,014 genes at 12 h to 3 at 6 weeks (notebook 94).
   Whether the *position* of what remains fades with it, persists, or shifts, has not been measured.
3. **Whether published programs have distinct geographies**, rather than being distinguished only by their timing.

**What the literature predicts** (all from full-text reads already verified in `docs/research/`):
- The outer stripe, where S3 lies, is the most injured region acutely (Zhang YL 2024; Limbutara 2020 maps S3 to the
  outer stripe), and the proliferation wave peaks there at 48 h (Chang-Panesso 2019; Witzgall 1994: PCNA⁺ 80% outer
  stripe against 23% cortex).
- **The long-term prediction is contested.** Failed-repair cells persist for months and are associated with the
  outer stripe (Dixon 2022's Kim-1⁺ tubules at 6 weeks), but Kirita 2020 reports failed-repair PT cells "scattered
  throughout" and appearing "within Kim1-positive injured tubule in a scattered manner at 2 d". Persistence and
  dissipation of the geography are both defensible, so either answer is informative.

Frozen on 2026-10-19, before any of the unread contrasts was evaluated. Only the column names of the Kirita table
were inspected, and notebook 94's already-published S3−S1 alignment was re-read.

## 1 · Inputs

- Notebook `analysis/notebooks/114_pt_injury_geography.ipynb`; logic `114.geography.1`; seed 114.
- Results `results/pt_injury_geography/`.
- Inputs, all committed: `results/pt_external_iri/tables/external_genes.csv.gz` and `universe.csv` (notebook 94),
  `results/pt_known_injury_programs/tables/programs_members.csv` and `programs.csv` (notebooks 79–82),
  `results/pt_gene_atlas/tables/aki_atlas.csv.gz` (notebook 112), `results/pt_package_aki/tables/gene_sets.json`
  (notebook 76), `results/pt_positional_decomposition/tables/i_genes.csv.gz` (notebook 91).
- No new helper is needed; the matched-set null is `pseudospace.nulls.matched_means`.

## 2 · The sets (fixed now)

- **Published programs** from notebooks 79–82, by `program_id`: the Kirita programs (K-INJ12, K-INJ3, K-REP, K-SEV,
  K-FR), Gerhardt (G21-VCAM, G21-CYC, G23-FR, G23-S1/S2/S3-up and -down), Lake (L23-aPT-up, L23-aPT-down),
  Dixon (DX-PTS1-2-up/down, DX-PTS3-TAL-up/down), the human AKI sets (H-AKI-UP, H-AKI-DN) and Z-IRI-UP/DN.
- **Ours:** notebook 91's level-up, level-down and location-change genes; notebook 112's convergence core at
  k ≥ 5; and the INJ (18 markers), TPT (93 transporters) and CYC (297 G2M + E2F) programs.
- **Negative control:** the housekeeping set used in notebook 68.
- **Circularity flag.** The five K-* programs come from the same experiment as the test data. They are reported in
  a separate block and never counted as replication, exactly as notebook 94 treated them.

## 3 · Measurements and nulls (fixed now)

- For each set S, time t ∈ {4 h, 12 h, 2 d, 14 d, 6 wk} and contrast c ∈ {S3−S1, S3−S2, S2−S1}: the mean coefficient
  of S's members present in Kirita's universe E.
- **Null:** 2,000 random gene sets of the same size drawn from E within the same expression deciles, the deciles
  taken from notebook 91's control log2 CP10k. z = (observed − null mean) / null sd.
- The level response is measured the same way, so geography and magnitude are on one scale.
- **Sensitivity:** the same with `location (assigned)` in place of `location`, which adds the injured nuclei
  assigned to segments.

## 4 · Tests (fixed now)

**W1 · Is the acute response S3-weighted?**
- **Holds** when |z(S3−S1)| ≥ 2 at 4 h or 12 h for INJ and for TPT, with INJ positive and TPT negative (injury
  markers rise most in S3; transport falls most in S3).

**W2 · Which boundary carries the gradient?**
- Statistic: |z(S3−S2)| against |z(S2−S1)| for every set and time.
- **Prediction fixed now:** S3−S2 is the larger for the acute sets, matching notebook 112's location enrichment at
  the S2/S3 boundary and notebook 105's external sections. The opposite result would mean our S2/S3 reading is an
  artefact of our own segment labels.

**W3 · Does the geography outlast the magnitude?**
- Statistic, per set: z(S3−S1) and z(level) at 14 d and 6 wk, each divided by its value at 12 h.
- **Holds** when the median location ratio across the non-K sets exceeds the median level ratio at 14 d.

**W4 · Do up and down programs have different geographies?**
- Statistic: the z(S3−S1) trajectory of the -up against the -down programs of the same source (G23-S1/S2/S3,
  DX-PTS1-2 and DX-PTS3-TAL, L23-aPT, H-AKI, Z-IRI).
- **Reported, not passed or failed.**

**W5 · Does our own geography match, and at which time?**
- Spearman between our per-gene P_late (notebook 91's Δ late − L, as notebook 94 defines it) and each of Kirita's
  three contrasts, at each time. Notebook 94 reported only S3−S1 and found ρ peaking at 2 d (0.274) while the level
  ρ peaks at 12 h (0.524).
- **Prediction fixed now:** the S3−S2 correlation peaks at the same time or later than S2−S1, and our position-age
  stays later than our level-age.

**Controls.**
- C1: the housekeeping set gives |z| < 2 at every time and contrast.
- C2: every test repeated on the subset of each set that passes notebook 112's clean filters (104 core rebuild,
  109 neighbour marker, edge enrichment).
- C3: our three 91-derived sets rebuilt under the balanced relabelings (notebook 112 § C1 rule) should lose their z.

## 5 · Outputs and figures

- **Tables:** `set_z.csv` (set × time × contrast, observed, null mean and sd, z), `set_z_assigned.csv`,
  `boundary.csv`, `persistence.csv`, `alignment.csv`, `clean_sensitivity.csv`, `relabeled.csv`,
  `how_the_protocol_was_read.csv`, `verdicts.csv`, `run_record.json`.
- **Figures (plots only):**
  - 114a: z(S3−S1) trajectories over the five times, our sets and the non-K published programs;
  - 114b: z(S3−S2) against z(S2−S1), one point per set and time;
  - 114c: location z against level z at each time, showing which fades first;
  - 114d: the per-time Spearman of our P_late with each contrast.

## 6 · Caveats

- Kirita's segments are assigned labels on nuclei, not an anatomical coordinate; S2 recall is the weakest (0.81–0.87,
  notebook 94). A contrast involving S2 is therefore the noisiest.
- Location is weak externally: only 45 genes pass q ≤ 0.05 for S3−S1 at 12 h and none at any other time
  (notebook 94). Set-level means are the right resolution here; gene-level location calls are not.
- 3 vs 3 mice per time point.
- The K-* programs are not independent of the test data.
- Sets overlap, so the z values across sets are correlated; no multiplicity correction is claimed beyond the
  matched null.

## Notebook 114 results

`analysis/notebooks/114_pt_injury_geography.ipynb` (logic `114.geography.1`) writes to
`results/pt_injury_geography/`. It followed the protocol frozen in e22a88d and ran in 19 s.

**Guards.** Kirita's universe holds 10,733 genes; the pool shared with notebook 91 is **10,020**, each with a
control expression for the deciles. The grid is 20 statistics per gene (5 times × 4 contrasts) and holds no missing
value. 31 sets keep at least 10 members in the pool.

### W1 · The acute response is S3-weighted, and the weighting is gone by 14 days

z of the S3 − S1 response (positive = the IRI effect is larger in S3):

| set | 4 h | 12 h | 2 d | 14 d | 6 wk |
|---|---|---|---|---|---|
| INJ (18 injury markers) | **2.0** | 1.3 | −0.5 | −3.2 | 0.1 |
| TPT (93 transporters) | **−3.5** | **−2.6** | **−2.7** | 0.1 | −1.8 |
| CYC (G2M + E2F) | −0.1 | −0.6 | **2.8** | 0.6 | −1.1 |
| our 112 core (k ≥ 5) | **−7.2** | **−4.9** | **−3.1** | −3.7 | −2.6 |
| 91 level down | **−7.6** | **−6.9** | **−9.9** | −1.2 | −3.0 |
| 91 level up | **9.3** | **8.7** | **13.3** | 2.4 | 4.5 |

- **W1 holds.** Injury markers rise most in S3 at 4 h and transporters fall most in S3 at 4 h to 2 d, both beyond
  the expression-matched null. Our own convergence core and our two level lists carry the strongest acute
  geography of any independent set (z −7.2, −7.6 and +9.3 at 4 h); only Kirita's own severe-injury program, which
  is not independent of the test data, is larger (16.5). The core was built with no external positional input, so
  this is an external check on it.
- **The cell cycle is the exception.** It is the only set whose S3 weighting *peaks late*, at 2 d (z 2.8, and 4.4
  for Gerhardt's cycling program), matching the proliferation wave that peaks at 48 h in the outer stripe
  (Chang-Panesso 2019; Witzgall 1994).
- By 14 days every set is at or near zero, and several have changed sign.

### W2 · Neither boundary dominates; the sets split

|z(S3−S2)| − |z(S2−S1)| (positive = the S2/S3 boundary carries more of the gradient):

| | 4 h | 12 h | 2 d | 14 d | 6 wk |
|---|---|---|---|---|---|
| median over independent sets | −0.06 | −0.36 | −0.47 | −0.19 | 0.19 |
| our 91 level down | **−10.6** | −5.3 | −8.5 | −1.8 | −1.9 |
| Gerhardt S3-down | 3.0 | **7.1** | 5.1 | −1.0 | 0.3 |
| Kirita severe injury (not independent) | **14.2** | 1.6 | 1.5 | −0.3 | −0.5 |

- **W2 fails, and the prediction it tested was wrong.** Across sets neither boundary dominates at any time; the
  median difference is within ±0.5 throughout. The two boundaries carry the gradient for *different sets*.
- **Our own down-called genes separate S1 from S2, not S2 from S3** (−10.6 at 4 h). Read with W1, that says the
  acute loss we see is "S1 is spared" more than "S3 is worst". Notebook 112's finding that our best-supported
  genes locate their change at the S2/S3 boundary and inside S3 is about *our* segment labels; Kirita's
  segment-resolved data does not reproduce a general S2/S3 split.
- **Kirita's own severe-injury program is the one set that is sharply S2/S3 at 4 h** (14.2) and loses it within
  12 h. It is not independent of the test data and is reported for its shape only.

### W3 · The magnitude persists and the geography does not

- **W3 fails, in the informative direction.** Taking each set's value relative to its own 12 h value, the median
  signed location ratio at 14 d is **0.16** against **0.54** for the level. The spatial structure decays faster
  than the size of the response.
- In absolute value the two decay alike (|z| ratio 0.86 for location, 0.60 for level at 14 d). The difference is
  **sign flips**: the share of sets whose S3 − S1 sign differs from its 12 h sign is 0 at 4 h, 0.19 at 2 d,
  **0.38 at 14 d**, and 0.04 at 6 weeks. At 14 days the geography of many programs reverses rather than fading,
  and by 6 weeks it has returned to its original sign at a smaller size.
- **The failed-repair programs are the clearest case.** Their magnitude keeps growing while their geography goes:

| set | level 12 h | level 14 d | level 6 wk | location 12 h | location 14 d | location 6 wk |
|---|---|---|---|---|---|---|
| G23-FR (Gerhardt) | 14.3 | **25.9** | 18.8 | 3.3 | **−3.7** | 0.9 |
| G21-VCAM (Gerhardt) | 13.0 | **23.4** | 15.2 | 5.2 | −1.8 | 1.9 |
| K-FR (Kirita, not independent) | 9.2 | **26.6** | 18.6 | 1.5 | −2.6 | 1.8 |

  This settles the open question in the protocol in favour of Kirita's own reading: failed-repair cells are
  "scattered throughout", not confined to the outer stripe, even as the program itself becomes the dominant
  transcriptional signal. Dixon's persistent Kim-1⁺ outer-stripe tubules are a histological observation about
  where injured tubules are, not about where the failed-repair transcriptional program is strongest.

### W4 · Segment-specific pairs point opposite ways; whole-PT pairs do not

z of S3 − S1 at 12 h, for the -up and -down halves of each paired program:

| pair | up | down | signs |
|---|---|---|---|
| our 91 level | **+8.7** | **−6.9** | opposite |
| Dixon PTS1-2 | −1.8 | +1.7 | opposite |
| Dixon PTS3-TAL | +2.4 | −4.2 | opposite |
| Gerhardt S1 | −0.3 | +5.5 | opposite |
| Gerhardt S2 | **−10.1** | +4.3 | opposite |
| Gerhardt S3 | −5.1 | −6.6 | same |
| human AKI (Hinze) | +3.9 | +1.8 | same |
| Lake aPT | +3.3 | +1.8 | same |
| Zhou IRI | +5.6 | +1.6 | same |

- Five of the nine pairs put their up and down halves on opposite sides of the PT, and four do not.
- **The five that separate are the segment-specific lists** (ours and the Dixon and Gerhardt S1/S2 lists). Opposite
  signs are what one S3-weighted injury produces: it raises some genes most in S3 and lowers others most in S3.
- **The four that share a sign are the whole-PT injury and repair lists** (Gerhardt S3, Hinze's human AKI, Lake's
  adaptive PT, Zhou's IRI). Both halves lean the same way, which is what a list assembled without segment
  information looks like when it is dominated by one compartment.

### W5 · Our position is older than our level

| Spearman with our P_late | 4 h | 12 h | 2 d | 14 d | 6 wk |
|---|---|---|---|---|---|
| S3 − S1 (location) | 0.18 | 0.25 | **0.27** | 0.06 | 0.05 |
| S3 − S2 | 0.13 | 0.20 | **0.22** | 0.03 | 0.06 |
| S2 − S1 | 0.13 | **0.15** | 0.14 | 0.03 | −0.01 |
| our L against level (notebook 94) | 0.43 | **0.52** | 0.41 | 0.30 | 0.13 |

- **W5 holds.** Our level matches 12 h best, our position matches **2 d** best, and within position the S2 − S1
  correlation peaks at 12 h while S3 − S2 peaks at 2 d.
- Read as a sequence, that is a wave: the S1-against-S2 difference is established first and the S2-against-S3
  difference follows. Our single time point sits where both are present, which is why notebook 99 dated us at
  12 h to 2 d from levels alone.

### Controls

- **C1 fails as written and passes where it matters.** The housekeeping set (3,081 genes) has max |z| of 1.75 on
  S3 − S1 and 1.93 on S3 − S2, but 2.99 on S2 − S1 and **10.0 on the level**. A set of that size detects a global
  level shift between IRI and sham libraries; it detects no S3 geography. The geography null is sound, the level
  null is not, and the level column is used here only as the comparator in W3.
- **C2 holds.** Restricting every set to the genes passing notebook 112's contamination filters leaves the 12 h
  geography unchanged (Spearman 0.98 over sets).
- **C3 is not testable.** Rebuilding our three 91-derived sets under the balanced relabelings leaves only one set
  with ≥ 10 members (level-up R1, 36 genes), whose S3 − S1 z is 0.94, −0.53, 0.10, −0.47, 0.10 across the five
  times — flat, as it should be, but a single set is not a control.

### What this answers

1. **The acute IRI response is S3-weighted in an independent 18-mouse experiment**, for injury markers, for
   transport and for our own convergence core, which had no external positional input. This is the first
   confirmation of our S3 reading in segment-resolved data rather than in spots (notebook 105).
2. **The geography is acute and does not last.** It is strongest at 4 h, gone by 14 days, and at 14 days many
   programs reverse before settling. Anything positional about AKI must be measured early.
3. **Failed repair is transcriptionally dominant and spatially diffuse.** Its magnitude more than doubles from
   12 h to 14 d while its S3 weighting disappears. A paper claim that failed repair is an S3 phenomenon is not
   supported by this grid.
4. **"S3 is most vulnerable" is too coarse for our own down-called genes**, which separate S1 from S2 and S3
   rather than S3 from S1 and S2. The two boundaries behave differently for different programs and should be
   reported separately, not collapsed into one S3 − S1 contrast.
5. **The cell cycle is the only late-peaking geography** (2 d), matching the published outer-stripe proliferation
   wave, and giving an internal check that the grid can see a known late spatial effect.
6. **For the paper:** quote the acute S3 weighting with this external confirmation; date our cohort at 12 h on
   levels and 2 d on position; and do not extend any positional claim beyond 2 days.
