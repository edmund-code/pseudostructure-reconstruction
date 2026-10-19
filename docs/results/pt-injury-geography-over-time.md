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
