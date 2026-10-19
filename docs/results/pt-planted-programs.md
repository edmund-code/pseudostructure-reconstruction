# From called sets to programs, with a planted truth (notebook 102)

**The question.** Notebook 100 found that every calling method spills a planted pathway difference into the sets that
share its genes. With 3 planted sets at 1 log2, the positional F-test calls about 200 other sets, almost all sharing
planted genes. So the number of called sets is not the number of findings, and the calls must be grouped into
programs before they are read.

Three groupings are in use:
- notebook 71's redundancy clusters (RED);
- notebook 75's overlap programs (OVL);
- notebook 93's driver-gene networks (NET).

Notebook 93 compared them without a truth. Here the truth is planted. Which grouping turns the positional F-test's
calls back into the planted programs, one program per planted set, with the fewest programs that contain no planted
set?

Protocol frozen on 2026-10-18, before any planted call set was grouped. Notebook 100's sensitivity and false-call
counts had been seen; no grouping of its calls had been computed.

## 1 · Labels and inputs

- Notebook `analysis/notebooks/102_pt_planted_programs.ipynb`; logic `102.planted_programs.1`; seed 100. Notebook
  100's seed is kept, so the planted data are identical.
- Results `results/pt_planted_programs/`.
- No new module. It uses `planting`, `pathway_pipelines`, `pathway_decoys.overlap_programs`, and
  `pathway_networks.driver_jaccard` and `leiden_programs`.
- **Replicates:** notebook 100's two arms with identical seeds and planting:
  - the frozen arm: 15 sets × 4 scenarios × 2 amplitudes × 8 replicates;
  - the low-dose arm: 3 sets, the same grid.
- **Method:** the positional F-test only, notebook 100's code. Its planted genes are refitted and spliced, and this
  time the per-gene z and q are kept.
- **Redundancy clusters:** notebook 71's `redundancy_clusters_aki_vs_control.csv` (membership Jaccard ≥ 0.5,
  connected components; its `representative` column).

## 2 · Groupings (each as first defined)

| Grouping | Definition |
|---|---|
| RED | the calls' redundancy-cluster labels (notebook 71) |
| OVL | `overlap_programs` on the calls. Contributing genes are the top 10% of z (notebook 75). Average linkage on 1 − overlap coefficient, cut 0.5. |
| NET | `driver_jaccard` on the calls' driver genes, then `leiden_programs` (threshold 0.1, resolution 1, seed 1234). Isolates are kept as one-call programs. Driver genes follow notebook 93's F rule: members with gene q ≤ 0.05; if fewer than 3, members in the top 10% of z; if still fewer, all members. |
| NONE | each call is its own program (the reference) |

## 3 · Metrics, per replicate and grouping

- **Programs:** the number of programs among the calls.
- **True programs:** programs that contain ≥ 1 planted set.
- **False programs:** programs with no planted set.
- **Separation:** the number of distinct programs holding the called planted sets, divided by the number of called
  planted sets. It is 1 when each planted set sits in its own program. The planted sets were drawn with Jaccard < 0.1
  and no shared planted gene, so they are distinct programs by construction.
- **Containment:** the share of all calls that sit in a true program. A spilled set grouped with its planted source is
  absorbed into a true program.
- **Baseline false programs:** the programs the grouping makes from the 29 baseline calls of this partition (no
  planting), reported once as the floor.

## 4 · Readings (fixed now)

- **Best grouping per arm and scenario (A = 1):** the grouping with the fewest false programs whose mean separation is
  ≥ 0.8.
  - If none reaches 0.8, the grouping with the highest separation is reported, and the arm is marked "no grouping
    separates".
  - Groupings within 1 false program of each other tie.
- The same is reported at A = 0.5, and descriptively for the low-dose arm.
- **The readout for the paper:** after the best grouping, how many false programs remain per planted program. This
  converts notebook 100's set-level spillover into program-level terms.

## 5 · Outputs and figures

- **Tables:** `replicate_programs.csv` (per replicate × grouping: programs, true, false, separation, containment),
  `summary.csv`, `readings.csv`, `baseline_programs.csv`, `how_the_protocol_was_read.csv`, `run_record.json`.
- **Figures (plots only):**
  - 102a: false programs against separation per grouping, both arms, A = 1 and 0.5;
  - 102b: programs per replicate by grouping, against the planted 15 or 3.

## 6 · Caveats

- Only the positional F-test's calls are grouped. Notebook 93's driver rule exists for clustering too, but would need
  clustering's per-list overlap genes for every replicate.
- The planted programs are disjoint by construction. Real programs overlap, so separation here is an upper bound.

## Notebook 102 results

`analysis/notebooks/102_pt_planted_programs.ipynb` (logic `102.planted_programs.1`) writes to
`results/pt_planted_programs/`. It followed the protocol frozen in 81b1dd0.
- G1–G2 passed.
- G4 passed: the rebuilt F-test calls equal notebook 100's on all 8 checked replicates, with the same planted sets.
- 128 planted replicates (2 arms × 4 scenarios × 2 amplitudes × 8).
- The first run's figure call used a renamed matplotlib argument. It was fixed before execution, and the results are
  unchanged.

**Baseline (29 calls, no planting):** NONE 29 programs, RED 27, OVL 20, NET 7.

**At A = 1** (means over 8 replicates; separation 1 = each called planted set in its own program):

| Arm · scenario | Calls | Programs NONE / RED / OVL / NET | False programs RED / OVL / NET | Separation RED / OVL / NET | Calls inside true programs, NET |
|---|---|---|---|---|---|
| 15 sets · shape | 690 | 690 / 437 / 112 / 9.5 | 422 / 97 / **2.6** | 1.00 / 1.00 / 0.46 | 0.96 |
| 15 sets · located | 660 | 660 / 406 / 108 / 9.6 | 391 / 93 / **1.8** | 1.00 / 1.00 / 0.52 | 0.98 |
| 15 sets · bump | 689 | 689 / 438 / 116 / 9.0 | 423 / 101 / **1.8** | 1.00 / 1.00 / 0.48 | 0.95 |
| 3 sets · shape | 244 | 244 / 155 / 53 / 7.9 | 152 / 50 / **5.3** | 1.00 / 1.00 / **0.88** | 0.66 |
| 3 sets · located | 196 | 196 / 140 / 46 / 9.1 | 137 / 43 / **6.3** | 1.00 / 1.00 / **0.96** | 0.66 |
| 3 sets · bump | 234 | 234 / 163 / 52 / 8.6 | 160 / 49 / **5.9** | 1.00 / 1.00 / **0.92** | 0.63 |

The F-test calls no planted set in the level scenarios, so those rows hold the baseline calls only.

**Readings (frozen rule: the fewest false programs with separation ≥ 0.8).**
- **Low dose (3 sets, the dose of one or a few programs):** NET wins every scenario at both amplitudes. Hundreds of
  called sets become 8–9 programs.
  - The planted programs stay separate (0.83–1.0).
  - About 5–7 false programs remain. That is the baseline's own 7 NET programs, so NET adds almost nothing beyond the
    baseline.
- **Heavy dose (15 sets):** OVL wins by the rule. NET collapses about 690 calls into 9–10 programs with only 2–3 false
  ones, but it merges the 15 planted programs (separation 0.46–0.52). RED and OVL keep planted sets apart but leave 93–420
  false programs.

**What this answers.**
1. **Counting called sets overstates findings by one to two orders of magnitude.** Every grouping helps, and RED, the
   grouping used so far, helps least: it removes about 40% of the spilled calls.
2. **Driver-gene networks (NET) are the best way to read a call list at realistic doses.** They recover one program
   per planted program and add almost no false programs beyond the baseline. This upgrades notebook 93's verdict on
   NET, which was "mixed, settings-sensitive" without a truth.
3. **NET's limit:** when many strong programs change at once, it merges distinct programs. OVL keeps them apart at the
   cost of about 6 false programs per planted program. Read NET's programs beside OVL's when a contrast is very strong,
   as the AKI condition is.
4. **For the paper:** report positional F-test calls as NET programs, with OVL as the finer view. Never quote raw set
   counts as numbers of findings.
