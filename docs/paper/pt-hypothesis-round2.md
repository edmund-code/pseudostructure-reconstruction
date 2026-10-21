# Hypothesis round 2: what the cohort supports, on statistics that survived the audits

2026-10-21. Six data-driven screens plus four physiology-led tracks, each attacked twice, then a
prioritiser that was told to rank by what a lab can test rather than by significance. The round was
designed around the one statistic that survived the species audit — the **pure-compartment
contrast**, which uses no curve, no registration and no centre of mass — so the amplitude-as-position
failure that destroyed the previous round cannot act on the primary result.

The round's best feature is that it removed three of the project's own five surviving reversals,
including the one ranked first.

## 0 · The screen, re-derived from scratch by me

Statistic, per specimen: `E = mean over pure-S3 bins of log1p(rate per 10k) − mean over pure-S1 bins`,
with pure bins taken at ≥ 0.90 label purity from `composition.csv` and a ≥ 20-unit guard. A gene is
a **reversal** when both specimens of each arm agree in sign, the arms disagree, and the weakest
specimen's |E| clears a threshold τ. The gene axis was asserted identical between the statistics
table and the binned-rate cache (11,055 rows, `True`).

I recomputed the whole screen independently, including the nulls, and **chose a different operating
point than the round did**:

| τ | cross-species reversals | mouse-vs-mouse null | human-vs-human null | mixed-arm nulls |
|---|---|---|---|---|
| 0.10 | 190 | 2 | 37 | 0, 1 |
| 0.15 | 105 | 0 | **2** | 0, 0 |
| **0.20** | **56** | **0** | **0** | **0, 0** |
| 0.25 | 32 | 0 | 0 | 0, 0 |
| 0.30 | 19 | 0 | 0 | 0, 0 |

The round reported τ = 0.15 as "the smallest threshold at which both within-arm nulls are empty".
That is true only *after* a six-filter stack, one filter of which is derived from an external
single-cell atlas. On the raw statistic the human-section null still fires twice at τ = 0.15, so
**τ = 0.20 is the honest operating point**: both within-arm nulls are empty with no filtering at all.

**The result at τ = 0.20 with internal-only filters** (detection ≥ 0.25 both species, not an
axis-construction gene, edge mass < 0.25 both arms): **47 genes**, and **47 of 47** pass an
orientation-free re-check — Spearman of the gene's profile against bin index across all 20 bins, no
purity selection, opposite-signed in all four specimens.

```
Abhd17b Abr Acadsb Acbd4 Acmsd Aco1 Acot11 Agpat3 Atp1a1 Atp6v1b2 Cat Cox7a1 Ctsh
Cxcl14 Cyfip2 Dbp Dcxr Ddit4 Ddr1 Defb1 Enpp6 Hipk3 Hpd Kif12 Lap3 Mfap3l Miox Mme
Msrb1 Odc1 Pck1 Pitpnc1 Rbpms Reep6 Rnf128 Sgk1 Sh3gl2 Slc16a12 Slc16a9 Slc17a1
Slc17a3 Slc38a10 Slc5a11 Tmem252 Tmem59 Ugt3a1 Vnn1
```

Two caveats on this list that are not optional. **Miox and Cox7a1 are axis-construction genes** under
the other table's polarity and are carried as context only (open item 34). And the list has **not**
had a PT-identity restriction applied, which is the filter that mattered most below.

## 1 · What the round removed, including its own front-runner

| Gene | Round-1 status | Verdict now | What killed it |
|---|---|---|---|
| **Ugt3a1** | reversal #1 | **removed** | The only many:many ortholog in the set, and its mouse paralog Ugt3a2 is absent from the shared panel, so paralog substitution cannot be excluded internally. 94% of the human effect is attributed to thick-limb admixture, and the sign flips to +0.392 extrapolated to zero differential. Human distal bins hold 5–24 counts; detection_human 0.488. |
| **Hao2** | "the one gene genuinely invisible at the ends" | **removed** | On PT-identity units the mouse sign **flips** (−0.411 → +0.561/+0.626), retention 0.047 — its "mouse proximal weighting" was thin-limb units sitting inside the PT-S3 label. The committed fitted statistic had called it conserved all along, and that disagreement is now explained. |
| **Abhd14a** | reversal | **removed** | detection_human 0.256–0.261, at the gate, with a human atlas neighbour ratio of 2.80. Independently flagged by the admixture screen (F2c). |
| **Abhd17b** | reversal | **removed** | detection 0.351/0.351, both near the floor; minimum \|E\| 0.208. |
| **Dcxr** | reversal | **kept, promoted** | Survives PT-identity restriction (retention 0.861), is injury-independent (frac_inj −0.055), compartment-local detection ≥ 0.73 in all eight compartments. |
| **Slc17a3** | reversal | **kept, demoted and reframed** | frac_inj +0.737 — mouse injury alone reproduces 74% of the species shift. And the mouse shape is a **mid-axis peak**, so the mouse fact is exclusion-from-S1, not concentration-in-S3. |

Also killed, each with the measurement: the oxalate/glyoxylate trio with Hao2 (Agxt detection_mouse
0.001, Slc26a6 0.089, Slc5a8 detection_human 0.072 — untestable regardless); the Slc26a1/Sat1
inversion (control distal detection 0.116 on 54/60 counts, the same floor that banned Slc5a1); the
Acy1⇄Acy3 "baton pass" (no positional content — member shares barely move, and both are same-sign in
both arms); and Slc6a18→Slc6a19 as a species handoff (same sign in both species).

**Ugt3a1's removal is the round's most valuable output**, because the project had it ranked first and
the literature pass had independently called its human limb unmeasured. Both were right: nobody has
measured it, and we cannot either.

## 2 · What the cohort cannot answer, with the number that proves it

**The collaborators' SGLT2-inhibitor question is unanswerable on this panel, and now it has a
powered null rather than an absence.** Every gene in the SGLT2→SGLT1 handoff and the distal escape
route fails: Slc5a1 detection_human 0.061 (mouse neighbour ratio 7.12), Slc2a1 0.110/0.098,
Slc2a9 detection_mouse 0.083, Abcg2 detection_human 0.060, Slc16a1 detection_mouse 0.140,
Slc5a8 detection_human 0.072; Slc22a11/OAT4 is absent from the shared space entirely; and Slc5a2
itself is an axis-construction gene, so SGLT2's own position on this axis is partly circular.

But the testable half gives a real answer: **SGLT2 → GLUT2 is identically positioned in the two
species** — entry/exit matching index −0.393 mouse against −0.396 human, ΔD = −0.003, against an
amplitude surrogate that by itself generates a median |ΔD| of 0.383 with 48.6% of 4,000 draws
exceeding the observation. The two genes' human/mouse amplitude ratios are 0.533 and 0.531,
identical. So the proximal glucose reabsorptive unit is complete and proximally weighted in **both**
species with no species difference in entry/exit matching; the distal escape route is unmeasured; and
**a positional explanation for species differences in SGLT2i response is not supported by this
cohort.** That is a more useful thing to tell the collaborators than silence, and it names the assay
that would answer it: one that actually detects SGLT1 in human PT, scored on straight versus
convoluted profiles.

**Whether any human effect is donor-level or species-level.** The two human sections are two sections
of one block, and they agree **worse than two separate animals**: 44 gene pairs differ between the
human sections at |Δz| ≥ 0.3 against 4 between the two mice. Every "4-of-4 specimen agreement" in
this project is 2 biological mouse replicates plus 2 technical human sections. This is the single
most important sentence in the round and it belongs in the paper's limitations verbatim.

## 3 · The ranked hypotheses

Eleven, each with a named kill condition. The three strongest, with their status:

1. **PCK1 — human PT makes PEP at the straight segment.** Mouse E −0.354/−0.255, human +0.743/+0.864;
   detection 0.98/0.99, the best-detected gene in the list; PT-identity retention 0.808. The
   mechanism is the dissociation: the gluconeogenic trunk (Fbp1, Aldob) stays proximal in both
   species and does not move, so the human straight segment would be doing cataplerotic carbon
   disposal rather than completing glucose synthesis. **Falsified if** human PCK1 is flat or
   convoluted-weighted in Umod-negative, Lrp2-positive profiles, or if Fbp1 moves distally with it.
   **A hazard the round did not weigh:** the published rodent record (abstract-only) reports PEPCK
   *highest* in S3, which is the opposite of our mouse limb. That disagreement must be resolved
   before PCK1 is published as a reversal — and the stain already proposed is what resolves it.
2. **SELENOP — injury partitions tubules rather than dimming a field.** Dispersion rises 12.5× while
   the mean moves 0.84×, at 42 counts per unit with detection 0.998, and the depth profile is flat in
   both arms. Abundance-matched comparators at similar dispersion ratios all show a 1.6–3.4× mean
   fall. **Falsified if** per-profile CV in IRI is ≤ 15% or unimodal, or if per-profile CV falls as
   1/√cells, which would make it cell-level averaging rather than a tubule-level state.
3. **CROT — ischaemia deletes the positional induction, not the setpoint.** Mouse E +1.839/+1.989,
   human −0.024/−0.068, with a co-deleted S3 identity module (Eci3, Adhfe1, Acox2, Haao, Hgd, Pipox,
   Dao, Amacr). **Falsified if** more than ~10% of IR outer-stripe PT profiles retain control-intensity
   CROT, or if cortical PT loses it as much as the outer stripe.

Then: PDZK1 and the apical urate transportasome (moderate); DCXR (moderate, and now a *validation* —
see `literature-novelty-ledger.md`); LGMN front-loading lysosomal proteolysis (moderate); **KIF12**
(moderate — a new reversal, mouse −0.471/−0.503 against human +0.384/+0.435, detection 0.87/0.86,
PT-identity retention 0.943, and a direct HNF1B target, which gives it a named upstream cause);
SLC17A3 reframed as tissue-state-sensitive (speculative); ACMSD and the NAD branch point
(speculative); ZFP36 and the loss of post-transcriptional brakes (moderate); CSF1R (contradicted by
literature — its signal may be interstitial macrophages, and it carries a mandatory
macrophage-abundance channel).

**A surprise neither track promoted: QPRT also reverses**, found while checking ACMSD — mouse
+0.127/+0.185 against human −0.296/−0.434 on PT-identity units, and independently on the committed
fitted statistic, with detection 0.656/0.567 and balanced probes. With HAAO conserved distal in both
arms and ACMSD reversed, the NAD branch point now has **three** genes pointing the same way instead
of one. That is the most interesting unexamined thing in the round.

## 4 · The round's own verdict on its design, which should be quoted rather than paraphrased

**The tracks were not ten independent instruments.** Seven of the ten used the same two-end contrast
on the same cache through different hand-curated gene lists, so the six-track "convergence" on the
urate module is close to arithmetic — one measurement reported five times. That is precisely the
error the convergence ladder was just caught making (`aki-axis-thesis-audit.md` §4), arrived at
independently, and it is the reason the round's own prioritiser discounted it.

**The labels are not anatomy.** Every segment call in this dataset is a pass-1 Leiden partition on
expression named from a nine-gene marker panel, so those nine genes are circular for any label-based
contrast. Hence the panel's scoring rule, which is non-negotiable and correct: **score per tubule
cross-section, assigned by anatomy alone** — glomerulus-attached convoluted profiles versus
medullary-ray and outer-stripe straight profiles — with no marker-derived segment call entering the
scoring. Report each species' own straight/convoluted ratio, never a cross-species magnitude.

**The compartment is not the cell type.** PT-S3-labelled units are 26.6% / 26.2% proximal tubule in
mouse and 48.2% / 31.8% in human; purity in `composition.csv` is purity in the *label* (its columns
sum to exactly 1.000) and constrains nothing about cell type. This is why PT-identity retention is a
required filter rather than a bonus: it killed Hao2 at 0.047 and is the reason to trust Pck1 (0.808),
Dcxr (0.861) and Kif12 (0.943).

**The anchor is Slc4a4 and the banned anchor is Slc34a1**, both as already recorded. Slc4a4's
human/mouse amplitude ratio (0.428) sits squarely inside the cohort-wide attenuation band, which is
exactly what an anchor must do — it marks position and not a species effect.
