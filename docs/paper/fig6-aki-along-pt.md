# Fig. 6: control against AKI along the proximal tubule — protocol and results

Fig. 6 is the paper's summary of the control vs ischaemia–reperfusion (IR) comparison. It shows what the shared
pseudostructure axis lets us say about injury along the PT that a comparison of S1/S2/S3 means cannot: where along the
PT each response sits, how it changes inside a segment, and where paired transporters hand off. It reuses notebook 76's
frozen programs and verdicts (`docs/results/pt-package-aki.md`) and adds the analyses below.

## Frozen protocol (2026-11-02, before any fit)

**Why these analyses.**
- The abstract claims "injury-associated shifts in transporter switch points and pathway responses within segments".
  Neither has been measured on the package fit: the poster's handoffs were hand-picked, and notebook 76 tested
  programs only.
- A comparison along a shared axis is only meaningful if injured tubules are placed by their position and not by
  their injury. The handoff (section 6, item 11) asks for this check.
- Notebook 76 and the audit (`aki-axis-thesis-audit.md`) set the limits: the S3 weighting of injury is visible in
  segment means alone, and notebook 76's transport-loss and within-segment program claims failed specificity under the
  two relabelings. Nothing here reopens those verdicts; this protocol asks different questions.

### 1 · Data

- **Fit.** The control/AKI package fit (`runs/aki`), all four sections, and only units with a reference label PT-S1,
  PT-S2 or PT-S3.
- **Counts.** Raw counts from `02_qc/units_qc.h5ad`; each unit's library is its total over non-mitochondrial genes.
- **Contrasts.** For any per-mouse quantity Q:
  - condition: mean(Q of IR2A2, IR2A4) − mean(Q of Ctrl1A2, Ctrl1A4);
  - the four pairs: Q(IR) − Q(control);
  - relabeling 1: mean(IR2A2, Ctrl1A2) − mean(IR2A4, Ctrl1A4) (the A2/A4 factor);
  - relabeling 2: mean(IR2A2, Ctrl1A4) − mean(IR2A4, Ctrl1A2).
- **Specific** means both relabelings' |value| ≤ 0.25 × the condition's |value| (notebook 76's rule).

### 2 · Placement of injured tubules (the confound check)

**Question.** Are IR tubules placed on the axis by genes that do not respond to injury?

- **Injury-responsive genes.** From whole-PT pseudobulk (log2 counts per 10k, summed counts), a gene is responsive when
  all four pair differences have the same sign and each is at least 0.5 log2 in size. Genes detected in fewer than 5%
  of PT units in every section are not assessed (they cannot drive the integration).
- **Injury-blind axis.** One rerun of stages 03–04 with the Fig. 3d machinery (`fig3_robustness.py`: R harmony,
  4 threads, the run's own configuration). The responsive genes are withdrawn from the integration HVG pool through
  the `exclude_from_harmony_hvg` hook, exactly as a gene-half rerun withdraws half the genes.
- **Readouts.**
  - Per section, Spearman ρ between the reference position and the injury-blind position, over PT units.
  - The same within each PT label.
  - Reference for the noise of a rerun: the mean ρ between the gene-half axes A and B of the five Fig. 3d splits, per
    section.
- **Pass.** In each IR section, ρ(reference, injury-blind) is at least the lower of the two control sections' ρ minus
  0.05, both over all PT units and within S3.
- **If it fails,** IR placement depends on injury genes. Then the figure says so, and sections 3–5 are repeated on the
  injury-blind axis, which becomes the axis of the figure.

### 3 · Gene responses along the PT

**Axis.** As Fig. 4, extended to the IR sections:
- Each gene's position for a unit is the mean of the five Fig. 3d gene-half axes built without that gene; a unit needs
  at least 3 of 5.
- **Scale x: control quantile.** A control unit's x is its rank among the control PT units of both mice (Fig. 4's x
  exactly). An IR unit's x is its mean position mapped through the control units' empirical distribution (linear
  interpolation), clipped to [0, 1]. The share of IR units clipped at each end is reported.
- **Boundaries and cores:** Fig. 4's, from the control units (cores more than 0.06 from b12 and b23).
- **Spline knots:** at control x quantiles (Fig. 4's knots); the same basis in all four mice.

**Genes.**
- Mitochondrial genes are left out.
- **Screened:** detected in at least 10% of PT units in each mouse of at least one condition. This admits genes that
  injury switches on or off.
- **Neighbour filter:** Fig. 4's rule (dropped if AL, DCT or CD is more than twice the PT rate), applied within each
  condition's sections; a gene failing it in either condition is dropped.

**Models.** Per gene and mouse, NB GLMs with a log(library) offset: curve (the spline), clusters (the S1/S2/S3 factor)
and bulk, with the dispersion estimated under clusters + curve (Fig. 4's code). Curves are evaluated on Fig. 4's grid
of 201 points, in log2 counts per 10k.

**Response.** A gene's AKI response along the PT, Δ(x), is the condition contrast of its curves. The segment view of
the same gene is the condition contrast of its cluster means, one value per segment.

**Category A: the response changes inside a segment.** For a segment core, W is the mean of Δ over the core's last
third minus its mean over the first third. A gene is called in a segment when all hold:
1. |W(condition)| ≥ 1 log2.
2. All four pair W have the condition's sign and |W| ≥ 0.5.
3. W is specific (both relabelings ≤ 0.25 × the condition).
4. **Raw counts within the label.** In each mouse, the change in log2 rate (summed counts, pseudo-count 0.5) between
   the core's first and last thirds, over units carrying that segment's label. Its condition contrast is ≥ 0.5 log2
   with the sign of rule 1, and all four pair contrasts share that sign. (This is Fig. 4's post-hoc check, adopted here
   before any fit.)
5. **Support.** Each third of the core holds at least 25 units of that label in every mouse.

- **Calibration.** A decoy shuffles x within each label within each mouse (seed 0), refits the curves and applies
  rules 1–5. The estimated false-call share is (decoy calls) / (real calls).
- **If that share exceeds 0.10,** Category A is reported as not calibrated and the figure shows no Category A list.
- **Ranking for examples:** the smallest |pair W|.

**Display (descriptive, no claim):** a gene × x map of Δ(x) for genes whose |Δ| reaches 1 log2 somewhere on the core
grid with all four pairs agreeing in sign there, ordered by hierarchical clustering of Δ(x), beside their segment-view
values.

### 4 · Transporter handoffs

- **Pairs (frozen elsewhere, before any AKI handoff on this fit):** notebook 101's mouse handoff pairs confirmed by
  Chen 2021 microdissection, 50 of 52 (`pt_transporter_handoffs/tables/pairs.csv`). They were defined on the human
  vs mouse fit, from control mice only.
- **Curves:** each gene's per-mouse curve from section 3; pair members outside the screened universe are fitted the
  same way for this purpose only.
- **Profiles and crossing:** `pseudospace.handoffs`, as notebook 101. Each curve is min–max scaled over the grid points
  inside that mouse's 2.5–97.5% range of x; the crossing x\* is the first point where early minus late turns from
  positive to non-positive.
- **Eligible:** a crossing in both control mice on this fit.
- **Shift.** Δx\* is the condition contrast of x\*. A pair is **shifted** when:
  - all four pair Δx\* share a sign;
  - |Δx\*| ≥ 0.05;
  - Δx\* is specific;
  - the 95% interval of the condition Δx\* from 200 bootstraps excludes 0. Each bootstrap resamples units within each
    mouse and refits the two curves.
- **Lost:** no crossing in either IR mouse, against a crossing in both controls.
- **Level, reported beside each pair:** the condition contrast of each member's curve averaged over x\* ± 0.1 of the
  controls.
- **Poster pairs.** The poster's lactate (Slc5a12 → Slc5a8) and glucose (Slc5a2 → Slc5a1) pairs are not in notebook
  101's list. They are reported the same way, marked as previously shown, and carry no claim.

### 5 · Programs and pathways along the PT

**Programs (no new verdicts).**
- These are notebook 76's frozen programs: transport T_PT, cell cycle C, and the injury markers R.
- Their profiles are recomputed on this figure's axis, as the mean over members of the condition Δ(x), with pairs and
  relabelings shown, beside the per-segment values.
- **Claims are notebook 76's:** only the cell-cycle level (A4) and the positive control (A2) are claimable. A
  recomputed statistic whose sign differs from notebook 76's is reported.

**Pathways within segments (a screen).**
- **Sets:** notebook 76's 1,495 (`gene_sets.json`), restricted to screened genes, keeping sets with at least 10
  members.
- **Statistic:** a set's W in a segment is the mean of its members' gene W (condition), with the per-segment value
  reported beside it.
- **Null:** 1,000 random gene sets of the same size from the screened universe. Two-sided empirical p, BH across sets
  × 3 segments.
- **A call needs all of:**
  - q ≤ 0.05;
  - all four pair set-W share the sign;
  - set-W is specific.
- **Calibration:** the same screen run on each relabeling contrast in place of the condition. The calibration share
  is (relabeling calls, mean of the two) / (condition calls).
- **If that share exceeds 0.10,** pathway calls are reported as not calibrated and carry no claim.
- **Known limit:** random sets ignore gene–gene correlation. That makes the null anti-conservative, so the relabeling
  calibration is the binding check.

### 6 · What would be a null

| Analysis | Null outcome | What happens then |
|---|---|---|
| Placement | ρ check fails | figure moves to the injury-blind axis (section 2) |
| Category A | no call, or calibration > 0.10 | figure shows the Δ(x) map only, with no within-segment claim |
| Handoffs | no shifted or lost pair | "handoff positions are preserved after IR" is the result |
| Pathways | no call, or calibration > 0.10 | the abstract's "pathway responses within segments" is withdrawn |

No threshold changes after seeing results. Any rule added afterwards is labelled post hoc here.

### 7 · Known limits

- Two mice per condition; the relabelings are the only null for between-mouse differences, and they carry the A2/A4
  factor.
- **Time point.** The IR mice match about 12 h–2 d after IRI (notebook 99); one time point.
- **Deep S3 is thinned in IR.** Notebook 76 found that the IR mice nearly empty the ninth control decile. The support
  rule (rule 5) and the clipping count report it.
- **Positions within segments are noisy** (notebook 72), which flattens real within-segment changes; Category A is
  conservative.
