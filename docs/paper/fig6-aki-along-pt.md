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

## Results

### Placement of injured tubules (2026-11-02)

`analysis/paper/fig6_placement.py`, output in `results/paper_figures/fig6/placement/`. To add the arm, the rerun
script `fig3_robustness.py` gained an `exclude` condition, which withdraws a named gene list from the HVG pool.

**Injury-responsive genes.**
- 10,404 genes were assessed. 3,812 are injury-responsive: 2,896 up and 916 down.
- All 3,812 were withdrawn from integration.
- The injury-blind rerun kept 33,375 units (the reference keeps 32,939), and its segment order ρ is 1.0.

Spearman ρ of each section's PT positions:

| Section | all PT | S1 | S2 | S3 | gene half A vs B, all PT |
|---|---|---|---|---|---|
| Ctrl1A2 | 0.938 | 0.855 | 0.847 | 0.683 | 0.902 |
| Ctrl1A4 | 0.971 | 0.891 | 0.900 | 0.673 | 0.946 |
| IR2A2 | 0.898 | 0.682 | 0.821 | 0.657 | 0.822 |
| IR2A4 | 0.914 | 0.689 | 0.812 | 0.732 | 0.839 |

**Verdict: passes.**
- **All PT.** The floor is 0.888, and the IR sections reach 0.898 and 0.914.
- **S3.** The floor is 0.623, and the IR sections reach 0.657 and 0.732.
- Without 3,812 injury genes, the IR tubules keep their order nearly as well as the control tubules do.

**Not in the pass rule: within S1, IR agrees less.** ρ is 0.68–0.69 in IR against 0.86–0.89 in controls. Within-S1
positions in IR therefore lean partly on injury genes. The IR gene-half agreement within S1 (0.40–0.43) is also low,
but so is the controls' (0.47–0.55).

### Screen (2026-11-02)

`analysis/paper/fig6_screen.py`, output in `results/paper_figures/fig6/screen/`.

**Inputs.**
- PT units: Ctrl1A2 2,780, Ctrl1A4 2,745, IR2A2 3,937, IR2A4 5,551.
- 32 gene-specific axes; on each, b12 is 0.29–0.30 and b23 is 0.67–0.68.
- 5,379 genes were screened. The neighbour filter dropped 2,934; it is applied within each condition, and injured PT
  loses so much transport expression that many genes now exceed twice the PT rate in a neighbour group.

**Category A, the response changes inside a segment: not calibrated (null).**
- There are 34 calls, one per gene: S1 8, S2 12, S3 14; 13 go up.
- The decoy (positions shuffled within label) gives 4 calls. The false-call share is 0.118, above the frozen 0.10.
- So, by the protocol, the figure shows no Category A list.
- Not used for any decision: the symmetric relabeling calibration gives 1 and 3 calls, a share of 0.059.

**Display map.** 2,220 screened genes reach |Δ| ≥ 1 log2 somewhere in a core, with all four pairs agreeing in sign
there.

**Handoffs: 11 of 52 pairs shift; none is lost.**
- All 52 pairs (50 from notebook 101, 2 from the poster) cross in both controls on this fit.
- The 11 shifted pairs are:

| Pair | Control x\* | AKI x\* | Δx\* (95% CI) | What moves |
|---|---|---|---|---|
| Slc16a1 → Slc16a4 (MCT1 → MCT4) | 0.63 | 0.21 | −0.42 (−0.55, −0.26) | late partner loses its S3 rise; MCT1 stays on late |
| Slc16a1 → Slc16a9 | 0.62 | 0.31 | −0.30 (−0.45, −0.27) | late partner loses its S3 rise |
| Slc16a12 → Slc16a9 | 0.62 | 0.32 | −0.30 (−0.39, −0.24) | late partner loses its S3 rise |
| Slc6a19 → Slc6a13 (B0AT1 → GAT2) | 0.53 | 0.30 | −0.23 (−0.28, −0.20) | late partner loses its S3 rise |
| Slc46a3 → Slc46a1 | 0.66 | 0.44 | −0.22 (−0.44, −0.07) | late partner loses its S3 rise |
| Slc3a2 → Slc3a1 (4F2hc → rBAT) | 0.51 | 0.30 | −0.20 (−0.23, −0.16) | rBAT loses its S3 rise |
| Slc6a8 → Slc6a13 | 0.52 | 0.34 | −0.19 (−0.29, −0.14) | late partner loses its S3 rise; the early one falls later |
| Slc22a6 → Slc22a2 (OAT1 → OCT2) | 0.65 | 0.50 | −0.14 (−0.15, −0.13) | small; neither half-point moves much, the crossing does |
| Slc39a5 → Slc39a8 (ZIP5 → ZIP8) | 0.45 | 0.62 | +0.17 (+0.11, +0.22) | ZIP5 stays on further into S2 |
| Slc39a5 → Slc39a14 (ZIP5 → ZIP14) | 0.53 | 0.72 | +0.20 (+0.08, +0.24) | ZIP5 stays on further into S2 |
| Slc5a9 → Slc5a3 | 0.46 | 0.74 | +0.27 (+0.22, +0.35) | the two controls disagree on where Slc5a9 falls (0.36 and 0.95) |

- The "what moves" column is descriptive. It compares the early profile's last fall below 0.5 and the late profile's
  first rise above 0.5 (`pseudospace.handoffs.half_points`) between conditions.
- **Unshifted pairs, for comparison.**
  - Cystine (Slc7a9 → Slc7a13) keeps its position: Δx\* +0.07, not specific.
  - Its late partner Slc7a13 is 4.3 log2 lower around the handoff. The handoff stays where it was, at about a
    twentieth of the level.
  - The poster's glucose (Slc5a2 → Slc5a1) and lactate (Slc5a12 → Slc5a8) pairs do not move (−0.02 and +0.01).
  - So the poster's "the cystine handoff is lost" is a level statement, not a position one.

**Post-hoc check of the handoff calls (not in the protocol; added after the screen).**
- The protocol calibrated Category A and the pathways, but not the handoffs.
- **Method.** The shift rule was applied symmetrically, with each relabeling as the contrast and the other two splits
  as its specificity check. The bootstraps are those already drawn.
- **Result.** Relabeling 1 gives 0 shifted pairs and relabeling 2 gives 0, against 11 for the condition.

**Pathways within segments: 122 calls, none under the relabelings (with a post-hoc correction, below).**
- 1,070 sets were tested in each segment, 3,206 set × segment tests in all. 122 are called:

| Segment | Calls | Rising | Falling |
|---|---|---|---|
| S1 | 89 | 12 | 77 |
| S2 | 24 | 12 | 12 |
| S3 | 9 | 0 | 9 |

- **Calibration.** Relabeling 1 gives 0 calls, and so does relabeling 2. None is even at q ≤ 0.05, against 338 for the
  condition. The calibration share is 0.
- **Effect sizes are modest.** |W| has a median of 0.13 log2 and a maximum of 0.45 among the calls. The largest
  relabeling is at most 0.25 × the condition (median 0.17).
- **Themes.**
  - S1, falling: DNA replication and cell-cycle sets (Synthesis of DNA −0.16, DNA replication −0.16). The AKI rise of
    the replication program is largest at the start of S1 and eases along it.
  - S1, rising: integrin signalling (+0.25) and gluconeogenesis (+0.21).
  - S2, rising: DNA strand elongation (+0.42), mismatch repair (+0.40), G1/S-specific transcription (+0.32). The
    replication response climbs through S2.
  - S2, falling: SREBP-driven cholesterol synthesis (−0.14 to −0.17) and plasma lipoprotein remodelling (−0.31).
  - S3, falling: ether lipid metabolism (−0.45) and oxidative phosphorylation / respiratory electron transport
    (−0.20 to −0.22). These fall further through S3.
- **Caveat on S1.** The first third of the S1 core starts at x = 0, where splines are least constrained. Both
  conditions share the basis, and fewer than 0.1% of IR units are clipped at the start, so this is not an end pile-up.
  Still, the S1 calls are the least robust kind.

**Post-hoc correction to the pathway statistic (added after the first run; labelled).**
- **What went wrong.** In the first run (kept in `paper_figures/fig6/screen_run1/`), the screen called 0 sets under
  any split. The set statistics there reached 34–6,921 log2, which is impossible.
- **Cause.** About a dozen genes have no count in one third of a core in some mouse; for example, Akr1c18 in S1 and
  Cyp27b1 in S3. Their splines dive without bound there, giving within-core changes up to 6,921 log2. These few genes
  dominated both the observed set means and the random-set null.
- **Rule.** A gene's W in a segment is treated as unsupported when, in any mouse, the first or last third of that
  core holds no count of it. Such genes leave that segment's set means and its random-set pool. This removes 19, 17
  and 30 gene × segment entries in S1, S2 and S3.
- **Scope.** The rule touches only the pathway statistic and the program summaries. Category A, the display map and
  the handoffs reproduce the first run exactly.
- **Disclosure.** Because the rule was written after a null was seen, these calls are reported as post hoc.

**Programs on this axis (notebook 76's verdicts stand).**
- **Signs agree with notebook 76.**
  - Transport T_PT: S1 +0.08, S2 −0.03, S3 −0.25. The S3 relabelings are −0.12 and −0.12, so it is not specific.
  - Cell cycle C: S1 −0.10, S2 +0.25, S3 −0.10. The S2 relabelings are 0.07 and 0.04, 0.27×, so it is not specific.
- These are set means over a fixed list, without the random-set comparison of the screen. The screen can still call
  narrower replication sets in S2, where those sets change more than the whole C program does.
- The neighbour filter leaves only 3 of R's 18 markers (Havcr1 among them), so R is not shown.

### Figure (2026-11-02)

`analysis/paper/fig6_panels.py`.
- **a.** Ctrl1A2 and IR2A2, PT tubules coloured by position.
- **b.** Density along the control-quantile axis, per section.
- **c.** Placement check.
- **d.** Havcr1, Krt20, Mki67 and Slc34a1. A mouse's curve is drawn only where 16 of its 20 bins have counts, and it
  is clipped to its bin range ± 0.25.
- **e.** C and T_PT profiles with pairs, relabelings and segment means.
- **f.** The 2,220-gene display map beside the segment view.
- **g.** Scaled profiles of three handoffs: one moving upstream, one moving downstream, and cystine, which does not
  move.
- **h.** All 52 handoff positions, and the 11 that move, per mouse.
- **i.** All 122 pathway calls, each centred on its segment mean.

**Presentation choices.**
- **Examples in g.** They were chosen after the screen as one of each kind.
- **Panel d.** Slc7a13 was swapped for Slc34a1 so no gene appears twice.

### Figure, second version (2026-11-02, after the user's review)

The user found the handoff and pathway panels hard to read. Both were redrawn in the layout of the control vs AKI
poster.
- **g, handoffs that move.** The 11 shifted pairs, each from its control switch to its AKI switch, with each mouse
  shown. Protein names are added. The control-vs-AKI scatter of all 52 pairs is dropped.
- **h, five handoffs in poster layout.** Each row gives the substrate, the genes and the proteins, with the control
  and AKI curves side by side and the segment bands behind them. The control switch is dashed and the AKI switch
  dotted. Bars show each partner's level in AKI.
  - **The y axis is scaled expression, not the poster's log counts.** The switch is defined on scaled profiles (notebook
    101), so the marked switch is exactly where the drawn curves cross.
  - On log counts, a level drop moves the crossing too. That is how the poster showed cystine as "lost" and moved
    lactate and glucose. On scaled profiles, none of those three pairs moves (Δx\* +0.07, +0.01 and −0.02).
  - The level bars carry the cystine result: Slc7a13 is 4.3 log2 lower.
  - **Rows (chosen after the screen):** MCT1 → MCT4 (largest shift), B0AT1 → GAT2 and 4F2hc → rBAT (upstream shifts,
    amino acids), ZIP5 → ZIP14 (downstream shift) and cystine (no shift, level loss).
- **i, four pathways in poster layout.** Each shows the S1/S2/S3 means (called segments outlined) beside the PT drawn
  as a tube. The tube is coloured by the pathway's AKI response along it, from weakest to strongest within the
  pathway, and a circle marks the strongest change.
  - The four were chosen after the screen as one readable representative per theme of the 122 calls:
    - S phase: proliferation, falling within S1;
    - integrin cell-surface interactions: rising within S1 and falling within S2, a peak at the S1/S2 boundary;
    - cholesterol biosynthesis: falling within S2;
    - oxidative phosphorylation: falling within S3.
  - Bile acid metabolism was tried first and dropped, because its within-S2 change (−0.08) is not visible.
  - The example genes are the three members whose change in the called segment is largest in the call's direction.

### Main figure vs Extended Data (2026-11-02, user: the main figure keeps the biggest findings)

- **Main Fig. 6** (`make_fig6.py`):
  - a. Programs along the PT.
  - b. The 11 handoffs that move.
  - c. Five handoffs in the poster layout.
  - d. Four pathways on the nephron schematic.
- **Extended Data** (`make_ed_fig_aki.py`):
  - a. Control and AKI sections coloured by position.
  - b. Coverage along the axis.
  - c. Placement check.
  - d. Injury and transport genes (Havcr1, Krt20, Mki67, Slc34a1).
  - e. The 2,220-gene response map beside the segment view.
- **Panel outputs** of `fig6_panels.py` are now named by content (`f6_<name>.png`), not by letter.

### Third version (2026-11-02, user: b and c were redundant; the nephron schematic was for the poster only)

- **Main Fig. 6.**
  - a. Programs along the PT, stacked.
  - b. Pathway heatmap.
  - c. Five handoffs in the poster layout.
- **Extended Data f.** The 11 moving handoffs, control to AKI, moved here from the main figure.
- **Pathway heatmap (main b).**
  - **Rows.** One row per distinct called pathway, 33 in all.
    - Calls in a segment that share 0.3 or more of their screened genes (Jaccard) are drawn once.
    - The representative is chosen in library order (Hallmark, Reactome, then KEGG), then by p.
    - A pathway called in two segments (integrin cell-surface interactions, in S1 and S2) is one row.
    - KEGG "Human Diseases" sets (for example Alzheimer disease and small cell lung cancer) are not drawn. In the PT
      they stand in for respiratory-chain or cell-cycle members. They remain in `pathways.csv`.
  - **Blocks.**
    - The AKI response along the PT.
    - The same pathways' S1/S2/S3 means.
    - The response centred on its mean within each segment core, with the transition windows left blank. This block
      shows the change the segment means cannot express.
  - **Layout.** Called segments are outlined, and rows are grouped by segment. The row list is written to
    `panels/f6_pathways_rows.csv`.
  - All of these are presentation choices made after the screen.
