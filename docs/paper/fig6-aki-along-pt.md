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
| Slc16a1 → Slc16a4 (MCT1 → MCT5) | 0.63 | 0.21 | −0.42 (−0.55, −0.26) | late partner loses its S3 rise; MCT1 stays on late |
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
  - **Rows (chosen after the screen):** MCT1 → MCT5 (largest shift), B0AT1 → GAT2 and 4F2hc → rBAT (upstream shifts,
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

## Follow-up: what each result is, and which are notable (2026-11-02, post hoc, descriptive)

The user asked to investigate each pathway and gene and to pick the most notable results.
`analysis/paper/fig6_interpret.py` writes three tables to `screen/`: `interpret_pathways.csv`, `interpret_genes.csv`
and `interpret_handoffs.csv`. Nothing here is a new claim test. These are checks on calls already made, and readouts
of what drives them.

### Pathway calls: four robustness checks

A call counts as **robust** when all four checks pass:
1. The change has the same sign on raw counts of units with that segment's label, in all four pairs.
2. It keeps its sign when the first and last 5% of the PT are trimmed.
3. With its top three genes removed, the set still beats random sets (p ≤ 0.05).
4. It is not carried only by end effects (this is the trimming check, applied to S1 and S3).

95 of the 122 calls pass: S1 77 of 89, S2 12 of 24, S3 6 of 9. They reduce to three signals.

- **S2, replication machinery climbs through S2 in AKI.**
  - Sets: DNA strand elongation +0.42, mismatch repair +0.40, lagging-strand synthesis +0.34, G1/S-specific
    transcription +0.32, mitotic prometaphase +0.15.
  - On raw counts within the S2 label the change is +0.37 to +0.47 log2.
  - Leading genes: Gins1, Mcm6, Mcm7, Pold2, Rfc4, Pcna, Rpa2, Cdt1, E2f1, Tk1, Tyms.
  - **The control mice are flat** (W −0.01 to −0.04), so injury creates this gradient.
- **S3, the respiratory chain falls through S3 in AKI.**
  - Sets: KEGG oxidative phosphorylation −0.22, Reactome respiratory electron transport −0.21.
  - The signal is broad: the top three genes carry 11% of it, and dropping them leaves p 0.001.
  - Raw counts give −0.23, and trimming keeps −0.22.
  - Leading genes: Ndufb3, Ndufb5, Ndufa5, Cox7a2, Atp6v0e2.
  - The controls are flat to slightly rising (+0.04 to +0.09), so injury creates this one too.
- **S1, replication and proteasome sets ease along S1 in AKI. This is weak.**
  - 43 of the 77 robust S1 calls are led by proteasome subunits (Psmc1, Psmb4, Psmf1, Psmd8), and the rest by
    replication genes (Pole2, Gins1, Ccne2, Mcm6).
  - Trimming the first 5% of the PT halves every one of them (median −0.07 log2), so most of the effect sits at the
    very start of the PT.
- **Not robust.**
  - **End effects:** ether lipid metabolism (S3, −0.45 falls to −0.06 when trimmed) and gluconeogenesis (S1).
  - **Two or three genes carry the call:** integrin and platelet sets, cholesterol biosynthesis and SREBP, plasma
    lipoprotein, scavenger receptors, KRAS-down, G alpha (q).
  - So two of the examples used in the second figure version (integrin and cholesterol) were weak.

### Single genes (Category A calls)

The Category A list missed calibration (false-call share 0.12). These calls are therefore named, checked examples,
not a list. Splitting each change into its control part and its AKI part shows two kinds of event.

- **Injury creates a within-segment gradient (controls flat).**
  - S2: Havcr1 (KIM-1) +1.0 in AKI against −0.3 in controls.
  - S3: Flvcr2 −1.8, Acnat1 −2.1, Anpep (CD13) −1.0, Fmo5 −1.1, Esr1 (S1) −0.9.
- **Injury erases a healthy within-segment gradient (AKI flat).** All in S2:

| Gene | Control | AKI |
|---|---|---|
| Cyp4a14 | +1.29 | −0.41 |
| Cyp4a10 | +1.31 | +0.04 |
| Cyp51 | +1.32 | 0.00 |
| Agt | +1.53 | +0.37 |
| Hpd | +1.17 | −0.16 |
| Cyp24a1 | −2.09 | +0.06 |
| Spns3 | −2.72 | −0.46 |
| Steap2 | −0.93 | +0.08 |

  - Hsd17b2 is reduced rather than erased (+2.66 to +1.63).
  - Cyp4a14 and Hsd17b2 are Fig. 4 Category 1 genes, and Cyp24a1 is a Fig. 4 Category 2 example. Injury erases the
    healthy within-PT structure that Fig. 4 shows.
- **Artifacts to drop.** Fgb and Serpina10 (S1) are barely detected in controls (whole-PT about −6 to −7 log2 per 10k).
  Their "rise within S1" is control noise, and the AKI curves are flat.
- **Near the neighbour threshold** (1.8–2.0× in IR): Hpgd, Tcf24, Tcte2, Snx16, Acot11. Acot11 also fails the trim.

### Is erasure systematic? Gradient retention

**Design.**
- Genes are chosen by a within-segment change of at least 0.75 log2 in one control mouse.
- The comparison is how much of it the other control keeps, against each IR mouse.
- Both directions of the choice are run.
- **Why this design.** Choosing in one control and reading the other avoids regression to the mean.

**Result** (median kept fraction, from `interpret_retention.csv`):

| Segment | Held-out control | IR mice |
|---|---|---|
| S1 | 0.43, 0.52 | 0.23–0.32 |
| S2 | 1.16, 0.84 | 0.49–0.56 |
| S3 | 0.58, 0.35 | 0.12–0.32 |

The IR mice keep less in all 12 comparisons.

**Can noisier IR positions explain it?** Noise would flatten every gene equally.
- **S2.** Genes whose whole-PT level barely changes (|Δ| < 0.25) keep 0.66–0.68 in IR against 0.87–0.98 in the held-out
  control. Injury-responsive genes keep 0.45–0.54 against 0.82–1.19.
- **So** about a quarter of the S2 flattening is common to all genes. That is the bound on placement noise, and it
  matches the slightly lower IR placement agreement.
- **The injury-responsive genes lose about twice as much.**
- **S1 and S3.** The two groups do not separate, so there flattening cannot be told apart from noise.
- **Positive control.** IR positions do carry within-segment structure: injury builds new gradients in S2 and S3,
  above.

### Handoffs: what moves them

- **Raw counts and segment means agree** that the upstream shifts come from late, S3-high partners losing expression in
  injured S3:
  - GAT2 (Slc6a13) S3 −2.9 log2;
  - rBAT (Slc3a1) −1.8, which also falls further within S3 on raw counts (−0.26 in controls, −1.08 in IR);
  - PCFT (Slc46a1) −2.1;
  - OCT2 (Slc22a2) −1.8;
  - MCT5 and MCT9 (Slc16a4, Slc16a9) −3.4 and −2.7.
- **MCT1 is the exception.** MCT1 (Slc16a1), an early partner, is induced 27-fold in injured S3 (+4.75 log2).
- **ZIP5 → ZIP14 and ZIP8 move downstream** because ZIP5 stays on further into S2.
- **Drop Slc5a9 → Slc5a3.** Slc5a3 is 23× higher in thick limb than in PT, so it is neighbour signal; the handoff
  screen applied no neighbour filter. The controls also disagree on where Slc5a9 falls. 10 shifted pairs remain.

### Literature check and artifact checks (2026-11-03)

A literature agent read full texts where it could; abstract-only sources are marked. It also checked direction in
the Liu 2017 IRI whole-kidney bulk time course (GSE98622; 48 h IRI against the study's shams, which exist only at 4 h,
24 h and 12 months) and the normal segment levels in Chen 2021 microdissection (GSE150338). Except for Slc7a13 and
Havcr1, the bulk fold changes quoted below are our reanalysis, not results stated by Liu 2017. **Re-checked
2026-11-03 by a full-text reader of every cited paper (`citation_checks/D_results_literature.md`); the rows below
carry its corrections.** Neither dataset is segment-resolved for shape. Both files are in the session scratchpad,
not the repo.

**Correction.** Slc16a4 is **MCT5** (HGNC). MCT4 is Slc16a3. Labels are fixed.

**Artifact checks** (run in this session on S3-labelled units, first against last third of the S3 core):
- **Late-S3 respiratory-chain fall is not thick-limb bleed.**
  - IR S3 units carry less TAL signal than controls (Slc12a1 −3.0 and −2.3 against +0.4 and +1.3 log2 per 10k).
  - TAL signal rises slightly toward late S3 in IR (Umod +0.5 to +0.7), which would raise OXPHOS, yet OXPHOS falls
    (−0.15 to −0.19 on summed counts).
  - The controls' slight OXPHOS rise within S3 may be TAL bleed (Slc12a1 +1.1 to +1.3 toward late S3).
- **MCT1 induction in injured S3 is in PT units.**
  - In IR S3 units, Slc16a1 does not track immune markers (Ptprc, Lyz2, Cd74; ρ −0.04 to 0.10) or endothelial ones
    (Emcn, Pecam1; ρ 0.05 to 0.09).
  - It tracks PT and injury markers weakly (Slc34a1, Lrp2, Havcr1; ρ 0.05 to 0.26).
  - It is detected in 13% and 32% of IR S3 units, so it is patchy.

**Verdicts by finding.**

| Finding | What is known | Positional aspect | Status |
|---|---|---|---|
| Replication and KIM-1 climb through S2 | Proliferating S3 cells after IRI (Humphreys 2008, which counts S3 only and starts at 2 d; Humphreys 2011 pools outer medulla with inner cortex; Witzgall 1994, rat); injury spreads from the corticomedullary junction (Xuanyuan 2025, maps only) | Never measured within a segment; none of these compares S1/S2 with S3 | **New**. Caveat: S-phase transcripts can mean endoreplication. Polyploid cells are mostly cortical S1/S2 (De Chiara 2022), and at 2 d 13.2% of tubular cells are in S/G2/M against 47.1% Ki-67+, about a quarter (Lazzeri 2018). Say "replication program", not "division". |
| Respiratory chain falls through S3 | Respiratory-chain proteins fall after IRI, though PGC-1α protein and cortical ATP rose at 24–72 h (Funk 2012); PGC-1α falls after IRI (Tran 2016, whose OXPHOS panels compare PGC-1α mice, not IRI with sham); zone-level ATP loss (Fujii 2019). Counterpoint: Arnoux 2025's electron-transport gene score is not reduced in PTS1–PTS3 after IRI | No within-S3 gradient reported | **New** (TAL bleed checked, above) |
| Handoffs move upstream as S3 partners fall | PT carrier loss in injured states (Li 2022; Wu 2024; Xuanyuan 2025: healthy-S3 neighbourhood absent at 12 h–2 d); OCT2 and OAT down after IRI; direction confirmed in Liu 2017 bulk (Slc7a13 −4.4, Slc3a1 −1.4, Slc22a2 −1.2, Slc6a13 −1.0, Slc46a1 −1.0) | No study reports where a paralog switch sits or that one moves | **New**. Caveat: dedifferentiation at 12 h–2 d is reversible, so this may be a transient state |
| Cystine handoff keeps its position, AGT1 down about 20× | AGT1 is apical S3 and male-only (Nagamori 2016); b0,+AT is S1-high (earlier work cited there; its own data show b0,+AT equal in both sexes) | Position kept although the level falls | An internal control for the handoff readout |
| MCT1 induced about 27× in injured S3 | MCT1 is basolateral and S1-high (Yanase 2008, abstract only); S3 is normally near-silent (0.2 TPM, Chen 2021); bulk induction about 3–5× in our reanalysis of Liu 2017 and of Balzer 2022's supplementary data (n = 2 per injured group; neither paper mentions Slc16a1); in Yu 2025 the knockout is MCT4 (Slc16a3) and MCT1 is supported only by inhibitors in HK-2 cells | Not reported; about 10× any induction in these reanalyses | Notable but needs ISH or protein; PT-intrinsic by the marker check |
| ZIP5 stays on further into S2 | ZIP8 and ZIP14 mRNA in the outer stripe by ISH (Fujishiro 2012, abstract only); ZIP8 knockout protects in IRI (Cai 2025, which never mentions ZIP14; its S3 panel reanalyses published data) | ZIP5 has no published PT localisation | New, rests on our data alone |
| Cyp24a1 extends across S2 (control S1-restricted, AKI flat and high) | Normally S1-restricted (85.9 / 2.4 / 0.0 TPM, Chen 2021; proximal tubules of the outer cortex by RNAscope, Fuchs 2024); FGF23 rises in AKI; IRI bulk up (our reanalysis of Liu 2017) | Expansion into S2 not reported | **New**. Contradicted in folic-acid AKI (Egli-Spichtig 2018, down 17-fold), so model-specific |
| Agt flattened within S2 | Agt mRNA is in S3 and the outer stripe by ISH, with S1–S2 carrying reabsorbed protein (Matsusaka 2012); by microdissection it is highest in PTS2, the medullary-ray straight tubule (11.7 / 322.1 / 217.1 TPM, Chen 2021); cortical Agt mRNA falls after ischaemia (Kontogiannis 1998: rat, 60-min bilateral, whole-cortex Northern); our Liu 2017 reanalysis −1.0 | Not positional before | Known direction |
| Cyp4a10 and Cyp4a14 flattened within S2 | Yan 2022 shows Cyp4a10/14 lower in a Ptpn6 mutant than in wild type, both after IRI, not IRI against sham | Not positional | Direction after IRI has no source as cited |
| Cyp51 flattened within S2 | Naito 2009 shows ischaemia inducing Hmgcr, not Cyp51; our Liu 2017 reanalysis gives Cyp51 −2.2 at 48 h | — | No conflict on the same gene; direction from our reanalysis only |
| Flvcr2 falls through S3 | Choline transporter; choline depleted after ischaemia (Tran 2016) | New | Interpretable; neighbour ratio 1.7 in IR |
| Fgb, Serpina10, Spp2, Slc34a3, Cyp2d26 in S1; Esr1; Acaa2; Spns3; Steap2; Apob; Hpgd | Fgb marks an injured PT cluster carrying S1/S2 markers (Melo Ferreira 2021); Slc34a3 and Cyp2d26 fall in bulk, against a rise; the others are near the detection floor or unreported | — | **Drop.** Fgb and Serpina10 are control noise (above). |

### The most notable results, in order

1. **Transporter handoffs move upstream after IRI because injured S3 loses its late partners.**
   - 10 of 52 pairs move after the neighbour filter; none moves under either relabeling.
   - Cystine is an internal control: AGT1 falls about 20-fold, yet the handoff stays where it was.
   - New as a positional readout, with direction confirmed in our reanalysis of an independent bulk time course.
2. **Injury builds new gradients inside segments.**
   - The replication program and KIM-1 climb through S2.
   - The respiratory chain falls through S3.
   - The controls are flat for all of them, and both pass all four robustness checks.
   - New as within-segment measurements; injury in these mice is graded within segments, not just S3 against the
     rest.
3. **Injury erases the healthy within-segment structure.**
   - Gradients chosen in one control are kept in the other control but about halved in AKI (12 of 12 comparisons).
   - Injury-responsive genes lose about twice as much as level-unchanged ones, which bound the placement noise.
   - **Examples with literature support:**
     - Cyp24a1 spreads from S1 across S2.
     - Agt is suppressed within S2 (by microdissection it is normally highest in the straight tubule).
     - Cyp4a14 and Cyp4a10 (the PPARα programme) are flattened; the direction after IRI still needs a source.
   - This links to Fig. 4: Cyp4a14, Hsd17b2 and Cyp24a1 are among the genes Fig. 4 shows.
4. **MCT1 is switched on in injured S3, which normally lacks it** (about 27-fold, patchy). Notable, but it needs
   validation by in-situ hybridisation or protein before it can be a headline.

**Figure.** Main Fig. 6 shows results 1–3: a, built gradients (Havcr1, the DNA-replication set, oxidative
phosphorylation, Anpep); b, erased gradients (Cyp4a14, Agt, Cyp51, Cyp24a1); c, gradient retention; d, the handoff rows.
Result 4 goes in the text and Extended Data.
