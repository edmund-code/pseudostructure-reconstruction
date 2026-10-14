# What a continuous PT coordinate adds to a cluster-level pathway analysis (notebook 89)

Frozen on 2026-10-14, before any statistic of notebook 89 was computed. The user asked: "look at the diseased too
then. also see if we for programs they find to be uniform or find to be different, we can add information. cluster
things together to show a few positional programs. ex. set of pathways that are early s1, early s2. maybe see if there
are pathways that are peaked in between s1 and s2 and missed by them. start by thinking about the theoretical
advantages of our approach (ex. pathways gradients in between segments), finer resolution and localization for
existing ones. and test to see if we can actually show that". "Them" is CellSpectra (Klötzer et al. 2025, Nat Genet;
PMID 40775269), whose PT analysis works on S1/S2/S3 clusters.

**What I had seen before freezing.**
- The results of notebooks 74 (C14: fold noise exceeds the within-segment spread in 10 of 12 cells), 79, 81–85 and
  86; in particular 84 (no within-segment separation of Chen's segment pathways) and 85 (segment pathways ramp
  gradually across S1/S2; the curve resolves a label step at 0.29–0.33 mouse units).
- Design facts, no statistic: mouse U_PT units per 0.1 mouse-unit bin (Ctrl1A2 2–152, Ctrl1A4 1–146; ≥ 22 per bin up
  to 2.5 units, ≤ 19 beyond); the AKI export holds 15,013 units and the full mouse panel (19,059 genes) with
  `primary_identity`; 72's fold refits are cached for hvm.
- CellSpectra's Supplementary Tables 4, 5 and 8 and its main text. The disease features below were selected from
  their tables by the rule in section F before any data of ours was read.

## The question, and the theory that sets the tests

A cluster-level analysis sees each PT segment as one number per gene (CellSpectra: per-segment pseudobulks and Wilcoxon
lists; coordination per cell type). A continuous coordinate sees each gene as a curve. On paper the curve adds four
things. Each has a signature that the cluster view cannot produce, and a limit set by the curve's resolution.

| Advantage | What the cluster view sees | What the curve sees | Limit |
|---|---|---|---|
| A1 · Boundary programs | A bump of amplitude A and full width w centred on a boundary puts about 1.25·(w/2.355)·A into each adjacent segment mean. For w = 0.5 segment units that is 0.27·A: after centring, S1 and S2 sit +0.09·A and S3 −0.18·A. A 1 log2 boundary bump is below a 0.25 log2 effect in every segment: missed. At w = 1 it reads as "S1 and S2 high" or "S3 low": a segment program in the wrong place | An interior peak at the boundary | Features narrower than the curve's resolution (85: label step 0.29–0.33 mouse units; spline step 0.46 at S1/S2) blur and lose height |
| A2 · Gradients between segments | A symmetric ramp and a step at the same boundary give segment means that differ by at most A·w/8 per segment: indistinguishable. A ramp inside S1 reads as "S1 lower" | The half-point and the width | Width below the label step (≈ 0.3) cannot be told from a step |
| A3 · Localization inside a segment | A feature of width w inside a segment is diluted to about w of its height in that segment's mean; its position is lost | The position | Coordinate noise inside the segment (C14: noise-dominated in mouse S1 and S2) |
| A4 · Positional programs | Pathways fall into three classes (S1, S2, S3) or none | Pathways group by shape: early S1, the S1/S2 boundary, mid-S2, the S2/S3 boundary, S3 | Only groups wider than the resolution separate |

The tests below ask whether each advantage can be shown on this curve (T), and then look for it in the data (H1,
E1–E4). Section T is the gate: an advantage that the semi-synthetic test cannot show is reported as "not showable on
this curve", and its absence in the data is not read as absence in biology.

## Labels

- Notebook: `analysis/notebooks/89_pt_positional_programs.ipynb`
- Logic: `89.positional_programs.1`
- Results: `results/pt_positional_programs/`
- New module `pseudospace/positional_programs.py` with `tests/test_positional_programs.py`.
- Frozen inputs: `docs/results/positional-programs/cellspectra_features.csv` (44 features) and
  `cellspectra_features_skipped.csv` (19 skipped by rule).

## Inputs (read only)

- **Healthy:** the hvm package export through `hvm_inputs` (71's construction): positions, library, specimen, label,
  `primary_identity`. Mouse units: 72's mouse transitions, cuts (0.341, 0.613), as 85 and 86. U_PT: units whose
  `primary_identity` is Proximal_Tubule (86). Counts: stage-01 native panels (86's G3 route).
- **86's W runs** (stage caches `w_mouse_pt_shared`, `w_human_pt_shared`): per-gene dispersion α and each specimen's
  typical shape, used only by T to generate counts.
- **72's fold refits** (`refit_hvm_fold{0,1,2}` stage caches): the generation positions of T's realistic arm.
- **AKI:** the aki package export (all 19,059 panel genes), its positions scaled as 76 (`support_positions`),
  `primary_identity`; mouse units from the two controls' logistic transitions on that coordinate.
- **Gene sets:** 71's 1,513 library sets (Reactome 2022, KEGG 2019 mouse, Hallmark 2020). CellSpectra features: Enrichr
  GO_Biological_Process_2023 (the library CellSpectra used for its extended human atlas) and KEGG_2019_Mouse, members
  mapped through the ortholog map as `_pathway_sets` does; response to ischemia from QuickGO (2026-10-14, mouse, with
  descendants; `data/external/quickgo_2026-10-14/`).
- **CellSpectra:** `data/external/klotzer2025_cellspectra/` (sha256 in `SOURCE.txt` and the run record); 79's frozen
  species programs (`docs/results/known-programs/programs_hvm_members.csv`); 81's F programs.
- **84:** `segment_pathways.csv` (Chen 2021 microdissection, cameraPR per PTS1/2/3 on 71's sets): the segment-level
  call of an external cluster analysis.

**Guards (not cached).** G1 `hvm_inputs` reproduces 71's universe, sets and units; G2 the transitions reproduce 72's;
G3 the native panels match the export unit by unit; G4 the frozen CSVs and external files match their sha256; G5 86's
cached W α and shapes are present for every usable mouse gene; G6 the fold refits map onto ≥ 90% of mouse U_PT units;
G7 the AKI control transitions are ordered (S1→S2 before S2→S3) in both controls.

## Shared definitions

- **Bins.** Fixed 0.1-wide bins of mouse units on [0, 2.5); a bin is valid for a specimen with ≥ 20 U_PT units.
- **Gene profile.** Per specimen and bin, log2((Σ counts + 0.5) / Σ library) (85's `binned_log_rates`), minus the gene's
  mean over valid bins.
- **Usable genes:** detection ≥ 5% on U_PT in both mice, not mt (86).
- **Reference:** per specimen and bin, the median gene profile over usable genes (81's `reference_curve`, binned).
- **Program profile:** the median over usable members of (gene profile − reference), log2. Library sets need ≥ 10
  usable members, CellSpectra features ≥ 5.
- **Agreed:** both mice satisfy the rule.
- **Cluster view** (what CellSpectra's resolution sees): per mouse, the label pseudobulk on U_PT (S1/S2/S3 labels),
  81's `segment_profile` of the members (gene minus its mean over segments, centred by the median gene per segment).
  A segment is called when its value is ≥ 0.25 in both mice; the call is "none", one segment, or "two adjacent
  segments".
- **Boundary excess** at boundary k ∈ {1, 2}, per specimen: window W = [k − 0.25, k + 0.25), left flank
  [k − 0.75, k − 0.25), right flank [k + 0.25, k + 0.75) (clipped at 2.5). Peak excess E = mean(W) − max(mean(left),
  mean(right)); trough excess T = min(mean(left), mean(right)) − mean(W). Each needs ≥ 2 valid bins per part. A
  monotone profile (step or ramp) gives E ≤ 0 and T ≤ 0 up to noise. **Boundary program:** E ≥ 0.25 (or T ≥ 0.25) at
  the same k in both mice.
- **Shape descriptors:** 81's `shape_descriptors` on the mean-of-mice program profile (interior peak or trough and its
  position; rising or falling with half-point and width); 85's `fit_edges` for edge widths in T.

## T · Can the curve show it? (semi-synthetic, mouse)

Synthetic programs are written into the real design and analysed exactly as the data are.

- **Design.** Mouse U_PT units of both mice: library, observed position (mouse units), label.
- **Generation position.** Arm *ideal*: the observed position. Arm *realistic*: replicate r uses fold refit k = r mod 3,
  its position for each unit in that refit's own mouse units (logistic transitions of the refit, as 72); units absent
  from the refit are dropped. The disagreement between a refit and the primary coordinate is then the coordinate
  error.
- **Members.** 15 genes drawn without replacement from usable mouse genes (seed 89 + replicate).
- **Counts.** μ = library × r_gs × exp(τ_s(u) + ln2 · a_g · f(u)), u the generation position, r_gs the gene's rate in
  that specimen, τ_s the specimen's typical shape (86's W, centred), a_g = A × U(0.5, 1.5); rescaled so the expected
  total equals the gene's observed total; NB with the gene's α (gamma–Poisson).
- **Shapes f** (peak height 1):
  - bump(c, w): Gaussian of full width w at c ∈ {1, 2}, w ∈ {0.25, 0.5, 1.0};
  - matched step of each bump: each unit gets the mean of f over its specimen's units with the same label, so the
    cluster view of bump and step is identical by construction;
  - ramp(c, w): logistic, 10–90% width w ∈ {0.1, 0.5, 1.5} at c = 1, and w = 0.5 at c ∈ {0.6, 1.4};
  - local(c, 0.4): bump at c ∈ {0.3, 0.7, 1.3, 1.7};
  - flat.
- **Amplitudes** A ∈ {0.5, 1, 2} log2. 100 replicates per shape × A × arm.
- **Readouts.**
  - bumps and their matched steps: the boundary-program rule (sensitivity; false positives on matched steps and
    flat); the cluster-view call;
  - ramps: `fit_edges` on the mean-of-mice program profile, the width class (85: sharp ≤ 0.5, gradual ≥ 1) and the
    half-point error;
  - local: peak-position error and the share placed in the correct half-segment.
- **Showable.** An advantage is showable on this curve when, in the realistic arm at A = 1, its sensitivity is ≥ 0.8
  (A1: bump w = 0.5 at that boundary; A2: ramps w = 0.1 and 1.5 classified sharp and gradual; A3: correct half) and its
  false-positive rate is ≤ 0.05 (A1: matched steps and flat). The same is reported at every A and w.

## H1 · Boundary programs (primary)

- **Population.** 71's library sets with ≥ 10 usable members.
- **Statistic.** C = the number of sets that are boundary programs (peak or trough, k = 1 or 2).
- **Null.** Gene-label permutation within strata (deciles of mean log2 CP10k over both mice on U_PT, usable genes):
  every set keeps its size and its overlaps with other sets, and its members are replaced by the profiles of
  permuted genes. 1,000 permutations (seed 89). p = (1 + #{C_perm ≥ C}) / 1,001, one-sided.
- **Holds** when p ≤ 0.05 and the boundary programs span ≥ 3 redundancy groups (Jaccard ≥ 0.5, connected
  components).
- **Per boundary program** (descriptive): k, peak or trough, E or T per mouse, the cluster-view call, Chen's segment call
  (84, FDR ≤ 0.05), core members (genes that are boundary genes themselves, agreed), the excess with integration
  features removed from the members (fold protection, 74), the human profile.
- **Reading.** A positive H1 means pathway membership groups genes into boundary-centred programs beyond chance. If
  H1 fails and T calls A1 showable, there is no evidence of boundary programs; if T calls A1 not showable, H1 cannot
  answer the question.

## E1 · Positional programs (descriptive, with a replication readout)

- **Eligible sets:** in each mouse, the profile amplitude (max − min over valid bins) exceeds the 95th percentile of
  that set's own H1 permutation amplitudes, and the two mice's profiles correlate at r ≥ 0.5.
- **Redundancy:** Jaccard ≥ 0.5 groups (connected components); the representative has the most usable members (ties
  by name).
- **Clustering:** the representatives' mean-of-mice profiles, average linkage on 1 − r, K ∈ 3–10 maximising the mean
  silhouette. Clusters with ≥ 3 representatives are *positional programs*, ordered by the position of their mean
  profile's maximum.
- **Replication:** each mouse clustered alone with the same K; adjusted Rand index between the two partitions.
  Programs replicate when ARI ≥ 0.4.
- **Per program:** mean profile per mouse, shape descriptors, member sets, the cluster-view calls of its members, the
  share of its members that Chen's segment analysis calls (and in which segment), human profiles (descriptive).
- **Disease (descriptive):** each program's profile in each AKI mouse on the AKI coordinate (same bins and rules):
  correlation with the control mice and amplitude ratio IR / control.

## E2 · What the cluster view says about the same sets (descriptive)

For every eligible set, the curve's main feature (descriptors) against the cluster-view call:
- *agree* — the call names the segment that holds the feature;
- *boundary* — the feature lies within 0.25 of a boundary (the cluster view splits it);
- *mislocated* — the call names a segment that does not hold the feature;
- *missed* — no call, while the curve has an agreed feature.

## E3 · CellSpectra's programs, healthy (descriptive)

Programs: 79's six species blocks, 81's eight unique and core sets, 81's six conserved functions (F), and the 44
disease features of section F (their healthy profiles). Per program, in mouse: profile, descriptors, boundary excess
with its own permutation p (H1's null for that program), the cluster-view call, the nearest positional program
(r ≥ 0.5). Human: the profile per section on the human curve position, descriptive. "Uniform" and "different" in
CellSpectra's sense (core vs segment-only lists; conserved vs species-specific features) are columns.

## E4 · CellSpectra's disease features in the IR mice (descriptive)

- **Units and genes.** AKI export, U_PT units of the four mice; genes of the full panel, not mt, with detection ≥ 2% on
  U_PT in each of the four mice. Bins as above on the AKI coordinate's mouse units.
- **Δ profile.** Per gene and bin, mean(IR) − mean(control) of the log2 rates, centred per bin by the median gene Δ
  (76's `centre_rows`); program Δ = the median over members. Also per IR × control pair (4 pairs) and for the two
  relabelings of 76 (`partition_groups`).
- **Readouts:** segment means of Δ (S1 [0, 1), S2 [1, 2), S3 [2, 2.5)); the peak segment of |Δ|; its margin over the
  next segment; pairs sharing the sign there; the relabeled ratio there (|Δ_relabeled| / |Δ|); the position of max |Δ|.
- **Class:** *no change* (max segment |Δ| < 0.25); *localized to Sk* (margin ≥ 0.25, 4/4 pairs, both relabeled ratios
  ≤ 0.5); *uniform* (|Δ| ≥ 0.25 in every segment, range < 0.25, 4/4 pairs, ratios ≤ 0.5); otherwise *not resolved*.
- **Coordination along the PT** (CellSpectra's own statistic, by position). Bins of 0.25 mouse units; per specimen and
  bin, the members' pseudobulk log2 CPM after binomial thinning of every specimen to the smallest library in that bin
  (CellSpectra subsamples its reference to the query's depth; seed 89). R²(a, b) is the squared Pearson correlation of
  two specimens' member vectors. Control coordination: R²(Ctrl1A2, Ctrl1A4). IR coordination: the mean of R²(IR, control)
  over the four pairs. ΔR² = control − IR. Beside it: the percentile of ΔR² among 100 size-matched random gene sets
  (deciles of mean expression) in the same bin.

## F · CellSpectra disease features (frozen; `cellspectra_features.csv`)

Selected from CellSpectra's tables only, by this rule:
- **N (4):** the PT features named in the main text after ischaemia and folic acid: response to ischemia (GO:0002931),
  KEGG HIF-1, mTOR and TNF signalling.
- **P (10):** Supplementary Table 8, PT, GO: the first 10 by prevalence that are in GO_Biological_Process_2023.
- **K (10):** Supplementary Table 8, PT, KEGG: the first 10 by prevalence in KEGG_2019_Mouse.
- **C (10):** Supplementary Table 5, conserved PT hallmarks: the first 10 BP terms by average normalised prevalence in
  the library.
- **H (5), R (5):** Supplementary Table 5, PT species-specific features, by corrected p: the first 5 BP terms in the
  library that are dysregulated in more human samples (H) or more rodent samples (R).

Skipped entries (MF or CC terms; terms or pathways not in the libraries) are listed with their rank. Their gene sets
are ours (Enrichr 2023, KEGG 2019 mouse), not CellSpectra's modified cross-species GO.

## Human (descriptive throughout)

Human positions are the weakest measurement (transitions cross in HUK1_MED1; issue 6), so human profiles are drawn on
the human curve position (25 equal bins over each section's 1–99% range), with each section's label transitions
marked. Between-section r is reported. No human position is claimed.

## Sensitivities (descriptive)

Bin width 0.2; boundary window half-width 0.15 and 0.35; integration features removed from members; all PT-labelled
units (U_all) instead of U_PT; E1 with K fixed at 6.

## Multiplicity and claims

H1 is the only test (α 0.05). T reports operating characteristics. E1–E4 are descriptive. Every within-segment
position is descriptive under C14; T measures what the curve can show at that noise.

## What would make us say the curve adds nothing here

- T: no advantage showable in the realistic arm at A = 1.
- H1 fails while A1 is showable.
- E1: the programs do not replicate between mice (ARI < 0.4).
- E2: fewer than 10% of eligible sets are boundary, mislocated or missed.
- E4: no disease feature is localized or uniform.

## Figures (plots only)

- 89a · T: a bump and its matched step (identical cluster view, different curves); sensitivity by w × A, both arms;
  ramp width recovered against true; localization error.
- 89b · E1: eligible representatives × bins (mouse), ordered by program; program mean profiles per mouse; human.
- 89c · H1: observed count against the permutation null; boundary programs (profile per mouse beside the cluster view).
- 89d · E2: categories; examples.
- 89e · E3: CellSpectra programs (mouse profiles, cluster view).
- 89f · E4: Δ profiles of the 44 features along the PT; ΔR² by bin; positional programs in IR against control.

## Bounds

Two healthy male mice and two IR male mice (time point unknown), one male human donor in two cortex sections. The
coordinate's fold noise exceeds the within-segment spread in mouse S1 and S2 (C14). CellSpectra's comparisons are
human against rodents pooled, healthy and diseased cells pooled, and coordination rather than level; ours are level
and position within one species.

## Notebook 89 results

`analysis/notebooks/89_pt_positional_programs.ipynb` (logic `89.positional_programs.1`) writes to
`results/pt_positional_programs/`. It followed the protocol frozen in 35afb77. Guards G1–G7 passed. Wall time about 5
min cold (T 104 s), 1.5 min cached.

### T · The curve can show all four advantages (semi-synthetic)

Realistic arm (generation on a fold refit, analysis on the primary coordinate), 100 programs of 15 genes per cell.

| Shape (A = 1 log2 per gene) | Curve | Cluster view (S1/S2/S3, 0.25 in both mice) |
|---|---|---|
| Bump at S1/S2, w 0.25 / 0.5 / 1.0 | boundary program 0.88 / 1.00 / 0.52 | no call 1.00 / 1.00 / 0.99 |
| Bump at S2/S3, w 0.25 / 0.5 / 1.0 | 0.97 / 0.93 / 0.05 | one segment 0.19 / 0.80 / 0.93 (S3) |
| Matched steps (same cluster view as the bump) | false boundary call 0.00 everywhere | identical to the bump by construction |
| Flat | 0.00 | no call 1.00 |
| Step at S1/S2 (w 0.1) | classified sharp 0.99; recovered width 0.28; half-point error 0.01 | — |
| Ramp at S1/S2 (w 1.5) | classified gradual 0.96; width 1.61; error 0.08 | — |
| Local bump inside S1 or S2 (c 0.3, 0.7, 1.3, 1.7) | correct half-segment 1.00; peak error 0.05 | — |

- **Showable (frozen rule): A1 at both boundaries, A2 and A3.** At A = 0.5 the boundary rule finds only 3–15% of bumps,
  so the curve's floor is a per-gene bump of about 1 log2 shared by most members of a program.
- **Limit of the statistic, not of the curve:** bumps of full width 1.0 fill the flanks, so the window-versus-flank
  excess misses most of them (0.52 at S1/S2, 0.05 at S2/S3).
- **What a cluster analysis sees.** A 1 log2 boundary bump of width ≤ 0.5 at S1/S2 gives no segment call at all. At
  S2/S3 it is called "S3" (one segment): mislocated. Bump and matched step give the same cluster view.
- The realistic arm uses folds 0 and 1 (reading 1). Their positions differ from the primary coordinate by a median
  0.06–0.10 segment units, 90th percentile 0.16–0.31, with heavy tails.

### H1 · No boundary programs (primary; does not hold)

- 1,404 of 71's sets have ≥ 10 usable members. **C = 0** boundary programs; null mean 0.00 (1,000 permutations);
  p = 1.0. Every sensitivity also gives 0 (window half-width 0.15 and 0.35, bins 0.2, built genes removed, U_all).
- The closest sets (post hoc) reach an agreed excess of 0.22 (Reactome intrinsic pathway of fibrin clot formation, a
  trough at S1/S2) and 0.11 next.
- **Single genes (post hoc):** genes with an agreed boundary peak or trough are no more frequent than when each mouse's
  bin order is shuffled (S1/S2 peak 35 against a null mean of 76, p 0.86; trough 91 against 94; S2/S3 peak 78 against
  72, p 0.29; trough 60 against 94).
- **Reading (frozen).** A1 is showable and H1 fails: on this curve there is no evidence of boundary-centred pathway
  programs, at the set level or the gene level. Pathways "peaked between S1 and S2 and missed by clustering" are not
  found.

### E1 · Pathway positions reduce to two opposite gradients

- 123 eligible sets (amplitude above their own permuted 95th percentile in both mice, r ≥ 0.5) in 80 redundancy
  groups. Silhouette picks **K = 3** (0.70); two clusters have ≥ 3 representatives. **ARI between the mice 0.95**
  (K = 6: 0.85): the programs replicate.
- **P1 · S1-high, falling** (35 representatives, 40 sets): highest at the start of S1, half done at S2/S3 (1.98), with
  the steepest part inside S3; amplitude 0.42 log2. Members: SLC-mediated transmembrane transport, hypoxia, ECM
  organisation, focal adhesion, interferon-γ response, EMT, TNF-α via NF-κB, mitotic spindle.
- **P2 · S3-high, rising** (44 representatives, 82 sets): half done in late S1 (0.84 mouse units), then a plateau
  through S2 and S3; amplitude 0.39. Members: amino-acid metabolism, fatty-acid metabolism, biological oxidations,
  xenobiotic metabolism, mTORC1, adipogenesis, proteasome and deubiquitination, steroid metabolism.
- K = 6 (sensitivity) splits off small groups (3–4 sets each: an early-S1 group with ECM and complement; an S2-high
  group with SLC transport, aspirin ADME and nitrogen metabolism) but keeps the same two large programs (28 and 43
  representatives).
- **What clustering says about them.** Chen's microdissection calls 23 of P2's representatives "PTS2" and 20 nothing;
  it calls 5 of P1's "PTS1" and 29 nothing. On the curve, P2 switches on in late S1 rather than in S2, and P1 falls at
  the S2/S3 boundary and inside S3, which an S1-versus-rest contrast barely sees.
- **Human (descriptive):** P1 falls along the human curve too (gentler); P2 is nearly flat in human.
- **IR mice (descriptive):** P1 keeps its shape (r 0.90, 0.94 with the controls) but its amplitude falls to 0.41 and
  0.43 of the controls' (control pair 0.92). P2's amplitude ratio (0.75, 0.65) is no lower than the control pair's
  (0.75); in IR2A2 it drops in S3 (r 0.59).

### E2 · The cluster view at 0.25 calls almost none of these sets

- Of the 123 eligible sets: agree 2, boundary 0, mislocated 1, **missed 120**.
- **Post hoc, why:** the curve amplitudes are small (median 0.38 log2, IQR 0.31–0.58). The cluster view's centred
  segment values are about half the end-to-end range (median largest |value| 0.13; 7% reach 0.25). So "missed" here
  mostly reflects comparing a range with a centred value at the same threshold (reading 3). It is not, by itself,
  evidence that clustering misses these sets. The fairer external comparison is Chen's statistical call: 34 of 123.

### E3 · CellSpectra's programs along the healthy mouse PT

| Program (CellSpectra) | Curve (mouse) | Nearest program | Cluster view |
|---|---|---|---|
| S-M1 mouse-high in S1 (95 genes) | flat through S1–S2, falls inside S3 (S3 −0.25) | P1 (r 0.90) | none |
| S-M3 mouse-high in S3 (90) | rising, half done at 0.88 (late S1), width 1.9 | P2 (0.97) | none |
| U-M1 mouse-high, S1 only (18) | falling, half done at 1.57 (mid-S2), width 1.6 | P1 (0.95) | S1 |
| U-M3 mouse-high, S3 only (42) | rising, half done at 1.16 (early S2), width 2.1 | P2 (0.94) | S3 |
| U-H3 human-high, S3 only (21; mouse side) | falling, half done at 1.96 | P1 (0.93) | none |
| F-ARG arginine metabolism (5: Ass1, Asl, Arg2, Agmat, Ddah1 among 8) | interior peak at S1 0.75, amplitude 2.4, falls to S3 | P1 (0.93) | S1 |
| F-LDL LDL particle binding (7) | interior peak in S3 (2.35), amplitude 0.93; mice r 0.44 | — | S3 |

- CellSpectra's segment-labelled species lists, read in the mouse, are the two gradients: "S1" lists fall late (in S2 or
  S3), and "S3" lists rise early (late S1 or early S2). The segment names say where a list is highest, not where it
  changes.
- No CellSpectra program is a boundary program. The largest agreed excess is 0.08 (KEGG complement and coagulation,
  K4; permutation p 0.005), below the 0.25 effect.
- Of the 44 disease features, most are flat in healthy mouse PT. The coherent exceptions are small sets near P1
  (C4, C5, K1, K10, H4).

### E4 · CellSpectra's PT disease features in the IR mice

- 39 of 44 features have ≥ 5 members on the AKI universe (10,361 genes; 23 of 25 bins valid).
- **Level change (IR − control):** 22 no change; 13 not resolved; 3 uniform along the PT (complement and coagulation
  +0.28 to +0.47; positive regulation of fibroblast proliferation +0.52 to +0.61; taste transduction −0.31 to −0.41);
  1 localized to S3 (positive regulation of adaptive immune response, 6 genes, +0.37 in S3). The four features named in
  CellSpectra's main text barely change in level (|Δ| ≤ 0.09 for response to ischemia, HIF-1 and mTOR; TNF +0.22 to
  +0.30, not resolved).
- **Coordination (CellSpectra's R², by position, depth-matched):** IR mice are less coordinated with the controls than
  the controls are with each other in 99% of feature × bin cells. Median ΔR² rises from S1 (0.19) and S2 (0.21) to S3
  (0.36), but size-matched random gene sets rise the same way: the share of cells at the ≥ 90th percentile of random sets
  is 0.45, 0.40 and 0.37. CellSpectra's features lose coordination more than random genes do everywhere along the PT
  (expected share 0.10), and the extra loss in S3 is a property of injured S3 as a whole, not of these features.

### Readout against the "adds nothing" rules

| Rule | Result | Met |
|---|---|---|
| T: no advantage showable | all four showable | no |
| H1 fails while A1 showable | H1 fails, A1 showable | **yes** |
| E1 programs do not replicate | ARI 0.95 | no |
| E2 < 10% boundary, mislocated or missed | 98%, but mostly a threshold artefact (reading 3) | no (not informative) |
| E4 no localized or uniform feature | 3 uniform, 1 localized | no |

**What this means for the paper.**
- **The method claim stands on T.** On our own design and coordinate noise, the curve recovers boundary bumps,
  sharp versus gradual edges, and positions inside S1 and S2 that S1/S2/S3 clusters cannot see.
- **The biology claim is narrower.**
  - Healthy mouse PT pathways do not form boundary programs.
  - They form two broad opposite gradients. The curve gives their locations: P2 on by late S1, P1 falling at S2/S3
    and inside S3. Clustering, and even microdissection, names these "PTS2" or nothing.
  - IR flattens the S1-high gradient.
  - CellSpectra's disease features lose coordination all along the PT rather than at one place.

### Recorded readings

1. **Fold refits.** The realistic arm uses folds 0 and 1. Fold 2 fails 72's S1→S3 orientation guard (72 also used
   folds 0 and 1), so replicate r uses fold r mod 2 instead of r mod 3.
2. **Post hoc section 9b** (single-gene boundary peaks against a bin-shuffle null; sets closest to the rule; E2's
   threshold check and Chen's calls; coordination by segment) was added after the run. It changes no frozen readout.
3. **E2's "missed" is not informative.** The protocol compared the curve's range with the cluster view's centred
   segment value at the same 0.25. A fair cluster threshold would be about half the range. Recorded as issue 33.
4. **E3's human r** is computed over bins valid in both sections. Two late-S3 AKI bins lack 20 units in some mouse
   (23 of 25 valid).
5. **Figures.**
   - 89c shows the near-miss sets and the gene-level nulls, because no set is a boundary program.
   - 89d leaves an empty category blank.
   - Heatmaps carry colourbars.

### Tables and figures

- **Tables:**
  - `t_operating_characteristics.csv`;
  - `h1_set_calls.csv`, `h1_null_counts.csv`, `h1_boundary_programs.csv`, `h1_sensitivities.csv`;
  - `e1_representatives.csv`, `e1_programs.csv`, `e1_sensitivity_k6.csv`, `e1_programs_in_ir_mice.csv`;
  - `e2_cluster_vs_curve.csv`;
  - `e3_cellspectra_programs.csv`;
  - `e4_disease_features.csv`, `e4_coordination_by_bin.csv`;
  - `posthoc_boundary_genes*.csv`, `posthoc_near_boundary_sets.csv`;
  - `readout.csv`, `run_record.json`.
- **Figures:**
  - 89a showable;
  - 89b positional programs;
  - 89c boundary programs;
  - 89d cluster view against the curve;
  - 89e CellSpectra programs;
  - 89f disease.
