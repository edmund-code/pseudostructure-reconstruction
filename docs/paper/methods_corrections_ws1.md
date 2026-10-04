# Corrections to `docs/paper/methods_draft.md` from workstream 1

Line numbers refer to the file as read on 2026-10-04. Evidence for each correction is in
notebooks 35–36, `results/pt_reconstruction_v2/` and `SCRATCH/ws1/plan.md`.

1. **Pass-2 PCA used 611 genes, not 752** (lines 263–265).
   - Replace "species-intersected HVGs (752 genes), 50 PCs" with: "species-intersected HVGs
     (752 genes). Principal components (50) were computed on the 611 of these that were also among
     pass 1's 2,000 Seurat HVGs."
   - Cause: scanpy 1.11's `sc.tl.pca` defaults to `mask_var='highly_variable'`, and the pass-2
     object inherits pass 1's `highly_variable` column.
   - Pass 1 is unaffected. Its `highly_variable` column is created after PCA, so the pass-1 PCA used
     all 732 genes, as stated on line 215.

2. **The coordinate section is no longer a placeholder** (lines 281–284). Replace the [WS1] box
   with the decision:
   - Notebook 13's scFates coordinate is retained.
   - Alternatives were evaluated in notebook 36 under pre-specified rules: equal-depth refit, DPT,
     PT-only integration and a conserved-anchor coordinate.
   - Every candidate also had to keep the `T_spatial` screen specific. This gate was applied to all
     candidates after DPT's result, as a coordinator clarification. Say so.
   - Downstream counts that depend on the coordinate are unchanged.

3. **Input space** (line 289). Add "pass-2 Harmony dimensions computed from 611 PCA genes (see
   above)".

4. **Orientation** (lines 292–297). Add two points:
   - Two of the six orientation markers are not zonated in human external data: Gatm (human snRNA
     S3 − early = −0.10) and Cyp7b1 (+0.39). The marker axis therefore sets only the root and is not
     a cross-species check.
   - scFates' automatic tip rule sets interior nodes to 0, so an all-negative tip score would root
     the curve at an interior node. The per-species z-scored score avoids this here. Notebook 36
     guards it explicitly.

5. **Projection uncertainty** (line 299). After "`n_map = 1`" add:
   - "Twenty resampled mappings changed positions by SD ≈ 0.007 and gene `T_spatial` by ρ = 0.994."
   - "This is far smaller than the disagreement between coordinates built from disjoint gene thirds,
     so `n_map` understates positional uncertainty."

6. **Add a subsection "Coordinate validation"** after the DPT comparator (around line 318). Use the
   Methods paragraph in `SCRATCH/ws1/r1_text.md`:
   - fold-protected external concordance (P1/P1w);
   - injected species × position absorption;
   - equal-depth thinning;
   - count splitting;
   - registration and numerical stability;
   - the arm comparison.

7. **Numerical reproducibility** (Software section, around line 734). Add:
   - "The scFates curve is numerically sensitive. Multiplying by 1e4/library instead of dividing
     changes values by ≤ 9 × 10⁻¹⁶, yet brings median within-segment agreement down to 0.985 and
     shifts individual structures by up to 0.44."
   - "Unsorted CSR indices alone, with identical values, give an overall Spearman of 0.989 and a
     maximum shift of 0.48."
   - "Exact reproduction requires scanpy's normalisation and sorted CSR indices. Notebooks 35–36
     verify it before every use. DPT is unaffected."

8. **External zonation is label-based** (lines 654–658).
   - State that the within-species zonation comparison uses reviewed-segment means, not the
     coordinate.
   - Add: "A coordinate-based version (notebook 35) correlates each gene's within-specimen Spearman
     correlation with the coordinate against the same external contrasts, using a gene-fold-protected
     coordinate. It gives 0.75 / 0.53 for mouse (male snRNA / microdissection) and 0.25 for human
     (snRNA donor half B), against 0.73 / 0.51 / 0.27 from labels."
   - Note that notebook 35 splits external donors between anchor selection and evaluation.

9. **External species × position is label-based except for gene selection** (lines 659–664). Add
   the coordinate-based interaction from notebook 36: r_human − r_mouse against human snRNA half B
   minus male mouse snRNA. Spearman is 0.44 on all shared genes and 0.63 on the top 10% by
   `T_spatial` (SCF13).

10. **Trajectory method** (lines 619–626). Replace the [WS1] open item with notebook 36:
    - DPT on notebook 13's embedding, with the same universe, strata and pathway family: 96 calls,
      56/66 robust and 38/42 core retained, NCDR 0.40 (not specific).
    - The conserved-anchor coordinate: 42/66 robust and 30/42 core retained, on 7,167 structures,
      because it does not register the species.
    - Keep the notebook-03 DPT run as historical, or drop it.

11. **Limitations to add under Statistics** (around line 718). The coordinate's human-S1 ordering
    depends on read depth:
    - Equal-depth agreement is 0.36–0.46 in human S1.
    - The correlation with library size falls from −0.31 / −0.38 to −0.04 / −0.13 after thinning.
    - Its human–mouse placement moves by up to 0.29 when a third of the genes or half the reads are
      removed.

12. **Open item, line 257** ("[WS1/coordinator: report their concordance]"). Not addressed by
    workstream 1. Notebooks 35–36 use notebook 13's labels unchanged.
