# How published continuous-axis and Fan-lab methods papers write Results and draw figures

Compiled 2026-11-11 to restyle Fig. 5 and its Results text. Full texts were fetched from Europe PMC `fullTextXML` or,
where that failed, the PMC article HTML. Results sections and figure legends were read in full for the papers marked
FULL; the others were checked for section structure and legends only.

| Paper | Source | Read |
|---|---|---|
| STalign, Clifton 2023 *Nat Commun* (PMC10709594) | Europe PMC XML | FULL (Results, legends) |
| STdeconvolve, Miller 2022 *Nat Commun* (PMC9055051) | Europe PMC XML | FULL (Results) |
| STcompare, Clifton 2026 *Bioinformatics* (PMC13577157) | Europe PMC XML | FULL (Results, legends, Fig. 2 image) |
| GASTON, Chitra 2025 *Nat Methods* (PMC12370209) | PMC HTML, author manuscript | FULL (Results, legends, Fig. 3 image) |
| Genes2Genes, Sumanaweera 2025 *Nat Methods* (PMC11725504) | Europe PMC XML | headings and legends, Fig. 5 image |
| Halpern 2017 *Nature* liver zonation (PMC5321580) | Europe PMC XML | Fig. 2 image and legend |
| MERINGUE, Miller 2021 *Genome Res* (PMC8494224) | Europe PMC XML | FULL (Results) |
| tradeSeq, Van den Berge 2020 *Nat Commun* (PMC7058077) | Europe PMC XML | FULL (Results) |
| condiments, Roux de Bézieux 2024 *Nat Commun* (PMC10821945) | Europe PMC XML | FULL (workflow, simulations, TGF-β case) |
| Belayer, Ma 2022 *Cell Syst* (PMC9814896) | PMC HTML, author manuscript | FULL (Results to DLPFC genes) |
| CellRank, Lange 2022 *Nat Methods* (PMC8828480); Seurat spatial, Satija 2015 *Nat Biotechnol* (PMC4430369) | Europe PMC XML | headings and openings only (less related) |

## Text

- Headings name what the method enables: "STcompare enables discovery of spatially differential gene expression in
  acute kidney injury", "GASTON uncovers continuous and discontinuous variation". No bold run-in sub-headings and no
  bullet lists inside Results.
- Paragraph pattern: context sentence about the data or tissue, then "We therefore/thus applied X to …",
  "Specifically, …", then the measurement, then a closing "Thus/In this manner/Taken together, …" sentence.
- Expectation before result: "Because these are biological replicates, we hypothesized that …"; "For a good
  alignment, we would expect …"; "Indeed, we found …"; "However, …".
- Controls are introduced as questions: "To evaluate whether …, we …"; comparisons with existing approaches are run in
  Results ("To evaluate whether these methods can address the same tasks …").
- Numbers as "93% (384/415)", named tests with exact P values.
- Genes are italic symbols; the protein or function is named once in a relative clause.
- Biology is hedged: "We speculate …", "may suggest", "consistent with previous studies".

## Figures

- Panel letters only; no section banners and no claim-sentence titles (the claim lives in the legend title).
- Real tissue first (expression on the section), then expression against the learned coordinate with binned points
  and the fitted function, domain boundaries as dashed lines labelled by domain above (GASTON Fig. 3).
- Several genes as small multiples on shared, normalized axes (Halpern Fig. 2b; Genes2Genes Fig. 5).
- Benchmarks against existing approaches as one simple quantitative panel (GASTON AUPRC bars; STcompare volcano).
- Legends: bold title sentence, then one descriptive phrase per panel with n, normalization and statistics.

## Edits these papers prompted (2026-11-11)

- condiments separates where cells sit along a shared trajectory (differential progression) from how genes change
  along it, and argues for one common trajectory when conditions share topology. Added a sentence on where human PT
  tubules sit along the shared axis (they span 0.07–0.98, but half lie in the mouse-S1 range, against a third of mouse
  tubules), pointing to the coverage panel in Extended Data.
- Belayer flags when its comparison is generous to itself; tradeSeq says signal "might be diluted in cluster-based
  DE". Added that straight lines through three cluster means are one simple baseline and finer clustering may recover
  more of each shift.
- MERINGUE and Genes2Genes state the expected result before reporting it ("As expected", "Indeed"); the text already
  does this for the swap control and the clustering comparison.
