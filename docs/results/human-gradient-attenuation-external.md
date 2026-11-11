# Is the flatter human PT gradient biology? External check (post hoc, 2026-11-11)

Not for the paper text; a check of where the human/mouse amplitude difference (0.63x, `human-arm-attenuation`) comes
from. Script: `analysis/scripts/attenuation_external_check.py`. Measure: per-gene S3 − S1 log2 fold change; amplitude
ratio = SD of human folds / SD of mouse folds over shared, well-expressed genes, and the median per-gene |human| /
|mouse| fold on 193–227 genes zonated (|S3 − S1| ≥ 2) in an independent mouse source (Chen 2021 microdissection), which
avoids selecting on either compared dataset.

| Comparison | SD ratio | Median per-gene ratio (Chen-zonated genes) |
|---|---|---|
| Ours: human sections vs mouse sections | 0.50 | 0.31 |
| Lake 2023 human snRNA (cortex S3) vs census mouse snRNA | 0.48–0.52 | 0.26–0.28 |
| Lake human, medullary S3 vs census mouse | 0.73 | 0.38 |
| Cheval 2012 SAGE microdissection, human vs mouse (≥1 / ≥2 / ≥5 tags per 10k) | 0.91 / 0.84 / 0.76 | too few genes |
| Calibration: ours mouse vs census mouse snRNA | 0.89 | 0.87 |
| Calibration: ours human vs Lake human cortex | 0.90 | 0.85 |
| Calibration: census mouse snRNA vs Chen mouse microdissection | 0.69 | 0.63 |

Reading:
- Our human sections are about as steep as an independent human atlas (0.85–0.90), just as our mouse sections are about
  as steep as an independent mouse atlas (0.87–0.89). The flattening is not specific to our human measurement.
- An independent human atlas is about half as steep as an independent mouse atlas (0.48–0.52), the same ratio as ours.
- Medullary S3 in Lake is steeper than cortical S3 (0.26 → 0.38 on mouse-zonated genes), so cortex-only tissue misses
  part of the late PT, but human stays far flatter than mouse.
- Cheval's single-method microdissection is the only same-lab, same-method comparison; it is too sparse (pooled
  libraries, low tag counts) to decide, and noise pushes its SD ratio toward 1.

What remains confounded in every external pair: human donors are older (50–74 y) and come from nephrectomy or biopsy
tissue, mice are young; human S1/S3 labels come from clustering. A true species difference is the simplest reading
consistent with all three human sources.
