# Experiment 0005 - Compact milestone replies: verdict and evidence, tables only when asked

- **id**: `reply-style-compact`
- **scope**: auto (applied automatically)
- **status**: superseded — replaced by `caveman-reply-style` (experiment `0008-caveman-reply-style.md`)
- **opened**: 2026-09-21T14:23:28+00:00
- **superseded**: 2026-09-21
- **setting**: `reply_style = compact`

## Why

output + reasoning is the dominant billable term (764k tokens in the last session) and the burst events are milestone narratives

## Prediction (must be beaten, not asserted)

output tokens per request -30%

## Measurement basis

Last recorded work order: `notebook6-umap-diagnosis` - 0 rounds, 0.0k output, 0.0k cache-miss, $0.0000; gates PASS.

## Verification

Superseded before its deciding measurement landed, so this experiment has no verdict and its
prediction was never tested. The reply-style variable passed to `caveman-reply-style`, which carries
the measurement.

Do not run `python tools/agent_loop.py revert reply-style-compact` here: `revert` restores
`reply_style` to the baseline value (`narrative`), which would silently undo the replacement rather
than retire this experiment. Revert `caveman-reply-style` instead if the compressed reply style
regresses.
