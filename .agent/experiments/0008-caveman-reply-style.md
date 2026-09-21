# Experiment 0008 - Caveman reply style: compressed output with technical substance intact

- **id**: `caveman-reply-style`
- **scope**: auto (applied automatically)
- **status**: applied (auto scope)
- **opened**: 2026-09-21T19:35:21+00:00
- **setting**: `reply_style = caveman-full`

## Why

supersedes reply-style-compact before its verdict landed, so the reply-style variable keeps exactly one owner. Same target - output + reasoning is the dominant billable term - with a stricter rule set: drop articles, filler, hedging and pleasantries; keep code, commands, paths, error strings and technical terms byte-exact; forbid tool-call narration and decorative tables, which is where the milestone bursts came from. Sourced from the caveman skill (global install, ~/.reasonix/skills/caveman). Honest caveat: the upstream 65-75% figure is the author's own preliminary benchmark, since repudiated as not rigorous, and it covers visible output only - never hidden reasoning tokens. AGENTS.md stays uncompressed: it is cached, so its saving is ~free, and it is the standing instruction file.

## Prediction (must be beaten, not asserted)

output tokens per request -40%; no quality change

## Measurement basis

Last recorded work order: `simplify notebook 6 gene universe` - 0 rounds, 0.0k output, 0.0k cache-miss, $0.0000; gates PASS.

## Verification

The next work order's `.agent/log.jsonl` row decides this. Accept if the prediction holds and every gate stays green; otherwise run `python tools/agent_loop.py revert caveman-reply-style` and leave it retired.
