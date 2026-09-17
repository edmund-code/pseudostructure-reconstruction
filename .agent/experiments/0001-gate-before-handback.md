# Experiment 0001 - Run the gate suite before every hand-back, and only hand back green

- **id**: `gate-before-handback`
- **scope**: auto (applied automatically)
- **status**: applied (auto scope)
- **opened**: 2026-09-17T18:32:06+00:00
- **setting**: `handback_gate = gates-green-required`

## Why

the single largest source of wasted rounds in the last session: eight user-visible failures that the dry-run and column audit would have caught before hand-back

## Prediction (must be beaten, not asserted)

first-try gate-pass rate -> 1.0; user re-runs per work order -> 1

## Measurement basis

Last recorded work order: `none yet` - ? rounds, 0.0k output, 0.0k cache-miss, $0.0000; gates unrecorded.

## Verification

The next work order's `.agent/log.jsonl` row decides this. Accept if the prediction holds and every gate stays green; otherwise run `python tools/agent_loop.py revert gate-before-handback` and leave it retired.
