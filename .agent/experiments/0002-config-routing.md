# Experiment 0002 - Route mechanical rounds to the flash model, escalate review/planning

- **id**: `config-routing`
- **scope**: approval (awaiting your approval)
- **status**: awaiting your approval
- **opened**: 2026-09-17T18:32:06+00:00
- **setting**: `model_routing = tiered`

## Why

the harness supports it (subagent_model, planner_model, subagent_efforts, max_output_tokens); only the scientific/architectural rounds need the frontier model

## Prediction (must be beaten, not asserted)

cost per work order -30-50%

## Measurement basis

Last recorded work order: `none yet` - ? rounds, 0.0k output, 0.0k cache-miss, $0.0000; gates unrecorded.

## Verification

The next work order's `.agent/log.jsonl` row decides this. Accept if the prediction holds and every gate stays green; otherwise run `python tools/agent_loop.py revert config-routing` and leave it retired.

## Requires your approval

config.toml edits: default_model stays flash, review/planning escalate, and an output cap is set. Needs your approval because it changes every session.
