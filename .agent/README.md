# `.agent/` — the efficiency loop's state

Small, tracked, revertible. Nothing here writes to the harness's own files.

| File | What it is |
| --- | --- |
| `policy.toml` | Who may change what: behaviour-scoped changes apply automatically, standing instructions and harness config need approval |
| `champion.json` | The settings in force, the baseline measurement, what has been applied or retired |
| `log.jsonl` | One row per work order: rounds, output tokens, cache-miss tokens, cost, gate verdict |
| `gates.json` | The last gate run (written by `tools/agent_gates.py`) |
| `experiments/` | One file per candidate change: why, the prediction it must beat, and how it is verified |

## The loop, in four commands

```bash
git add -A <changed paths>                      # hygiene reads the index, so stage first
python tools/agent_gates.py                     # hard constraint: is this safe to hand back and commit?
python tools/agent_metrics.py --since <ISO>     # what did that work order cost?
python tools/agent_loop.py log "<work order>" --since <ISO>   # record it
python tools/agent_loop.py propose              # the single next candidate change
```

## Why the metrics are what they are

The harness telemetry (`~/.reasonix/stats/*.jsonl`) shows the re-sent conversation is ~99.8%
**cache hits**, so context size is nearly free. What is billed is the new content each round
(cache **miss**) and the **output + reasoning** stream. That is why the loop optimises rows in this
order: rounds and gates first (quality per round), output discipline second, cost third.

## Rules that keep it honest

- A change is accepted only if its prediction is beaten **and** every gate stays green.
- One change per experiment; a regression is reverted, never rationalised.
- Measurements come from new work orders, not the ones that inspired the change.
- No gate may be weakened to make a change look cheap.
