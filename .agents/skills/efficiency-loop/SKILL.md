---
name: efficiency-loop
description: Run the repository's agent efficiency loop - gate the work before hand-back, measure what the work order cost, then propose or apply exactly one behaviour change. Use at the start of a work order, before every hand-back, and at the end of a work order.
---

# Efficiency loop

The loop optimises three things **together**, with gates as a hard constraint: rounds per work order,
gate-pass rate on the first try, and cost. Telemetry says context is ~99.8% cached and nearly free, so
the billable terms are output tokens and cache misses, and the thing that actually costs the user is
rounds.

## Before hand-back (never skip)

```bash
git add -A <the paths you changed>   # the hygiene gate reads the index, so stage first
python tools/agent_gates.py
```

Hand back only when it prints `GATES: PASS`. A `WARN` on *working-copy revision* means a tracked
notebook's working copy differs from `HEAD` in its **source** — that copy was edited by hand. Never
regenerate such a file (a pack would discard that work); reconcile first, and ask before overwriting.

## At the end of a work order

```bash
python tools/agent_metrics.py --since <ISO start of the work order>
python tools/agent_loop.py log "<short work order name>" --since <ISO start>
python tools/agent_loop.py propose
```

Then, for the proposed candidate:

```bash
python tools/agent_loop.py apply <id>     # auto scope: applied and recorded
```

An `approval`-scope candidate is only ever written as an experiment awaiting the user's decision.

## Rules

- **Gates are the constraint, never the variable.** Cheapness bought by weakening a gate is a loss.
- **One change per experiment.** The next work order's row in `.agent/log.jsonl` decides it.
- **Revert on regression**: `python tools/agent_loop.py revert <id>` retires it permanently.
- **Measure on new work orders**, not the ones that inspired the change.
- **One cycle per work order**: batch edits, then pack/verify/commit once. Per-item hand-backs are the
  single largest avoidable cost.
