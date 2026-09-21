"""The efficiency loop: measure a work order, propose one change, apply it, re-measure.

There is no model training here - the weights are remote - so what evolves is the *scaffold*: reply
style, batching, tool-output discipline, reconnaissance routing, and how gates are used. The loop is
deliberately small and honest:

    report   - what the last work orders cost (rounds, output tokens, cache misses, dollars) + gates
    propose  - the single next candidate change, ranked by expected effect, not yet tried
    apply    - behaviour-scoped changes are applied and recorded; anything touching standing
               instructions or harness config is written as an *experiment awaiting approval*
    revert   - a change that did not help is retired, and stays retired
    log      - record the current window's measurement into .agent/log.jsonl

Hard constraints, from .agent/policy.toml: gates must pass, one change per experiment, a regression is
reverted rather than argued away, and a working copy that diverges from HEAD in source is never
regenerated.

Usage:
    python tools/agent_loop.py report
    python tools/agent_loop.py propose
    python tools/agent_loop.py apply batch-rounds
    python tools/agent_loop.py revert reply-style-compact
    python tools/agent_loop.py log "review round 3 fixes" --since 2026-09-17T09:00
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
from datetime import datetime, timezone

REPO = pathlib.Path(__file__).resolve().parent.parent
AGENT_DIR = REPO / ".agent"
CHAMPION_PATH = AGENT_DIR / "champion.json"
EXPERIMENTS_DIR = AGENT_DIR / "experiments"

# Ranked by expected effect on rounds and output tokens, with the prediction each change must beat.
CATALOGUE = [
    {
        "id": "gate-before-handback",
        "title": "Run the gate suite before every hand-back, and only hand back green",
        "scope": "auto",
        "setting": ("handback_gate", "gates-green-required"),
        "why": "the single largest source of wasted rounds in the last session: eight user-visible "
               "failures that the dry-run and column audit would have caught before hand-back",
        "prediction": "first-try gate-pass rate -> 1.0; user re-runs per work order -> 1",
    },
    {
        "id": "batch-rounds",
        "title": "One edit-pack-verify-commit cycle per work order, not per item",
        "scope": "auto",
        "setting": ("rounds_per_work_order", "one-cycle-per-work-order"),
        "why": "each round re-sends the conversation and produces output tokens; eight fixes delivered "
               "one at a time cost eight rounds and seven user re-runs",
        "prediction": "rounds per work order -40%, user re-runs per work order <= 1",
    },
    {
        "id": "quiet-tools",
        "title": "Quiet tooling by default: counts, windows and verdicts, never full dumps",
        "scope": "auto",
        "setting": ("tool_output", "counts-and-windows-only"),
        "why": "cache-miss tokens are the new content each round, and one unbounded grep of a notebook "
               "cost 16k tokens of base64 for no information",
        "prediction": "cache-miss tokens per round -50%",
    },
    {
        "id": "reply-style-compact",
        "title": "Compact milestone replies: verdict and evidence, tables only when asked",
        "scope": "auto",
        "setting": ("reply_style", "compact"),
        "why": "output + reasoning is the dominant billable term (764k tokens in the last session) and "
               "the burst events are milestone narratives",
        "prediction": "output tokens per request -30%",
    },
    {
        "id": "caveman-reply-style",
        "title": "Caveman reply style: compressed output with technical substance intact",
        "scope": "auto",
        "setting": ("reply_style", "caveman-full"),
        "why": "supersedes reply-style-compact before its verdict landed, so the reply-style variable "
               "keeps exactly one owner. Same target - output + reasoning is the dominant billable term "
               "- with a stricter rule set: drop articles, filler, hedging and pleasantries; keep code, "
               "commands, paths, error strings and technical terms byte-exact; forbid tool-call "
               "narration and decorative tables, which is where the milestone bursts came from. Sourced "
               "from the caveman skill (global install, ~/.reasonix/skills/caveman). Honest caveat: the "
               "upstream 65-75% figure is the author's own preliminary benchmark, since repudiated as "
               "not rigorous, and it covers visible output only - never hidden reasoning tokens. "
               "AGENTS.md stays uncompressed: it is cached, so its saving is ~free, and it is the "
               "standing instruction file.",
        "prediction": "output tokens per request -40%; no quality change",
    },
    {
        "id": "subagent-recon",
        "title": "Route wide reconnaissance through the explore/research subagents",
        "scope": "auto",
        "setting": ("recon", "subagent-first"),
        "why": "reconnaissance was the most token-hungry non-build phase: noisy reads enter the main "
               "context, while a subagent returns one distilled answer",
        "prediction": "large-output events per work order -80%",
    },
    {
        "id": "config-routing",
        "title": "Route mechanical rounds to the flash model, escalate review/planning",
        "scope": "approval",
        "setting": ("model_routing", "tiered"),
        "why": "the harness supports it (subagent_model, planner_model, subagent_efforts, "
               "max_output_tokens); only the scientific/architectural rounds need the frontier model",
        "prediction": "cost per work order -30-50%",
        "requires": "config.toml edits: default_model stays flash, review/planning escalate, and an "
                    "output cap is set. Needs your approval because it changes every session.",
    },
    {
        "id": "agents-md-lean",
        "title": "Split AGENTS.md: rules standing, rationale moved to a referenced doc",
        "scope": "approval",
        "setting": ("standing_instructions", "split-rules-from-rationale"),
        "why": "14.3 KB is re-sent every round; it is cached, so the gain is small but real, and the "
               "risk is that a mid-session edit invalidates the cache prefix - so do it between sessions",
        "prediction": "prompt tokens per round -3.5k; no quality change",
        "requires": "a patch to AGENTS.md plus a docs/ rationale file. Needs your approval: it is your "
                    "standing instruction file.",
    },
]


def load_champion() -> dict:
    return json.loads(CHAMPION_PATH.read_text())


def save_champion(champion: dict) -> None:
    champion["updated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    CHAMPION_PATH.write_text(json.dumps(champion, indent=1) + "\n")


def last_measurement() -> dict | None:
    log_path = AGENT_DIR / "log.jsonl"
    if not log_path.exists():
        return None
    rows = [json.loads(line) for line in log_path.read_text().splitlines() if line.strip()]
    return rows[-1] if rows else None


def write_experiment(entry: dict, champion: dict, status: str) -> pathlib.Path:
    EXPERIMENTS_DIR.mkdir(parents=True, exist_ok=True)
    number = len(list(EXPERIMENTS_DIR.glob("*.md"))) + 1
    path = EXPERIMENTS_DIR / f"{number:04d}-{entry['id']}.md"
    basis = last_measurement() or {}
    body = [
        f"# Experiment {number:04d} - {entry['title']}",
        "",
        f"- **id**: `{entry['id']}`",
        f"- **scope**: {entry['scope']} "
        f"({'applied automatically' if entry['scope'] == 'auto' else 'awaiting your approval'})",
        f"- **status**: {status}",
        f"- **opened**: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"- **setting**: `{entry['setting'][0]} = {entry['setting'][1]}`",
        "",
        "## Why",
        "",
        entry["why"],
        "",
        "## Prediction (must be beaten, not asserted)",
        "",
        entry["prediction"],
        "",
        "## Measurement basis",
        "",
        f"Last recorded work order: `{basis.get('work_order', 'none yet')}` - "
        f"{basis.get('requests', '?')} rounds, {int(basis.get('output_tokens', 0))/1e3:.1f}k output, "
        f"{int(basis.get('cache_miss_tokens', 0))/1e3:.1f}k cache-miss, "
        f"${basis.get('cost_usd', 0):.4f}; gates "
        f"{'PASS' if basis.get('gates', {}).get('passed') else 'unrecorded'}.",
        "",
        "## Verification",
        "",
        "The next work order's `.agent/log.jsonl` row decides this. Accept if the prediction holds and "
        "every gate stays green; otherwise run `python tools/agent_loop.py revert "
        f"{entry['id']}` and leave it retired.",
        "",
    ]
    if entry.get("requires"):
        body += ["## Requires your approval", "", entry["requires"], ""]
    path.write_text("\n".join(body))
    return path

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("report")
    sub.add_parser("propose")
    apply_parser = sub.add_parser("apply")
    apply_parser.add_argument("experiment_id")
    revert_parser = sub.add_parser("revert")
    revert_parser.add_argument("experiment_id")
    log_parser = sub.add_parser("log")
    log_parser.add_argument("work_order")
    log_parser.add_argument("--since", required=True, help="ISO start of the work order window")
    args = parser.parse_args()

    champion = load_champion()

    if args.command == "report":
        subprocess.run([sys.executable, "tools/agent_metrics.py", "--since",
                        (last_measurement() or {}).get("window", {}).get("since", "1970-01-01")], cwd=REPO)
        print("\nsettings in force:")
        for key, value in champion["settings"].items():
            print(f"  {key:24s} {value}")
        print(f"\napplied: {', '.join(champion['applied']) or 'none'}")
        print(f"retired: {', '.join(champion['rejected']) or 'none'}")
        print(f"superseded: {', '.join(champion.get('superseded', [])) or 'none'}")
        return 0

    if args.command == "propose":
        done = (set(champion["applied"]) | set(champion["rejected"])
                | set(champion.get("superseded", [])))
        candidate = next((entry for entry in CATALOGUE if entry["id"] not in done), None)
        if candidate is None:
            print("every candidate in the catalogue has been tried; the loop has nothing left to "
                  "propose until a new work order produces new measurements")
            return 0
        print(f"next candidate: {candidate['id']} ({candidate['scope']})\n  {candidate['title']}")
        print(f"  why: {candidate['why']}\n  prediction: {candidate['prediction']}")
        return 0

    if args.command == "log":
        command = [sys.executable, "tools/agent_metrics.py", "--record", args.work_order, "--since", args.since]
        return subprocess.run(command, cwd=REPO).returncode

    entry = next((item for item in CATALOGUE if item["id"] == args.experiment_id), None)
    if entry is None:
        print(f"unknown experiment id {args.experiment_id!r}; catalogue: "
              f"{', '.join(item['id'] for item in CATALOGUE)}")
        return 1

    if args.command == "apply":
        if entry["id"] in champion["applied"]:
            print(f"{entry['id']} is already applied")
            return 0
        key, value = entry["setting"]
        if entry["scope"] == "auto":
            champion["settings"][key] = value
            champion["applied"].append(entry["id"])
            path = write_experiment(entry, champion, "applied (auto scope)")
            save_champion(champion)
            print(f"applied {entry['id']}: {key} = {value}\n  experiment file: {path.relative_to(REPO)}")
            print("  the next work order's measurement decides whether it stays")
        else:
            path = write_experiment(entry, champion, "awaiting your approval")
            print(f"{entry['id']} needs approval before it changes anything.\n"
                  f"  experiment file: {path.relative_to(REPO)}")
            print(f"  requires: {entry.get('requires', '')}")
        return 0

    if args.command == "revert":
        if entry["id"] in champion["applied"]:
            champion["applied"].remove(entry["id"])
        if entry["id"] not in champion["rejected"]:
            champion["rejected"].append(entry["id"])
        champion["settings"][entry["setting"][0]] = champion["baseline_settings"].get(entry["setting"][0], "unset")
        save_champion(champion)
        print(f"reverted {entry['id']}: retired, setting restored to the baseline value")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
