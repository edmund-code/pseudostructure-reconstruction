"""Measure what a work order cost: rounds, output tokens, cache misses, dollars.

Reads the harness's own telemetry **read-only** (`~/.reasonix/stats/*.jsonl`, override with
`REASONIX_STATS_DIR`) and joins it to this repository's work-order log (`.agent/log.jsonl`). It writes
nothing outside `.agent/`, because the harness owns its own state files.

Why these columns: the telemetry shows the re-sent context is ~99.8% cache hits and therefore nearly
free, while the billable mass is (a) cache *misses* - the new content each round: tool output, the
user's message, the reply - and (b) output + reasoning tokens. So the metrics that matter are
per-request miss/output volume, the number of requests (rounds), and the cost those imply.

Usage:
    python tools/agent_metrics.py                      # last 7 days, per day
    python tools/agent_metrics.py --since 2026-09-16T22:00
    python tools/agent_metrics.py --record "review round 3 fixes" --since 2026-09-17T09:00
    python tools/agent_metrics.py --json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pathlib
import statistics
from datetime import datetime, timezone

REPO = pathlib.Path(__file__).resolve().parent.parent
AGENT_DIR = REPO / ".agent"
LOG_PATH = AGENT_DIR / "log.jsonl"


def stats_dir() -> pathlib.Path:
    return pathlib.Path(os.environ.get("REASONIX_STATS_DIR", pathlib.Path.home() / ".reasonix" / "stats"))


def load_records() -> list[dict]:
    """Every telemetry record on disk, newest file last. Unreadable lines are skipped, not fatal."""
    records = []
    for path in sorted(glob.glob(str(stats_dir() / "*.jsonl"))):
        for line in open(path, encoding="utf-8", errors="replace"):
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            record.setdefault("ts", "")
            records.append(record)
    return records


def _number(record: dict, key: str) -> float:
    try:
        return float(record.get(key) or 0)
    except (TypeError, ValueError):
        return 0.0


def summarise(records: list[dict]) -> dict:
    """The cost picture for a set of requests."""
    requests = len(records)
    output = sum(_number(r, "completion") for r in records)
    reasoning = sum(_number(r, "reasoning") for r in records)
    miss = sum(_number(r, "cache_miss") for r in records)
    hit = sum(_number(r, "cache_hit") for r in records)
    cost = sum(_number(r, "cost_amount") for r in records)
    return {
        "requests": requests,
        "output_tokens": int(output + reasoning),
        "cache_miss_tokens": int(miss),
        "cache_hit_tokens": int(hit),
        "cost_usd": round(cost, 5),
        "cost_per_request_usd": round(cost / requests, 5) if requests else 0.0,
        "miss_share": round(miss / max(hit + miss, 1), 5),
    }


def busiest(records: list[dict], count: int = 5) -> list[dict]:
    """The requests that actually cost money: biggest output bursts."""
    ranked = sorted(records, key=lambda r: -(_number(r, "completion") + _number(r, "reasoning")))
    return [{
        "ts": str(record.get("ts", ""))[:19],
        "model": record.get("model"),
        "output_tokens": int(_number(record, "completion") + _number(record, "reasoning")),
        "cache_miss_tokens": int(_number(record, "cache_miss")),
        "cost_usd": round(_number(record, "cost_amount"), 5),
    } for record in ranked[:count]]


def read_log() -> list[dict]:
    if not LOG_PATH.exists():
        return []
    rows = []
    for line in LOG_PATH.read_text().splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return rows


def append_log(work_order: str, since: str, until: str, summary: dict, gates: dict | None) -> dict:
    AGENT_DIR.mkdir(parents=True, exist_ok=True)
    row = {
        "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "work_order": work_order,
        "window": {"since": since, "until": until},
        "gates": gates or {},
        **summary,
    }
    with open(LOG_PATH, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--since", help="ISO timestamp; only requests at or after it are counted")
    parser.add_argument("--until", help="ISO timestamp; only requests before it are counted")
    parser.add_argument("--record", metavar="WORK_ORDER", help="append the measured window to .agent/log.jsonl")
    parser.add_argument("--json", action="store_true", help="print the summary as JSON")
    parser.add_argument("--history", action="store_true", help="show the recorded work-order history")
    args = parser.parse_args()

    if args.history:
        rows = read_log()
        if not rows:
            print("no work orders recorded yet (.agent/log.jsonl is empty)")
            return 0
        print(f"{'work order':38s} {'reqs':>5s} {'output(k)':>9s} {'miss(k)':>8s} {'cost($)':>8s}  rounds/order")
        for row in rows:
            print(f"{str(row.get('work_order'))[:38]:38s} {row.get('requests', 0):5d} "
                  f"{row.get('output_tokens', 0)/1e3:9.1f} {row.get('cache_miss_tokens', 0)/1e3:8.1f} "
                  f"{row.get('cost_usd', 0):8.4f}  {'PASS' if row.get('gates', {}).get('passed') else '-'}")
        return 0

    records = load_records()
    if not records:
        print(f"no telemetry under {stats_dir()} - set REASONIX_STATS_DIR if it lives elsewhere")
        return 1
    if args.since:
        records = [r for r in records if str(r.get("ts", "")) >= args.since]
    if args.until:
        records = [r for r in records if str(r.get("ts", "")) < args.until]

    summary = summarise(records)
    gates_path = AGENT_DIR / "gates.json"
    gates = json.loads(gates_path.read_text()) if gates_path.exists() else None

    if args.json:
        print(json.dumps({"summary": summary, "busiest": busiest(records), "gates": gates}, indent=1))
    else:
        print(f"window: {args.since or 'all'} -> {args.until or 'now'}   ({len(records)} requests)")
        print(f"  rounds (requests)      : {summary['requests']}")
        print(f"  output + reasoning     : {summary['output_tokens']/1e3:.1f}k tokens")
        print(f"  cache miss (new content): {summary['cache_miss_tokens']/1e3:.1f}k tokens "
              f"({100*summary['miss_share']:.2f}% of prompt volume)")
        print(f"  cost                   : ${summary['cost_usd']:.4f} "
              f"(${summary['cost_per_request_usd']:.5f}/request)")
        if gates:
            verdict = "PASS" if gates.get("passed") else "FAIL"
            print(f"  gates                  : {verdict} "
                  f"({gates.get('passed_count')}/{gates.get('check_count')} checks, "
                  f"{gates.get('ran_utc', '')[:19]})")
        print("  biggest output bursts (the actual cost driver):")
        for entry in busiest(records, 3):
            print(f"    {entry['ts']}  out {entry['output_tokens']/1e3:5.1f}k  "
                  f"miss {entry['cache_miss_tokens']/1e3:4.1f}k  ${entry['cost_usd']:.4f}")

    if args.record:
        until = args.until or datetime.now(timezone.utc).isoformat(timespec="seconds")
        row = append_log(args.record, args.since or "all", until, summary, gates)
        print(f"\nrecorded work order '{args.record}' in {LOG_PATH.relative_to(REPO)} "
              f"(rounds {row['requests']}, gates {'PASS' if (gates or {}).get('passed') else 'n/a'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
