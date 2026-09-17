"""One command, seconds, deterministic: is the working tree safe to hand back and safe to commit?

This is the loop's *hard constraint*. A candidate efficiency change is never accepted if it weakens a
gate: cheapness bought with quality is not a gain. Gates:

1. hygiene (staged)      - what a commit would record carries no private data, no executed notebooks
2. notebook name order   - no name used before it is defined, no conditional-only binding
3. notebook column audit - every frame column read by name exists in that frame
4. working-copy revision - WARN: the working notebook diverges from HEAD in its *source*
5. synthetic dry-run     - every code cell of the notebook executes against a synthetic dataset
6. pytest                - the repository's synthetic regression tests (skip with --fast)

Gate 4 is a warning, not a failure, because editing a notebook is normal. It exists so that
*overwriting* one cannot happen by accident: a working copy with source changes that are not in HEAD
was authored by hand, and regenerating the file would discard that work.

Usage:
    python tools/agent_gates.py                 # all gates
    python tools/agent_gates.py --fast          # skip the dry-run and pytest
    python tools/agent_gates.py --json
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import pathlib
from datetime import datetime, timezone

REPO = pathlib.Path(__file__).resolve().parent.parent
AGENT_DIR = REPO / ".agent"
NOTEBOOK_DIR = REPO / "analysis" / "notebooks"
PRIMARY_NOTEBOOK = NOTEBOOK_DIR / "06_human_mouse_spatial_rewiring.ipynb"


def run(command: list[str], timeout: int = 300) -> tuple[int, str]:
    """Run a gate. A timeout is a FAIL, not a crash: an unbounded gate is not a gate."""
    try:
        process = subprocess.run(command, cwd=REPO, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return 124, f"timed out after {timeout}s"
    return process.returncode, ((process.stdout or "") + (process.stderr or "")).strip()


class Report:
    """Stream each gate's verdict as it finishes, and keep the machine-readable summary."""

    def __init__(self, as_json: bool):
        self.as_json = as_json
        self.checks: list[dict] = []

    def add(self, name: str, status: str, detail: str) -> None:
        self.checks.append({"name": name, "status": status, "detail": detail})
        if not self.as_json:
            marker = {"pass": "PASS", "fail": "FAIL", "warn": "WARN"}[status]
            print(f"  [{marker}] {name:24s} {detail[:110]}", flush=True)


def _cells(path: pathlib.Path) -> list[str] | None:
    try:
        return ["".join(cell["source"]) for cell in json.load(open(path))["cells"]]
    except Exception:
        return None


def working_copy_divergence() -> tuple[str, str]:
    """Which tracked notebooks differ from HEAD in source (outputs ignored)."""
    affected = []
    for path in sorted(NOTEBOOK_DIR.glob("*.ipynb")):
        relative = str(path.relative_to(REPO))
        code, head_text = run(["git", "show", f"HEAD:{relative}"])
        if code != 0:
            continue
        try:
            head = ["".join(cell["source"]) for cell in json.loads(head_text)["cells"]]
        except json.JSONDecodeError:
            continue
        work = _cells(path)
        if work is None:
            continue
        differing = sum(1 for a, b in zip(work, head) if a != b) + abs(len(work) - len(head))
        if differing:
            affected.append(f"{path.name} ({differing} cells)")
    if not affected:
        return "warn", "working copies match HEAD in source"
    return "warn", ("working copies diverge from HEAD in source: " + ", ".join(affected) +
                    " - do NOT regenerate these files without reconciling first")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fast", action="store_true", help="skip the dry-run and pytest")
    parser.add_argument("--with-slow-tests", action="store_true",
                        help="include the tests CI marks slow (default: the fast lane only)")
    parser.add_argument("--json", action="store_true", help="print the result as JSON")
    args = parser.parse_args()
    python = sys.executable

    notebooks = [str(p.relative_to(REPO)) for p in sorted(NOTEBOOK_DIR.glob("*.ipynb"))]
    report = Report(args.json)

    code, output = run([python, "tools/check_repository_hygiene.py", "--staged"])
    report.add("hygiene (staged)", "pass" if code == 0 else "fail",
               output.splitlines()[-1] if output else "clean")

    code, output = run([python, "tools/check_notebook_stage_order.py", *notebooks])
    report.add("notebook name order", "pass" if code == 0 else "fail",
               output.splitlines()[-1] if output else "clean")

    code, output = run([python, "tools/check_notebook_columns.py",
                        str(PRIMARY_NOTEBOOK.relative_to(REPO))])
    unresolved = [line for line in output.splitlines() if "UNRESOLVED" in line or "collision" in line]
    report.add("notebook column audit", "pass" if code == 0 else "fail",
               " | ".join(unresolved) if unresolved else "clean")

    status, detail = working_copy_divergence()
    report.add("working-copy revision", status, detail)

    if not args.fast:
        code, output = run([python, "tools/notebook_dryrun.py"], timeout=600)
        summary = [line for line in output.splitlines() if line.startswith("dry-run:")]
        report.add("synthetic dry-run", "pass" if code == 0 else "fail",
                   summary[-1] if summary else output.splitlines()[-1][:120])

        # The fast lane the repository's CI uses; --with-slow-tests opts into the rest.
        pytest_command = [python, "-m", "pytest", "tests", "-q"]
        if not args.with_slow_tests:
            pytest_command += ["-m", "not slow"]
        code, output = run(pytest_command, timeout=900)
        summary = [line for line in output.splitlines()
                   if "passed" in line or "failed" in line or "error" in line]
        report.add("pytest", "pass" if code == 0 else "fail",
                   summary[-1][:120] if summary else output.splitlines()[-1][:120])

    checks = report.checks

    failures = [c for c in checks if c["status"] == "fail"]
    warnings = [c for c in checks if c["status"] == "warn"]
    result = {
        "ran_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "passed": not failures,
        "check_count": len(checks),
        "passed_count": len(checks) - len(failures) - len(warnings),
        "checks": checks,
    }
    AGENT_DIR.mkdir(parents=True, exist_ok=True)
    (AGENT_DIR / "gates.json").write_text(json.dumps(result, indent=1) + "\n")

    if args.json:
        print(json.dumps(result, indent=1))
    else:
        print(f"GATES: {'PASS' if result['passed'] else 'FAIL'} "
              f"({result['passed_count']}/{result['check_count']} passed, {len(warnings)} warning)")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
