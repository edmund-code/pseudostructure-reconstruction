#!/usr/bin/env python3
"""Fail when the clean repository contains private artifacts or notebook output.

By default the tracked files are inspected as they exist in the working tree. With ``--staged`` the
content that a commit would record is inspected instead (the index), so a contributor can keep live
notebook outputs in the working tree -- which is the normal state while running a notebook -- without
either committing them or having them deleted to satisfy this gate.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

MAX_TRACKED_BYTES = 25 * 1024 * 1024
FORBIDDEN_SUFFIXES = {".h5ad", ".geojson", ".tif", ".tiff", ".ndpi", ".pt", ".pth", ".ckpt"}
# Constructed to keep this checker from matching its own policy literals.
FORBIDDEN_PATH_MARKERS = ("/" + "home/", "/" + "Users/")
TEXT_SUFFIXES = {".py", ".md", ".yaml", ".yml", ".csv", ".toml", ".ipynb"}


def tracked_files(root: Path) -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True
    )
    return [root / name.decode() for name in result.stdout.split(b"\0") if name]


def staged_bytes(root: Path, rel: Path) -> bytes | None:
    """The content the index holds for this path, or None when it is not staged."""
    result = subprocess.run(
        ["git", "show", f":{rel.as_posix()}"], cwd=root, capture_output=True
    )
    return result.stdout if result.returncode == 0 else None


def check_content(rel: Path, size: int, raw: bytes) -> list[str]:
    failures: list[str] = []
    if rel.suffix.lower() in FORBIDDEN_SUFFIXES:
        return [f"private/binary artifact tracked: {rel}"]
    if size > MAX_TRACKED_BYTES:
        failures.append(f"tracked file exceeds {MAX_TRACKED_BYTES // (1024 * 1024)} MiB: {rel}")
    if rel.suffix == ".ipynb":
        try:
            import nbformat
            notebook = nbformat.reads(raw.decode("utf-8", errors="replace"), as_version=4)
        except Exception as exc:  # pragma: no cover - reported as a hygiene failure
            failures.append(f"invalid notebook {rel}: {exc}")
        else:
            for index, cell in enumerate(notebook.cells):
                if cell.get("outputs") or cell.get("execution_count") is not None:
                    failures.append(f"executed notebook cell tracked: {rel} cell {index}")
                    break
    if rel.suffix.lower() in TEXT_SUFFIXES:
        text = raw.decode("utf-8", errors="ignore")
        if any(marker in text for marker in FORBIDDEN_PATH_MARKERS):
            failures.append(f"machine-specific path tracked: {rel}")
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged", action="store_true",
                        help="inspect the index (what a commit would record), not the working tree")
    args = parser.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    failures: list[str] = []
    for path in tracked_files(root):
        rel = path.relative_to(root)
        if args.staged:
            raw = staged_bytes(root, rel)
            if raw is None:            # staged as deleted, or never added
                continue
            size = len(raw)
        else:
            if not path.exists():
                continue
            size, raw = path.stat().st_size, path.read_bytes()
        failures.extend(check_content(rel, size, raw))
    if failures:
        print("Repository hygiene failed:", *[f"- {item}" for item in failures], sep="\n")
        return 1
    print("Repository hygiene passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
