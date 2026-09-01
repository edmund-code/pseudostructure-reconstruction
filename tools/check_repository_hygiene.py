#!/usr/bin/env python3
"""Fail when the clean repository contains private artifacts or notebook output."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

MAX_TRACKED_BYTES = 25 * 1024 * 1024
FORBIDDEN_SUFFIXES = {".h5ad", ".geojson", ".tif", ".tiff", ".ndpi", ".pt", ".pth", ".ckpt"}
# Constructed to keep this checker from matching its own policy literals.
FORBIDDEN_PATH_MARKERS = ("/" + "home/", "/" + "Users/")


def tracked_files(root: Path) -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True
    )
    return [root / name.decode() for name in result.stdout.split(b"\0") if name]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    failures: list[str] = []
    for path in tracked_files(root):
        rel = path.relative_to(root)
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            failures.append(f"private/binary artifact tracked: {rel}")
            continue
        if path.stat().st_size > MAX_TRACKED_BYTES:
            failures.append(f"tracked file exceeds {MAX_TRACKED_BYTES // (1024 * 1024)} MiB: {rel}")
        if path.suffix == ".ipynb":
            try:
                import nbformat
                notebook = nbformat.read(path, as_version=4)
            except Exception as exc:  # pragma: no cover - reported as a hygiene failure
                failures.append(f"invalid notebook {rel}: {exc}")
                continue
            for index, cell in enumerate(notebook.cells):
                if cell.get("outputs") or cell.get("execution_count") is not None:
                    failures.append(f"executed notebook cell tracked: {rel} cell {index}")
                    break
        if path.suffix.lower() in {".py", ".md", ".yaml", ".yml", ".csv", ".toml", ".ipynb"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if any(marker in text for marker in FORBIDDEN_PATH_MARKERS):
                failures.append(f"machine-specific path tracked: {rel}")
    if failures:
        print("Repository hygiene failed:", *[f"- {item}" for item in failures], sep="\n")
        return 1
    print("Repository hygiene passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
