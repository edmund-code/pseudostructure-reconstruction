#!/usr/bin/env python3
"""Stage output-free copies of tracked notebooks WITHOUT modifying the working files.

Running a notebook fills its cells with outputs, and the hygiene gate rejects committed outputs. The
naive fix -- ``nbconvert --clear-output --inplace`` -- destroys the outputs a contributor may still be
reading, so this stages the cleaned content instead and leaves the working tree untouched:

    python tools/stage_notebooks_without_outputs.py            # clean the notebooks that need it
    python tools/stage_notebooks_without_outputs.py --list     # report, change nothing

Then commit (without ``git add`` on those paths, which would re-stage the dirty working copies) and
check the result with ``python tools/check_repository_hygiene.py --staged``.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def has_outputs(path: Path) -> bool:
    import nbformat
    notebook = nbformat.read(path, as_version=4)
    return any(cell.get("outputs") or cell.get("execution_count") is not None
               for cell in notebook.cells)


def clean_copy(path: Path) -> Path:
    """An output-free copy of `path` in a temporary directory (working file untouched)."""
    import nbformat
    notebook = nbformat.read(path, as_version=4)
    for cell in notebook.cells:
        if cell.get("cell_type") == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
    handle = tempfile.NamedTemporaryFile(suffix=path.suffix, delete=False)
    handle.close()
    copy = Path(handle.name)
    nbformat.write(notebook, copy)
    return copy


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="report without staging")
    parser.add_argument("paths", nargs="*", type=Path,
                        help="notebooks to stage (default: every tracked notebook carrying outputs)")
    args = parser.parse_args(argv)

    if args.paths:
        candidates = [p if p.is_absolute() else REPO_ROOT / p for p in args.paths]
        missing = [p for p in candidates if not p.exists()]
        if missing:
            print("not found:", *missing, sep="\n  ")
            return 1
    else:
        names = subprocess.run(["git", "ls-files", "-z", "*.ipynb"], cwd=REPO_ROOT,
                               check=True, capture_output=True).stdout.split(b"\0")
        candidates = [REPO_ROOT / name.decode() for name in names if name]
    dirty = [path for path in candidates if has_outputs(path)]
    if not dirty:
        print("No tracked notebook carries outputs; nothing to stage.")
        return 0
    for path in dirty:
        print(f"notebook with outputs: {path.relative_to(REPO_ROOT)}")
    if args.list:
        return 0
    for path in dirty:
        copy = clean_copy(path)
        try:
            blob = subprocess.run(["git", "hash-object", "-w", str(copy)], cwd=REPO_ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip()
            # --add so a notebook that is new to the index can be staged too, not only a tracked one.
            subprocess.run(["git", "update-index", "--add", "--cacheinfo", "100644", blob,
                            path.relative_to(REPO_ROOT).as_posix()],
                           cwd=REPO_ROOT, check=True)
        finally:
            copy.unlink(missing_ok=True)
        print(f"staged output-free: {path.relative_to(REPO_ROOT)} (working file left as is)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
