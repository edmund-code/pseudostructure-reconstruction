"""Slide manifest: ``slide_id,image_path,geojson_path,split`` (extra columns preserved).

Required columns are validated up front and every path is checked for existence, so a
typo fails at load with a clear message instead of halfway through patch extraction.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

REQUIRED = ("slide_id", "image_path", "geojson_path", "split")
VALID_SPLITS = ("train", "val", "test")


def load_manifest(path: str | Path, *, check_paths: bool = True,
                  splits: Optional[tuple[str, ...]] = None) -> list[dict]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Manifest not found: {path}\n"
            f"Expected a CSV with columns: {','.join(REQUIRED)}"
        )
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"Manifest {path} has no rows.")

    missing = [c for c in REQUIRED if c not in rows[0]]
    if missing:
        raise ValueError(
            f"Manifest {path} is missing required column(s): {missing}. "
            f"Found: {sorted(rows[0])}"
        )

    out = []
    seen: set[str] = set()
    for i, r in enumerate(rows, start=2):  # line 1 is the header
        sid = (r.get("slide_id") or "").strip()
        if not sid:
            raise ValueError(f"{path}:{i}: empty slide_id")
        if sid in seen:
            raise ValueError(f"{path}:{i}: duplicate slide_id '{sid}'")
        seen.add(sid)
        split = (r.get("split") or "train").strip()
        if split not in VALID_SPLITS:
            raise ValueError(f"{path}:{i}: split '{split}' not in {VALID_SPLITS}")
        row = dict(r)
        row["slide_id"], row["split"] = sid, split
        for key in ("image_path", "geojson_path"):
            p = Path(row[key].strip()).expanduser()
            if check_paths and not p.exists():
                raise FileNotFoundError(f"{path}:{i}: {key} does not exist: {p}")
            row[key] = str(p)
        out.append(row)

    if splits:
        out = [r for r in out if r["split"] in splits]
        if not out:
            raise ValueError(f"No manifest rows with split in {splits}")
    return out


def manifest_summary(rows: list[dict]) -> str:
    by: dict[str, list[str]] = {}
    for r in rows:
        by.setdefault(r["split"], []).append(r["slide_id"])
    return " | ".join(f"{k}: {len(v)} ({', '.join(v)})" for k, v in sorted(by.items()))
