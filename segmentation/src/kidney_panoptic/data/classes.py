"""Canonical class vocabulary and the mapping from raw annotation labels onto it.

The GeoJSON files carry historic label names (``Tubules``, ``Tubules2``, ``Tubules3``,
``Glomeruli``, ``Blood Cells``, ...). ``ClassSpec`` maps them onto the canonical
vocabulary declared in the config while REMEMBERING the original name per instance as a
``subtype``, so a later tubule-injury sub-classifier can be trained without re-extracting
patches.

Index 0 is always ``background``; it is a semantic label, never an instance.

Three label groups are handled specially and never become instances:
* ``background_names``  -> rasterized as TRUSTED BACKGROUND (free negatives).
* ``dense_roi_names``   -> mark a fully-annotated region (ignore = 0 inside).
* ``ignore_names``      -> dropped entirely.
"""
from __future__ import annotations

import fnmatch
from dataclasses import dataclass
from typing import Optional

BACKGROUND = "background"


@dataclass(frozen=True)
class ClassSpec:
    names: tuple[str, ...]              # canonical, ordered; names[0] == 'background'
    compact: tuple[str, ...]
    elongated: tuple[str, ...]
    mapping: dict                       # raw label -> canonical name
    background_names: tuple[str, ...]
    dense_roi_names: tuple[str, ...]
    ignore_names: tuple[str, ...]

    # ------------------------------------------------------------------ build
    @classmethod
    def from_config(cls, cfg: dict) -> "ClassSpec":
        names = tuple(cfg.get("classes") or [BACKGROUND])
        if names[0] != BACKGROUND:
            raise ValueError(f"classes[0] must be '{BACKGROUND}', got '{names[0]}'")
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate entries in classes: {names}")

        mapping = dict(cfg.get("class_map") or {})
        # Canonical names always map to themselves so a config can omit identities.
        for n in names:
            mapping.setdefault(n, n)
        unknown = {v for v in mapping.values() if v not in names}
        if unknown:
            raise ValueError(
                f"class_map sends labels to classes not in `classes`: {sorted(unknown)}. "
                f"Known classes: {list(names)}"
            )

        compact = tuple(cfg.get("compact_classes") or ())
        elongated = tuple(cfg.get("elongated_classes") or ())
        for group, label in ((compact, "compact_classes"), (elongated, "elongated_classes")):
            bad = [c for c in group if c not in names]
            if bad:
                raise ValueError(f"{label} contains unknown classes: {bad}")

        return cls(
            names=names,
            compact=compact,
            elongated=elongated,
            mapping=mapping,
            background_names=tuple(cfg.get("background_names") or ()),
            dense_roi_names=tuple(cfg.get("dense_roi_names") or ()),
            ignore_names=tuple(cfg.get("ignore_names") or ()),
        )

    # ------------------------------------------------------------- properties
    @property
    def n_classes(self) -> int:
        return len(self.names)

    @property
    def thing_names(self) -> tuple[str, ...]:
        """Instance-bearing classes (everything but background)."""
        return self.names[1:]

    @property
    def compact_indices(self) -> tuple[int, ...]:
        return tuple(self.names.index(c) for c in self.compact)

    def index(self, name: str) -> int:
        return self.names.index(name)

    # ---------------------------------------------------------------- lookups
    @staticmethod
    def _matches(label: str, patterns: tuple[str, ...]) -> bool:
        return any(label == p or fnmatch.fnmatch(label, p) for p in patterns)

    def role(self, raw_label: str) -> str:
        """Classify a raw annotation label: instance | background | dense_roi | ignore."""
        if self._matches(raw_label, self.ignore_names):
            return "ignore"
        if self._matches(raw_label, self.dense_roi_names):
            return "dense_roi"
        if self._matches(raw_label, self.background_names):
            return "background"
        return "instance"

    def canonical(self, raw_label: str) -> Optional[str]:
        """Canonical class for an instance label, or None if it is not an instance.

        Unmapped instance labels raise — a silently dropped structure class is a data
        bug that would otherwise show up only as unexplained recall loss.
        """
        if self.role(raw_label) != "instance":
            return None
        try:
            return self.mapping[raw_label]
        except KeyError:
            raise KeyError(
                f"Annotation label '{raw_label}' has no entry in class_map and is not a "
                f"canonical class. Known: {sorted(self.mapping)}. Add it to class_map, "
                f"background_names, dense_roi_names, or ignore_names."
            ) from None

    def canonical_index(self, raw_label: str) -> Optional[int]:
        name = self.canonical(raw_label)
        return None if name is None else self.names.index(name)
