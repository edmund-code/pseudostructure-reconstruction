"""YAML config loading with ``_base_`` inheritance and dotted CLI overrides.

A config may declare ``_base_: other.yaml`` (path relative to itself) to inherit
defaults; the child is deep-merged on top. CLI overrides are dotted ``key=value``
strings parsed as YAML scalars, so types survive (``train.lr=1e-4`` is a float,
``decode.merge_fragments=true`` is a bool).
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable

import yaml


class Config(dict):
    """dict with attribute access and a dotted ``get_path`` helper."""

    def __getattr__(self, key: str) -> Any:
        try:
            val = self[key]
        except KeyError as e:
            raise AttributeError(key) from e
        return Config(val) if isinstance(val, dict) and not isinstance(val, Config) else val

    def __setattr__(self, key: str, value: Any) -> None:
        self[key] = value

    def get_path(self, dotted: str, default: Any = None) -> Any:
        node: Any = self
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node


def _deep_merge(base: dict, override: dict) -> dict:
    out = deepcopy(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


def load_config(path: str | Path, overrides: Iterable[str] | None = None) -> Config:
    """Load a YAML config, resolving ``_base_`` inheritance then CLI overrides."""
    path = Path(path)
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    base_name = raw.pop("_base_", None)
    if base_name is not None:
        base_path = (path.parent / base_name).resolve()
        raw = _deep_merge(load_config(base_path), raw)

    # Record where the config came from so relative paths inside it can be resolved.
    cfg = Config(raw)
    cfg["_config_dir"] = str(path.parent.resolve())
    if overrides:
        apply_overrides(cfg, overrides)
    return cfg


def apply_overrides(cfg: dict, overrides: Iterable[str]) -> None:
    """Apply ``a.b.c=value`` strings in place; ``value`` is parsed as a YAML scalar."""
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"Override '{item}' is not of the form key.path=value")
        dotted, raw_val = item.split("=", 1)
        value = yaml.safe_load(raw_val)
        node = cfg
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
            if not isinstance(node, dict):
                raise ValueError(f"Cannot descend into non-dict at '{part}' for '{item}'")
        node[parts[-1]] = value


def resolve_path(cfg: dict, value: str | Path) -> Path:
    """Resolve a possibly-relative config path against the config file's directory."""
    p = Path(value)
    if p.is_absolute():
        return p
    base = cfg.get("_config_dir")
    return (Path(base) / p).resolve() if base else p.resolve()
