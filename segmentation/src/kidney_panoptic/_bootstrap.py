"""Put ``src`` on ``sys.path`` so scripts run without installation.

Every script under ``scripts/`` starts with::

    import _bootstrap  # noqa: F401

which works because the scripts directory is on ``sys.path[0]`` and a sibling
``_bootstrap.py`` shim there imports this module.
"""
from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

REPO_ROOT = _SRC.parent
