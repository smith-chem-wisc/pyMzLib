"""Run the ``>>>`` examples in the guides against recorded bridge output.

    python -m pytest --doctest-glob="*.md" docs/guides

A guide example written as a ``>>>`` session (a fenced ``pycon`` block) is executed, exactly like a
docstring example: ``PYMZLIB_BRIDGE`` points at the replay bridge, which answers each call from a
fixture recorded from the real bridge, and only if the call's arguments fit that recording. So the
output a guide shows is output pyMzLib actually produced, and a guide that drifts from the code
fails CI rather than misleading a reader.

The fixtures come from ``pkg/python/src/conftest.py`` - one replay table, one rule for which
recording may answer which call - so the guides and the docstrings cannot disagree about it.
Excluded from the site by ``exclude_docs`` in mkdocs.yml.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SOURCE = Path(__file__).resolve().parents[1] / "pkg" / "python" / "src" / "conftest.py"
_spec = importlib.util.spec_from_file_location("_pymzlib_docstring_conftest", _SOURCE)
_docstrings = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_docstrings)

_replay_bridge = _docstrings._replay_bridge
_doctests_use_the_replay_bridge = _docstrings._doctests_use_the_replay_bridge
