"""MkDocs hooks: render the per-verb fact tables before every build.

Wired up by ``hooks:`` in mkdocs.yml. The fragments are committed too, and CI separately runs
``render_spec_docs.py --check``, so a stale commit fails even though the build itself would have
rendered fresh tables. See docs/contributing/reference-facts.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import render_spec_docs  # noqa: E402


def on_pre_build(config, **kwargs):  # noqa: ARG001 - MkDocs hook signature
    for line in render_spec_docs.write():
        print(f"render_spec_docs: {line}")
