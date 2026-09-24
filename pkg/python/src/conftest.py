"""Run the docstring examples against recorded bridge output.

    pytest --doctest-modules pkg/python/src/pymzlib

Only doctests are collected from this directory, so everything here applies to them alone; the
test suite in pkg/python/tests never sees it. This file is outside the ``pymzlib`` package and is
not shipped in the wheel.

Each doctest runs with ``PYMZLIB_BRIDGE`` pointing at tests/replay_bridge.py, which answers a call
from a fixture recorded from the real bridge. Which fixtures a verb may answer from comes from the
``examples`` of its spec in docs/specs/, plus ``REPLAY_EXTRA`` below for verbs whose spec has not
been written yet. An example whose arguments fit no recording fails, and the replay bridge says
which recording it tried and why it did not fit.
"""

from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "pkg" / "python" / "tests" / "fixtures"
REPLAY = ROOT / "pkg" / "python" / "tests" / "replay_bridge.py"

#: Recordings for verbs with no spec yet, or recordings a spec does not list as an example.
#: When a verb's spec lands with these fixtures among its ``examples``, delete its entry here.
REPLAY_EXTRA = {
    "readers read-spectra": ["readers_spectra_ms2.json"],
    "readers read-features": ["readers_features_topfd.json"],
    "readers read-matches": [
        "readers_matches_casanovo.json",
        "readers_matches_mspathfinder.json",
        "readers_matches_mzid.json",
    ],
    "pride files": ["pride_PXD000001_files.json"],
    "pride search": ["pride_search_plasmodium.json"],
    "quant flashlfq": ["flashlfq_small.json"],
    "quant median-polish": ["median_polish_small.json"],
    "peptidoform fragments": ["peptidoform_P02768_small.json"],
    "sdrf read": ["sdrf_read_PXD000070.json"],
    "sdrf pool": ["sdrf_pool_two.json"],
}


def replay_table() -> dict:
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import spec_facts
    finally:
        sys.path.pop(0)
    table: dict = {}
    for spec in spec_facts.load_specs():
        for example in spec.get("examples") or []:
            found = spec_facts.find_fixture(example["fixture"])
            if found is None:
                raise pytest.UsageError(
                    f"{spec['_file']} names fixture {example['fixture']}, which is not under pkg/"
                )
            table.setdefault(spec["verb"], []).append(str(found))
    for verb, names in REPLAY_EXTRA.items():
        for name in names:
            path = str(FIXTURES / name)
            if path not in table.setdefault(verb, []):
                table[verb].append(path)
    return table


def launcher(directory: Path) -> Path:
    """An executable that runs the replay bridge, since ``bridge_path()`` wants a file to execute.

    POSIX gets a shell script. Windows gets a real .exe (pip's own console-script launcher): a
    .cmd would route every argument through cmd.exe, which reads the ``|`` in ``trypsin|P`` as a
    pipe. The .exe imports ``replay_bridge`` from tests/, which ``_replay_bridge`` puts on
    PYTHONPATH.
    """
    if sys.platform == "win32":
        from pip._vendor.distlib.scripts import ScriptMaker

        maker = ScriptMaker(None, str(directory))
        maker.executable = sys.executable
        maker.variants = {""}
        return Path(maker.make("replay-bridge = replay_bridge:main")[0])
    path = directory / "replay-bridge"
    path.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{REPLAY}" "$@"\n', encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture(scope="session")
def _replay_bridge(tmp_path_factory):
    directory = tmp_path_factory.mktemp("replay")
    table = directory / "table.json"
    table.write_text(json.dumps(replay_table(), indent=1), encoding="utf-8")
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("PYMZLIB_BRIDGE", str(launcher(directory)))
        mp.setenv("PYMZLIB_REPLAY_TABLE", str(table))
        mp.setenv("PYTHONPATH", str(REPLAY.parent), prepend=os.pathsep)
        yield


@pytest.fixture(autouse=True)
def _doctests_use_the_replay_bridge(_replay_bridge, doctest_namespace, tmp_path, monkeypatch):
    import pymzlib

    doctest_namespace["pymzlib"] = pymzlib
    # Examples that write files (out=, output_directory=) write into a scratch directory.
    monkeypatch.chdir(tmp_path)
    assert os.environ["PYMZLIB_BRIDGE"]
