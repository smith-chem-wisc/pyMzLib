"""A stand-in mzLib bridge that answers from recorded fixtures instead of running mzLib.

The docstring examples run in CI against this (pkg/python/src/conftest.py points
``PYMZLIB_BRIDGE`` at a launcher for it), so an example is executed Python, not decoration, and
needs no .NET, no network and no data file. pyMzLib cannot tell it from the real bridge: it is
invoked the same way, with the same arguments, and answers with the same envelope.

    replay_bridge.py readers read-spectra --path sliced_ethcd.mzML --limit 3

``PYMZLIB_REPLAY_TABLE`` names a JSON file mapping each wire verb to its candidate fixtures (the
conftest builds it from the specs' ``examples`` plus its own list for verbs with no spec yet). A
fixture is chosen only if it is **consistent with the call**, so an example cannot print output
recorded for different arguments:

* ``--path`` must name the file the fixture was recorded from (compared by file name), and any
  other ``--option value`` whose snake_case name is a top-level key of the fixture must equal it;
* ``--limit``/``--offset`` must reproduce the fixture's ``returned_count`` from its
  ``record_count`` (or ``row_count``), and a ``--flag`` with a ``<flag>_included`` key must match it;
* a ``--paths-stdin`` call fits only a bulk recording (``files`` + ``read_count``, BULK.md), and a
  one-path call only a one-document recording.

No match, or more than one, is answered as a usage error naming the candidates, so the doctest
fails and says why. Standard library only: it runs under whatever Python runs the tests.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import PurePath, PureWindowsPath


def envelope_error(message: str) -> dict:
    return {
        "ok": False,
        "data": None,
        "error": {"type": "usage", "message": f"replay bridge: {message}"},
    }


def parse(argv: list) -> tuple:
    verb_parts = []
    i = 0
    while i < len(argv) and not argv[i].startswith("--"):
        verb_parts.append(argv[i])
        i += 1
    options = {}
    while i < len(argv):
        name = argv[i][2:]
        if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
            options[name] = argv[i + 1]
            i += 2
        else:
            options[name] = True
            i += 1
    return " ".join(verb_parts), options


def base(path: str) -> str:
    return PureWindowsPath(path).name if "\\" in path else PurePath(path).name


def mismatch(data: dict, options: dict) -> str:
    """Why this recording cannot be the answer to these options, or "" if it can."""
    for name, value in options.items():
        key = name.replace("-", "_")
        if name in ("limit", "offset", "out"):
            continue
        if value is True:
            flag = data.get(f"{key}_included", data.get(key))
            if flag is not None and flag is not True:
                return f"--{name} given, but the recording has {key}={flag!r}"
            continue
        if key not in data or isinstance(data[key], (dict, list)):
            continue
        recorded = data[key]
        if key == "path" or key.endswith("_file"):
            if base(str(recorded)) != base(value):
                return f"recorded from {base(str(recorded))!r}, not {base(value)!r}"
        elif str(recorded) != value:
            return f"recorded with {key}={recorded!r}, not {value!r}"
    # BULK.md: a --paths-stdin call has its own envelope (files[], read_count, ...). Neither shape
    # answers for the other, or a one-path example would "fit" a bulk recording that has no path.
    bulk_recording = isinstance(data.get("files"), list) and "read_count" in data
    if ("paths-stdin" in options) != bulk_recording:
        return "a bulk (--paths-stdin) recording" if bulk_recording else "a one-document recording"
    # A filter the recording applied that the call did not ask for.
    if data.get("ms_order") is not None and "ms-order" not in options:
        return f"recorded with ms_order={data['ms_order']!r}, but the call has no ms-order"
    total = data.get("record_count", data.get("row_count"))
    if total is not None and "returned_count" in data and "out" not in options:
        offset = int(options.get("offset", 0))
        expected = max(0, int(total) - offset)
        if "limit" in options:
            expected = min(expected, int(options["limit"]))
        if int(data.get("offset", 0)) != offset or int(data["returned_count"]) != expected:
            return (
                f"recorded window offset={data.get('offset', 0)} returned={data['returned_count']} "
                f"of {total}, not what limit/offset ask for"
            )
    if data.get("peaks_included") and "peaks" not in options:
        return "recorded with peaks, but the call did not ask for them"
    return ""


def answer(argv: list) -> dict:
    verb, options = parse(argv)
    table_path = os.environ.get("PYMZLIB_REPLAY_TABLE")
    if not table_path:
        return envelope_error("PYMZLIB_REPLAY_TABLE is not set")
    with open(table_path, encoding="utf-8") as fh:
        table = json.load(fh)
    candidates = table.get(verb) or []
    if not candidates:
        return envelope_error(f"no fixture is recorded for '{verb}'")
    fits, reasons = [], []
    for fixture in candidates:
        with open(fixture, encoding="utf-8") as fh:
            data = json.load(fh)
        data = data["data"] if isinstance(data, dict) and "ok" in data and "data" in data else data
        why = mismatch(data, options) if isinstance(data, dict) else ""
        if why:
            reasons.append(f"{base(fixture)}: {why}")
        else:
            fits.append((fixture, data))
    if len(fits) == 1:
        return {"ok": True, "data": fits[0][1], "error": None}
    if not fits:
        return envelope_error(f"no recording of '{verb}' fits {options}: " + "; ".join(reasons))
    return envelope_error(
        f"'{verb}' {options} fits several recordings: " + ", ".join(base(f) for f, _ in fits)
    )


def main() -> int:
    result = answer(sys.argv[1:])
    sys.stdout.write(json.dumps(result))
    return 0 if result["ok"] else 2


if __name__ == "__main__":
    sys.exit(main())
