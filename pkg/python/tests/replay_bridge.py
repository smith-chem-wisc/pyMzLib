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
  one-path call only a one-document recording;
* a recording with a ``written`` key fits a ``--out`` call only if ``written`` is set, and a call
  without ``--out`` only if it is null;
* a recording listed in ``RECORDED_WITH`` fits only a call that sends the options it was made with;
* an optional input file in ``OPTIONAL_FILES`` (echoed as ``<name>_file``, null when absent) must
  be given by the call exactly when the recording used one;
* a recording with a key in ``GIVEN_ECHOES`` (``thresholds_given``) fits only a call that sends
  its option exactly when that key is true, and a numeric echo matches a call's option by value
  (``--ratio 1.0`` fits a recorded ``ratio`` of 1).

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


#: Wire options a verb echoes under another name, so a recording made with one value cannot answer
#: a call with another (peptidoform fragments echoes --max-mods as max_modifications).
ECHOED_AS = {
    "max-mods": "max_modifications",
    "max-isoforms": "max_modification_isoforms",
    "psms": "psm_file",
    "peptides": "peptides_file",
    "from": "source_format",
    "to": "target_format",
}


#: Recordings made with an option the call must send as well. The echo alone cannot say so: a call
#: that leaves the option out takes a default (mzLib's) that the replay does not know, and the echo
#: of a value that was given looks the same as the echo of a default. So a call without
#: --max-mods cannot be answered from the max-mods 1 albumin digest.
RECORDED_WITH = {
    "peptidoform_P02768_small.json": {"max-mods": "1"},
    "proteins_align_charged.json": {"end-gaps": "charged"},
    # library read echoes both, but a call that leaves them out takes the defaults, and the echo
    # would not stop a recording made with them from answering it.
    "library_read_fragments.json": {"fragments": True},
    "library_read_index_only.json": {"index-only-threshold": "0"},
    "stats_fit_mixed_lmm_ml.json": {"estimator": "ml"},
}


#: Optional input files a verb echoes as ``<name>_file``, null when the option was not given (stats
#: correlate --given). A recording made with one answers only a call that gives one, and the
#: reverse, or a plain correlation could print a partial one.
OPTIONAL_FILES = ("given",)


#: Echoes that say whether an option was given, when its value echo cannot (stats entrapment-fdp
#: echoes the thresholds it used, which are mzLib's default when --thresholds was absent).
GIVEN_ECHOES = {"thresholds_given": "thresholds"}


def same_number(recorded: object, value: str) -> bool:
    """A numeric echo against an option's text: JSON writes 1.0 as 1, Python sends repr(1.0)."""
    if isinstance(recorded, bool) or not isinstance(recorded, (int, float)):
        return False
    try:
        return float(value) == float(recorded)
    except ValueError:
        return False


def mismatch(data: dict, options: dict) -> str:
    """Why this recording cannot be the answer to these options, or "" if it can."""
    for key, name in GIVEN_ECHOES.items():
        given = data.get(key)
        if isinstance(given, bool) and (name in options) != given:
            return f"recorded {'with' if given else 'without'} --{name}"
    for name in OPTIONAL_FILES:
        key = f"{name}_file"
        if key in data and (name in options) != (data[key] is not None):
            return f"recorded {'with' if data[key] is not None else 'without'} --{name}"
    for name, value in options.items():
        key = ECHOED_AS.get(name, name.replace("-", "_"))
        if name in ("limit", "offset", "out"):
            continue
        if value is True:
            flag = data.get(f"{key}_included", data.get(key))
            if flag is not None and flag is not True:
                return f"--{name} given, but the recording has {key}={flag!r}"
            continue
        # An input file option is echoed back as <name>_file (--peptides -> peptides_file,
        # --responses -> responses_file): hold the call's file to the recording's, by file name.
        if key not in data and isinstance(data.get(f"{key}_file"), str):
            if base(data[f"{key}_file"]) != base(value):
                return f"recorded from {base(data[key + '_file'])!r}, not {base(value)!r}"
            continue
        if key not in data or isinstance(data[key], (dict, list)):
            continue
        recorded = data[key]
        if key == "path" or key.endswith("_file"):
            if base(str(recorded)) != base(value):
                return f"recorded from {base(str(recorded))!r}, not {base(value)!r}"
        elif str(recorded) != value and not same_number(recorded, value):
            return f"recorded with {key}={recorded!r}, not {value!r}"
    # BULK.md: a --paths-stdin call has its own envelope (files[], read_count, ...). Neither shape
    # answers for the other, or a one-path example would "fit" a bulk recording that has no path.
    # A recording is bulk when it carries BULK.md's per-input files[] (entries with a path) and no
    # top-level path. read_count is not required: a verb that refuses on_error="skip" (proteins
    # classify-peptides) has nothing to count. pride files also has a files list, of PRIDE files.
    files = data.get("files")
    # An empty list says nothing either way (PRIDE's answer for an unknown accession is files=[]),
    # so it counts as bulk only when BULK.md's read_count says so.
    bulk_recording = (
        isinstance(files, list)
        and (bool(files) or "read_count" in data)
        and all(isinstance(f, dict) and "path" in f for f in files)
        and "path" not in data
    )
    if ("paths-stdin" in options) != bulk_recording:
        return "a bulk (--paths-stdin) recording" if bulk_recording else "a one-document recording"
    # A recording that wrote a file answers only a call that asked for one, and the reverse: the
    # payload's `written` block is the evidence, so an out= example cannot print a recording that
    # wrote nothing.
    if "written" in data and ("out" in options) != (data["written"] is not None):
        return "a recording that wrote out=" if data["written"] is not None else "a recording without out="
    # proteins read: a filtered read (--accessions-stdin) and an unfiltered one never stand in for
    # each other, or an example would print rows its filter did not select.
    if "accession_filter_count" in data:
        if ("accessions-stdin" in options) != (data["accession_filter_count"] is not None):
            return "a recording with the other accession filter (filtered vs unfiltered)"
        # --sequences adds a sequence column: a recording with it and one without never stand in
        # for each other, or an example would print sequences it did not ask for (or lack them).
        if ("sequences" in options) != ("sequence" in (data.get("column_names") or [])):
            return "a recording with the other --sequences choice"
    # proteins align: --no-residue-map leaves residue_map null, and only then.
    if "residue_map" in data and ("no-residue-map" in options) != (data["residue_map"] is None):
        return "a recording with the other --no-residue-map choice"
    # A filter the recording applied that the call did not ask for.
    if data.get("ms_order") is not None and "ms-order" not in options:
        return f"recorded with ms_order={data['ms_order']!r}, but the call has no ms-order"
    if data.get("kit") is not None and "kit" not in options:
        return f"recorded for kit={data['kit']!r}, but the call asks for every kit"
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


#: Verbs whose input travels on stdin and is echoed in the recording, so a recording answers only a
#: call that sent the same input. Read only for these: other verbs may inherit a terminal's stdin.
STDIN_ECHO = {"quant flashlfq", "peptidoform convert", "proteins parse-accession", "proteins align",
              "stats combine-p", "stats fit-logistic", "stats fit-mixed"}


def stdin_mismatch(data: dict, stdin: str, verb: str = "quant flashlfq") -> str:
    """The call's stdin must be the recording's input.

    quant flashlfq: the same mzML runs (by file name). peptidoform convert and proteins
    parse-accession: the same lines, in the same order, as the recording's ``input`` column.
    proteins align: the same pairs, in the same order, as the recording's aligned rows with their
    gaps removed (mzLib refuses a '-' in an input sequence, so removing them is exact).
    stats combine-p: the same features in the same first-appearance order, the same number of
    tests, and weights and effects given or not as the recording's ``weighted`` and ``signed`` say.
    stats fit-logistic and fit-mixed: the same coefficients, in the same order, as the recording's
    ``tested_coefficients``.
    """
    if verb in ("stats fit-logistic", "stats fit-mixed"):
        sent = [line.strip() for line in stdin.splitlines() if line.strip()]
        recorded = list(data.get("tested_coefficients") or [])
        if sent != recorded:
            return f"recorded for coefficients {recorded}, not {sent}"
        return ""
    if verb == "stats combine-p":
        lines = [line.split("\t") for line in stdin.splitlines() if line.strip()]
        features = list(dict.fromkeys(fields[0].strip() for fields in lines))
        recorded = list((data.get("columns") or {}).get("feature") or [])
        if features != recorded:
            return f"recorded for features {recorded}, not {features}"
        if len(lines) != data.get("line_count"):
            return f"recorded for {data.get('line_count')} tests, not {len(lines)}"
        weighted = any(len(f) >= 3 and f[2].strip().lower() not in ("", "na", "nan") for f in lines)
        signed = any(len(f) == 4 for f in lines)
        if (weighted, signed) != (data.get("weighted"), data.get("signed")):
            return (
                f"recorded with weighted={data.get('weighted')} signed={data.get('signed')}, "
                f"not weighted={weighted} signed={signed}"
            )
        return ""
    if verb in ("peptidoform convert", "proteins parse-accession"):
        sent = [line for line in stdin.splitlines() if line.strip()]
        recorded = list((data.get("columns") or {}).get("input") or [])
        if sent != recorded:
            return f"recorded for inputs {recorded}, not {sent}"
        return ""
    if verb == "proteins align":
        sent = [tuple(line.split("\t")) for line in stdin.splitlines() if line.strip()]
        columns = data.get("columns") or {}
        recorded = [
            (s.replace("-", ""), t.replace("-", ""))
            for s, t in zip(columns.get("aligned_source") or [], columns.get("aligned_target") or [])
        ]
        if sent != recorded:
            short = [(s[:12], t[:12]) for s, t in recorded]
            return f"recorded for pairs starting {short}, not these"
        return ""
    sent = {base(line.split("	")[0].strip()) for line in stdin.splitlines() if line.strip()}
    recorded = {base(str(f.get("full_path", ""))) for f in data.get("spectra_files") or []}
    if sent != recorded:
        return f"recorded for runs {sorted(recorded)}, not {sorted(sent)}"
    return ""


def answer(argv: list, stdin: str = "") -> dict:
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
        for name, value in RECORDED_WITH.get(base(fixture), {}).items():
            if not why and options.get(name) != value:
                why = f"recorded with --{name} {value}, which the call does not send"
        if not why and verb in STDIN_ECHO and isinstance(data, dict):
            why = stdin_mismatch(data, stdin, verb)
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
    argv = sys.argv[1:]
    verb, _ = parse(argv)
    stdin = sys.stdin.read() if verb in STDIN_ECHO else ""
    result = answer(argv, stdin)
    sys.stdout.write(json.dumps(result))
    return 0 if result["ok"] else 2


if __name__ == "__main__":
    sys.exit(main())
