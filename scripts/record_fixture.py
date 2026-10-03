"""Record a bridge answer as a test fixture, so docs and doctests replay real output.

    python scripts/record_fixture.py readers_records_crux.json -- readers read-records --path <file> --limit 3
    python scripts/record_fixture.py my.json --stdin lines.txt -- quant median-polish --peptides <file>

Runs the bridge named by ``PYMZLIB_BRIDGE`` (or ``--bridge``) with the arguments after ``--``,
checks that it answered ``ok``, and writes the envelope's ``data`` to
pkg/python/tests/fixtures/<name>, which is the shape every recorded fixture has. The replay bridge
(pkg/python/tests/replay_bridge.py) then answers a call from it only if the call's arguments fit
the recording.

Fixtures are never written by hand: a hand-written fixture documents what someone believed the
bridge says, and the point of replaying them is that the docs show what it does say. Record from a
bridge built at ``code/mzlib.pin``, from data that ships with mzLib's own tests where possible, so a
reader can rerun the call.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "pkg" / "python" / "tests" / "fixtures"


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("name", help="fixture file name, e.g. readers_records_crux.json")
    ap.add_argument("--bridge", default=os.environ.get("PYMZLIB_BRIDGE"), help="bridge executable")
    ap.add_argument("--stdin", type=Path, help="a file whose text is sent on the bridge's stdin")
    ap.add_argument("--force", action="store_true", help="overwrite an existing fixture")
    ap.add_argument(
        "--compact", action="store_true", help="no indentation: for recordings over a few hundred KB"
    )
    argv = sys.argv[1:]
    split = argv.index("--") if "--" in argv else len(argv)
    ns = ap.parse_args(argv[:split])
    args = argv[split + 1 :]
    if not ns.bridge:
        ap.error("set PYMZLIB_BRIDGE or pass --bridge")
    if not args:
        ap.error("give the bridge arguments after --")
    target = FIXTURES / ns.name
    if target.exists() and not ns.force:
        ap.error(f"{target.name} exists; pass --force to re-record it")
    # Always send stdin, empty when not given: a verb that reads its input there must see end of
    # file, not wait on whatever terminal this script was started from.
    stdin = ns.stdin.read_text(encoding="utf-8") if ns.stdin else ""
    completed = subprocess.run(
        [ns.bridge, *args],
        input=stdin,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    try:
        envelope = json.loads(completed.stdout)
    except json.JSONDecodeError:
        print(completed.stdout[:2000], completed.stderr[:2000], sep="\n", file=sys.stderr)
        return 1
    if not envelope.get("ok"):
        print(json.dumps(envelope.get("error"), indent=1), file=sys.stderr)
        return 1
    text = (
        json.dumps(envelope["data"], separators=(",", ":"), ensure_ascii=False)
        if ns.compact
        else json.dumps(envelope["data"], indent=1, ensure_ascii=False)
    )
    target.write_text(text + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {target.relative_to(ROOT)} ({target.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
