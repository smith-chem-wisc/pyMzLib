"""Copy the per-verb specs from a bridge checkout into docs/specs/, and record where they came from.

    python scripts/sync_specs.py --from E:/CodeReview/bridge/design/verbs

The specs (one YAML file per wire verb: params, units, result fields, null meanings, errors,
caveats, citations, per-language spellings) are written once, in the bridge repository, so that
pyMzLib, mzLibRust and mzLibR cannot drift apart on facts. That repository is private, and this one
is public, so CI here cannot read it. The copy in docs/specs/ is what the reference renderer
(scripts/render_spec_docs.py), the docstring lint (pkg/python/tests/test_spec_docs.py) and the
doctest replay bridge read.

Never edit docs/specs/*.yaml by hand: a fact that is wrong is wrong in the bridge, and fixing it
there is what makes every binding's lint flag the pages that repeated it. Fix it there, then run
this script.

The copy is a mirror: specs that no longer exist in the source are deleted here. docs/specs/SOURCE
records the bridge commit, and says so when the source had uncommitted changes, so a reviewer can
tell a vendored spec that matches a bridge commit from one that matches someone's working tree.
"""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "docs" / "specs"


def git(source: Path, *args: str) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(source), *args],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return completed.stdout.strip()


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--from",
        dest="source",
        type=Path,
        required=True,
        help="the bridge repository's design/verbs directory",
    )
    args = ap.parse_args()

    source: Path = args.source.resolve()
    specs = sorted(source.glob("*.yaml"))
    if not specs:
        print(f"No *.yaml specs in {source}; is that bridge's design/verbs?", file=sys.stderr)
        return 1

    DEST.mkdir(parents=True, exist_ok=True)
    wanted = {p.name for p in specs}
    for stale in sorted(DEST.glob("*.yaml")):
        if stale.name not in wanted:
            stale.unlink()
            print(f"removed  {stale.name} (no longer in the bridge)")
    for spec in specs:
        target = DEST / spec.name
        before = target.read_bytes() if target.exists() else None
        shutil.copyfile(spec, target)
        if before != target.read_bytes():
            print(f"{'updated' if before is not None else 'added  '}  {spec.name}")

    commit = git(source, "rev-parse", "HEAD") or "unknown (not a git checkout)"
    dirty = git(source, "status", "--porcelain", "--", ".")
    lines = [
        "# Written by scripts/sync_specs.py. Do not edit; re-run the script.",
        "repository: trishorts/bridge (private)",
        "path: design/verbs",
        f"commit: {commit}",
        f"uncommitted_changes: {'yes' if dirty else 'no'}",
        f"synced: {dt.date.today().isoformat()}",
        f"specs: {len(specs)}",
    ]
    (DEST / "SOURCE").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(specs)} specs from {commit[:12]}{' (with uncommitted changes)' if dirty else ''}")
    print("Now run: python scripts/render_spec_docs.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
