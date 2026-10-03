"""Hold the guides to the promises the docs make about themselves.

* **Every Python block in a guide runs, or says why it does not.** A ``pycon`` block (a ``>>>``
  session) is executed in CI against the replay bridge (docs/conftest.py). A ``python`` block is
  allowed only with a visible ``title="Not run: <reason>"``, so a reader can tell an example that
  was checked from one that was not, and no block is skipped silently.
* **Every guide opens with a question -> function -> mzLib table**, so a reader can find the call
  for their question before reading the prose.
* **Every guide ends with what to cite**, rendered from the specs' DOIs.
* **No counts or mzLib versions in prose.** "38 formats" was once written in seven places and went
  stale in all of them with one mzLib release. A count belongs in executed output or a generated
  table; a version belongs on the Upgrading page or in the changelog.

Standard library only, so it runs in every test job.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs"
GUIDES = sorted((DOCS / "guides").glob("*.md"))
TEACHING_PAGES = [*GUIDES, DOCS / "getting-started.md", DOCS / "index.md", DOCS / "faq.md"]

FENCE = re.compile(r"^(\s*)(`{3,})(.*)$")
#: Languages that are not Python and need no execution: output, shell, data.
NOT_PYTHON = {"text", "bash", "console", "shell", "powershell", "tsv", "json", "yaml", "toml", "r", "rust"}


def blocks(path: Path):
    """(line number, language, info string) for each fenced block's opening fence."""
    inside = None
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = FENCE.match(line)
        if not match:
            continue
        _, ticks, info = match.groups()
        if inside is None:
            inside = ticks
            language = info.strip().split(" ")[0] if info.strip() else ""
            yield number, language, info.strip()
        elif ticks == inside and not info.strip():
            inside = None


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


@pytest.mark.parametrize("page", TEACHING_PAGES, ids=rel)
def test_every_python_block_runs_or_says_why_not(page):
    problems = []
    for number, language, info in blocks(page):
        if language == "pycon" or language in NOT_PYTHON:
            continue
        if language == "python":
            if not re.search(r'title="Not run: [^"]+"', info):
                problems.append(
                    f"{rel(page)}:{number}: a python block must run (write it as a pycon >>> "
                    'session) or carry title="Not run: <why>"'
                )
            continue
        problems.append(f"{rel(page)}:{number}: fenced block with no language ({info!r})")
    assert not problems, "\n".join(problems)


@pytest.mark.parametrize("guide", GUIDES, ids=rel)
def test_every_guide_opens_with_a_question_table(guide):
    head = "\n".join(guide.read_text(encoding="utf-8").splitlines()[:40])
    assert re.search(r"^\| You want[^|]*\|", head, re.MULTILINE), (
        f"{rel(guide)}: open with a '| You want to ... | Call | mzLib ... |' table in the first 40 lines"
    )


@pytest.mark.parametrize("guide", GUIDES, ids=rel)
def test_every_guide_says_what_to_cite(guide):
    text = guide.read_text(encoding="utf-8")
    assert f'--8<-- "docs/reference/_generated/cite.{guide.stem}.md"' in text, (
        f"{rel(guide)}: end with a '## Cite' section including "
        f"docs/reference/_generated/cite.{guide.stem}.md"
    )


#: Pages whose job is to record versions: the release history and how a release is cut.
VERSIONED_PAGES = {"docs/upgrading.md", "docs/contributing/releasing.md"}

COUNT = re.compile(
    r"\b(?:\d+|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|nineteen)"
    r"\s+(?:file types|file formats|formats|wire verbs|verbs|readers)\b",
    re.IGNORECASE,
)
VERSION = re.compile(r"\bmzLib\s+v?1\.0\.\d+")


def prose_lines(path: Path):
    """Lines outside fenced blocks and HTML comments: what a reader reads as the page's claims."""
    inside = None
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        match = FENCE.match(line)
        if match:
            ticks, info = match.group(2), match.group(3).strip()
            if inside is None:
                inside = ticks
            elif ticks == inside and not info:
                inside = None
            continue
        if inside is None and not line.lstrip().startswith("<!--"):
            yield number, line


PROSE_PAGES = [
    p
    for p in sorted(DOCS.rglob("*.md"))
    if "_generated" not in p.parts and "specs" not in p.parts and rel(p) not in VERSIONED_PAGES
] + [ROOT / "README.md", ROOT / "pkg" / "python" / "README.md"]


@pytest.mark.parametrize("page", PROSE_PAGES, ids=rel)
def test_no_counts_or_mzlib_versions_in_prose(page):
    problems = []
    for number, line in prose_lines(page):
        for pattern, what in ((COUNT, "a count of formats or verbs"), (VERSION, "an mzLib version")):
            for found in pattern.finditer(line):
                problems.append(f"{rel(page)}:{number}: {what} in prose: {found.group(0)!r}")
    assert not problems, (
        "Counts and versions go stale with the next mzLib release. Show the count as executed "
        "output or a generated table; put a version on docs/upgrading.md.\n" + "\n".join(problems)
    )
