"""Every README's coverage table must keep up with the package.

This exists because they did not. `readers` and FlashLFQ quantification both shipped, both got a
documentation guide, and both were listed on the docs site — while the repository README still
advertised two areas out of four, and the PyPI README still advertised **one**. Nobody noticed from
inside the repository, because every other surface was correct; it took someone opening the GitHub
landing page.

That asymmetry is the thing worth defending against. `docs/` is docs-as-code: built, linked, and
`--strict`-checked in CI, so a stale page tends to announce itself. The READMEs are checked by
nothing, and they are the first — often only — thing a prospective user reads. A capability that is
implemented, tested, documented and *invisible* is, from outside, not shipped.

Both READMEs are checked, not just the repository one. They are different files with different
audiences (GitHub and PyPI), and the PyPI one had drifted further precisely because it is the one
nobody working in the repo ever looks at.

Deliberately coarse: it asserts each public capability is *mentioned*, not how it is described.
Pinning the prose would fail on every wording change and teach people to edit the test until it
passes, which is worse than no test.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import pymzlib

#: The capability modules a user can reach. Anything here is a thing pyMzLib can do, so it belongs
#: on the front page. Everything else in ``__all__`` is plumbing — errors, the bridge locator, the
#: version — which no coverage table should list.
CAPABILITY_MODULES = ("flashlfq", "peptidoform", "pride", "readers")

#: Repo-root-relative. Both are front pages; they simply face different registries.
README_PATHS = ("README.md", "pkg/python/README.md")


def _repository_root() -> Path | None:
    """The checkout root, identified by a marker no installed copy carries.

    Walking up to the *nearest* README is what the first version of this test did, and it silently
    read ``pkg/python/README.md`` for both cases — passing or failing for the wrong file. Anchoring
    on a directory that only exists in a source checkout makes the wrong answer impossible rather
    than unlikely.
    """
    for directory in Path(__file__).resolve().parents:
        if (directory / ".github").is_dir() and (directory / "README.md").is_file():
            return directory
    return None


def _coverage_table(readme: str, name: str) -> str:
    """Just the 'What's covered' section, so a stray mention elsewhere cannot satisfy the test."""
    match = re.search(
        r"^##\s+What's covered\s*$(.*?)(?=^##\s)", readme, re.MULTILINE | re.DOTALL
    )
    assert match, f"{name} no longer has a \"What's covered\" section; this test needs updating"
    return match.group(1)


def test_the_capability_modules_are_all_exported():
    """Guards the list above against a module being renamed out from under this test."""
    for module in CAPABILITY_MODULES:
        assert module in pymzlib.__all__, f"pymzlib.__all__ no longer exports {module!r}"


@pytest.mark.parametrize("relative_path", README_PATHS)
def test_every_capability_appears_in_the_readme_coverage_table(relative_path):
    root = _repository_root()
    if root is None:
        pytest.skip("not a source checkout; there is no README to check")

    readme = (root / relative_path).read_text(encoding="utf-8")
    table = _coverage_table(readme, relative_path).lower()

    missing = [module for module in CAPABILITY_MODULES if module not in table]

    assert not missing, (
        f'{relative_path}: the "What\'s covered" table does not mention {missing}. '
        "A capability that ships without reaching the front page is invisible to everyone who has "
        "not read the source. Add a row, or drop the module from CAPABILITY_MODULES if it is no "
        "longer a public capability."
    )
