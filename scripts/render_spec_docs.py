"""Render the reference fact tables from the vendored verb specs.

    python scripts/render_spec_docs.py           # (re)write docs/reference/_generated/*.md
    python scripts/render_spec_docs.py --check   # exit 1 if that would change anything

One Markdown fragment per spec in docs/specs/, named after it: readers.read-spectra.yaml becomes
docs/reference/_generated/readers.read-spectra.md. Reference pages include a fragment with
pymdownx.snippets (``--8<-- "docs/reference/_generated/readers.read-spectra.md"``); the prose
around it stays hand-written.

The fragments hold only what the spec owns - parameters, units, ranges, result fields, what null
means, error kinds, caveats, the mzLib code it wraps, citations, per-language spellings and
``since`` - in the shared reference-page order every binding uses. Retyping those by hand is how
the bindings came to disagree about how many file types mzLib reads.

The fragments are committed, and CI runs --check, so a spec change without a re-render fails the
build rather than publishing stale facts. The MkDocs hook (scripts/mkdocs_hooks.py) also renders
on every build, so ``mkdocs serve`` is never behind the specs.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import spec_facts  # noqa: E402

MZLIB_BLOB = "https://github.com/smith-chem-wisc/mzLib/blob"
FIXTURE_BLOB = "https://github.com/smith-chem-wisc/pyMzLib/blob/main/pkg/python/tests/fixtures"
DASH = "—"


def cell(value: Any) -> str:
    """A value made safe for one Markdown table cell."""
    if value is None or value == "":
        return DASH
    if isinstance(value, bool):
        return f"`{str(value).lower()}`"
    text = " ".join(str(value).split())
    return text.replace("|", "\\|")


def code(value: Any) -> str:
    if value is None or value == "":
        return DASH
    return "`" + str(value).replace("|", "\\|") + "`"


def table(header: list[str], rows: Iterable[list[str]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return lines


def render_params(spec: dict[str, Any]) -> list[str]:
    params = spec.get("params") or []
    if not params:
        return ["This verb takes no parameters."]
    rows = []
    for p in params:
        name = spec_facts.py_param_name(p["name"])
        shown = code(name) + ("" if name == p["name"] else f" (wire `--{p['name']}`)")
        if p.get("required"):
            default = "required"
        elif p.get("default") is None:
            default = "absent"
        else:
            default = code(
                p["default"] if not isinstance(p["default"], bool) else str(p["default"]).lower()
            )
        rows.append(
            [
                shown,
                code(p.get("type")),
                default,
                cell(p.get("unit")),
                code(p.get("range")) if p.get("range") else DASH,
                cell(p.get("doc")),
            ]
        )
    return table(["Parameter", "Type", "Default", "Unit", "Range", "Meaning"], rows)


def field_rows(fields: list[dict[str, Any]]) -> list[list[str]]:
    rows = []
    for f in fields:
        if f.get("nullable"):
            null = "yes: " + cell(f.get("null_means"))
        else:
            null = "never"
        doc = cell(f.get("doc"))
        if f.get("present_when"):
            doc += f" *Present only with `{f['present_when']}`.*"
        rows.append([code(f["wire"]), code(f.get("type")), cell(f.get("unit")), null, doc])
    return rows


def render_returns(spec: dict[str, Any]) -> list[str]:
    result = spec.get("result") or {}
    header = ["Field", "Type", "Unit", "Null?", "Meaning"]
    out: list[str] = []
    envelope = result.get("envelope_fields") or []
    if envelope:
        out += ["Top-level fields:", ""] + table(header, field_rows(envelope))
    cols = result.get("columns")
    if isinstance(cols, list) and cols:
        out += ["", "Per-row fields (one value per record, in `columns` or each list entry):", ""]
        out += table(header, field_rows(cols))
    elif cols == "per-format":
        out += [
            "",
            "**The columns are per-format**: they are the file's own record fields, listed in "
            "`column_names` on every result. Every cell follows these rules:",
            "",
        ]
        out += [f"- {cell(rule)}" for rule in result.get("cell_rules") or []]
    return out


def render_errors(spec: dict[str, Any]) -> list[str]:
    rows = []
    for e in spec.get("errors") or []:
        kind = e.get("kind")
        exc = spec_facts.PY_EXCEPTION.get(kind)
        when = e.get("when")
        rows.append(
            [
                code(kind),
                f"[`{exc}`][pymzlib.{exc}]" if exc else DASH,
                "never raised by this verb" if when == "never" else cell(when),
            ]
        )
    return table(["Kind", "pyMzLib raises", "When"], rows)


def render(spec: dict[str, Any], commit: str) -> str:
    py = spec_facts.py_binding(spec)
    lines: list[str] = [
        f"<!-- GENERATED by scripts/render_spec_docs.py from docs/specs/{spec['_file']} "
        f"(bridge {commit[:12]}). Do not edit: fix the spec in the bridge, sync it, re-render. -->",
        "",
        f"> {cell(spec.get('summary'))}",
        "",
        f"Wire verb `{spec['verb']}`"
        + (f" · Python [`{py}`][{py}]" if py and spec_facts.shipped_in_python(spec) else "")
        + (
            f" · Python `{py}` (not yet shipped)"
            if py and not spec_facts.shipped_in_python(spec)
            else ""
        ),
        "",
        "### Wraps",
        "",
    ]
    for w in spec.get("wraps") or []:
        pin = str(w.get("pin") or "")
        link = f"{MZLIB_BLOB}/{pin}/{w['path']}" if pin and w.get("path") else None
        target = f"[`{w['symbol']}`]({link})" if link else f"`{w['symbol']}`"
        lines.append(f"- {target} in `{w.get('path', '?')}` at mzLib `{pin[:8] or '?'}`")
    lines += ["", "### Parameters", ""] + render_params(spec)
    lines += ["", "### Returns", ""] + render_returns(spec)
    lines += ["", "### Errors", ""] + render_errors(spec)

    caveats = spec.get("caveats") or []
    lines += ["", "### Caveats", ""]
    lines += [f"- {cell(c)}" for c in caveats] if caveats else ["None recorded."]

    examples = spec.get("examples") or []
    if examples:
        lines += ["", "### Example", ""]
        for ex in examples:
            lines += ["```python", ex.get("py", "").strip(), "```", ""]
            lines.append(
                f"Recorded from the live bridge as [`{ex['fixture']}`]({FIXTURE_BLOB}/{ex['fixture']}). "
                "The doctests in the API reference replay these recordings in CI."
            )
            lines.append("")
        lines.pop()

    if spec.get("performance"):
        lines += ["", "### Performance", "", cell(spec["performance"])]

    b = spec.get("bindings") or {}

    def spelling(lang: str) -> str:
        entry = b.get(lang) or {}
        if not entry.get("name"):
            return DASH
        text = code(entry["name"])
        if entry.get("options"):
            text += f" with `{entry['options']}`"
        return text

    lines += ["", "### Same verb in other bindings", ""]
    lines += table(
        ["Python (pyMzLib)", "Rust (mzLibRust)", "R (mzLibR)"],
        [[spelling("py"), spelling("rust"), spelling("r")]],
    )

    cites = spec.get("cite") or []
    lines += ["", "### References", ""]
    lines += (
        [f"- [doi:{c['doi']}](https://doi.org/{c['doi']}): {cell(c.get('for'))}" for c in cites]
        if cites
        else ["None: no publication is attached to this verb."]
    )

    since = spec.get("since") or {}
    lines += ["", "### Since", ""]
    lines += table(
        ["Wire protocol", "pyMzLib", "mzLibRust", "mzLibR"],
        [
            [
                cell(since.get(k)) if since.get(k) is not None else "not yet shipped"
                for k in ("protocol", "pymzlib", "mzlibrust", "mzlibr")
            ]
        ],
    )

    questions = spec.get("open_questions") or []
    if questions:
        lines += ["", '!!! warning "Not yet verified"', ""]
        lines += [f"    - {cell(q)}" for q in questions]

    return "\n".join(lines) + "\n"


def expected() -> dict[Path, str]:
    commit = spec_facts.source_commit()
    return {
        spec_facts.GENERATED / f"{spec_facts.fragment_name(s)}.md": render(s, commit)
        for s in spec_facts.load_specs()
    }


def write() -> list[str]:
    """Write every fragment and delete fragments whose spec is gone. Returns the changed paths."""
    want = expected()
    spec_facts.GENERATED.mkdir(parents=True, exist_ok=True)
    changed = []
    for path in sorted(spec_facts.GENERATED.glob("*.md")):
        if path not in want:
            path.unlink()
            changed.append(f"removed {path.name}")
    for path, text in want.items():
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            path.write_text(text, encoding="utf-8", newline="\n")
            changed.append(f"wrote {path.name}")
    return changed


def check() -> list[str]:
    want = expected()
    problems = []
    have = set(spec_facts.GENERATED.glob("*.md")) if spec_facts.GENERATED.exists() else set()
    for path in sorted(have - set(want)):
        problems.append(f"{path.name}: no spec in docs/specs/ any more")
    for path, text in want.items():
        if not path.exists():
            problems.append(f"{path.name}: missing")
        elif path.read_text(encoding="utf-8").replace("\r\n", "\n") != text:
            problems.append(f"{path.name}: differs from its spec")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--check", action="store_true", help="fail instead of writing")
    args = ap.parse_args()
    if args.check:
        problems = check()
        for p in problems:
            print(f"  ! {p}")
        if problems:
            print("The reference fact tables are stale. Run: python scripts/render_spec_docs.py")
            return 1
        print(f"{len(expected())} fragments match their specs")
        return 0
    for line in write():
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
