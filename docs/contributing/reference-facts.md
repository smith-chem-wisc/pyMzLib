# How the reference facts are generated, and why

pyMzLib, mzLibRust and mzLibR all project the same bridge, so they describe the same verbs. When
each binding wrote those descriptions by hand, they drifted. By September 2026 Python said mzLib
reads 32 file types and Rust said 31, the R help pages were hand copies of Python docstrings, and
no binding ran any of its examples.

The fix is to write the facts about each wire verb once, and have every binding render them and
check its own prose against them.

## The pieces

```text
bridge (private)                pyMzLib (this repository)
design/verbs/*.yaml  ── sync ──▶ docs/specs/*.yaml + SOURCE
                                    │
                                    ├─ scripts/render_spec_docs.py ─▶ docs/reference/_generated/*.md
                                    │                                   └─ included by docs/reference/<module>.md
                                    ├─ pkg/python/tests/test_spec_docs.py  (docstring lint)
                                    └─ pkg/python/src/conftest.py          (doctest replay bridge)
```

**The spec** (`<module>.<verb>.yaml`) holds the facts about one wire verb:

- its parameters, with their types, defaults, **units** and valid ranges;
- its result fields, with their **units** and **what null means**;
- its error kinds, and when each one happens;
- its caveats, the mzLib code it wraps (at the pin), DOIs, the verb's name in each binding, recorded
  example fixtures, and `since`.

The binding owns the idiom: its names, its prose, its guides. The spec format and its rules are in
the bridge's `design/verbs/README.md`.

**The vendored copy.** The specs live in the bridge repository, which is private, and this
repository's CI cannot read it. So `docs/specs/` carries a copy, and `docs/specs/SOURCE` records the
bridge commit it came from. Never edit that copy. A fact that is wrong is wrong in the bridge, so
fix it there and sync. Every binding then picks up the correction, and every binding's lint flags
the pages that repeated the error.

```bash
python scripts/sync_specs.py --from ../bridge/design/verbs    # wherever your bridge checkout is
python scripts/render_spec_docs.py
```

**The rendered tables.** `render_spec_docs.py` writes one fragment per spec, for example
`docs/reference/_generated/readers.read-spectra.md`. A fragment holds the wraps, parameters,
returns, errors, caveats, example, other-binding spellings, references and since sections, in the
order every binding's reference page uses. A page under `docs/reference/` includes each fragment
with `pymdownx.snippets`, under a heading the page writes itself:

```markdown
## `read_spectra`

--8<-- "docs/reference/_generated/readers.read-spectra.md"
```

The fragments are committed. The MkDocs hook `scripts/mkdocs_hooks.py` re-renders them on every
build, so `mkdocs serve` always shows the current specs. CI fails if a committed fragment differs
from what the specs render.

**The docstring lint** (`pkg/python/tests/test_spec_docs.py`). For every spec that names a
`bindings.py` function, it:

- imports that function;
- parses its Google docstring with griffe (the parser mkdocstrings uses);
- fails if a spec param is missing from the signature or from `Args:`;
- fails if a result field is missing from `Attributes:` on the returned class. When the function
  returns a list of row objects, the per-row fields are checked on the row class instead;
- fails if a param or field whose spec `unit` is not null does not mention that unit. "Skip this
  many" fails; "Scans to skip" passes.

A deliberate difference from the spec, such as a renamed field or a field that the Python shape
absorbs, must be declared with a reason in `PYTHON_DEVIATIONS` in `scripts/spec_facts.py`. The
lint fails on a declaration that names nothing in the spec, so stale entries cannot pile up. The
lint also checks that every verb the spec says pyMzLib ships has its fragment included on a
reference page.

**The executed examples.** `pytest --doctest-modules pkg/python/src/pymzlib` runs every docstring
example. `pkg/python/src/conftest.py` points `PYMZLIB_BRIDGE` at `pkg/python/tests/replay_bridge.py`,
a stand-in bridge that answers each call from a fixture recorded from the real one. The same
fixtures are the ones mzLibRust and mzLibR replay. A verb may answer from the fixtures its spec
lists under `examples`, plus any listed in `REPLAY_EXTRA` in that conftest (for verbs that have no
spec yet). The replay bridge answers only when a recording fits the call: the same file, the same
echoed options, and a `limit`/`offset` that reproduces the recorded window. A mismatch fails the
example and names each recording and why it did not fit. An example therefore cannot print output
recorded for other arguments.

An example that needs the network, a download or pandas keeps `# doctest: +SKIP`. The line above it
says why, as a `>>> #` comment, so a reader can tell a deliberate skip from a forgotten one.

## In CI

`docs.yml` runs, on every pull request:

| check | fails when |
|---|---|
| `render_spec_docs.py --check` | a fragment is stale or missing |
| `pytest tests/test_spec_docs.py` | a docstring lacks a spec param or field, or its unit |
| `pytest --doctest-modules pkg/python/src/pymzlib` | an example raises or prints something else |
| `mkdocs build --strict` | a fragment is not found, or a cross-reference breaks |

`links.yml` runs lychee over every external link weekly, including the mzLib source links and DOIs
that the fragments render. It does not run on pull requests: another site's outage is not a reason
to block a change.

## When you add or change a verb

1. The spec lands in the bridge first, or with your change (SPEC-1). Sync it here.
2. Render. Add the `--8<--` include to `docs/reference/<module>.md`. If the module has no page yet,
   create one and add it to the nav.
3. Write the docstrings so that the lint passes: every param in `Args:`, and every field in the
   result class's `Attributes:`, each with its unit.
4. Give the docstring an example that runs against the spec's recorded fixture. Remove any `+SKIP`
   that no longer applies.

See [Adding a capability](adding-a-capability.md) for the whole recipe.
