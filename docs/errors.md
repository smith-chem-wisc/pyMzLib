# Errors

Every error pyMzLib raises tells you whose problem it is, so you know whether to fix your call,
retry later, or report a bug.

| Exception | Whose problem | What to do |
|---|---|---|
| `UsageError` | **Your call**: a missing file, a malformed accession, a format without the view you asked for, an option without a value | Fix the call. Most are raised before any process starts. |
| `ServiceUnavailableError` | **A remote service**: PRIDE, UniProt or the Gene Ontology server was down, rate-limiting, timed out, or dropped the connection | Retry later. Not your bug, and not pyMzLib's. |
| `BridgeError` | **mzLib**: it read your input and failed, or rejected it. `error_type` is mzLib's own exception type, such as `MzLibException` | Read the message; it usually names the file and line. If it looks like a defect, [report it](https://github.com/smith-chem-wisc/pyMzLib/issues). |
| `BridgeTimeoutError` | **Time**: the call ran past the `timeout` you gave | Give it longer, or a smaller input. pyMzLib will not guess whether the bridge was slow or stuck. |
| `ProjectNotFoundError` | **The accession**: PRIDE has no project by that name, or it has no files | Check the accession for a typo. Raised by `pymzlib.pride` only. |
| `BridgeNotFoundError` | **The install**: no bridge for this platform, or `PYMZLIB_BRIDGE` points at nothing | Reinstall the wheel; in a source checkout, build the bridge. |

They nest, so you can catch as broadly or as narrowly as you like:

```text
PyMzLibError                 every pyMzLib error
├── UsageError               your call
├── BridgeError              mzLib failed; .error_type says how
│   └── ServiceUnavailableError   a remote service failed: retry
├── BridgeTimeoutError
├── BridgeNotFoundError
└── pride.ProjectNotFoundError    no such PRIDE project
```

A typical pattern retries outages and reports everything else:

```pycon
>>> import pymzlib
>>> try:
...     files = pymzlib.pride.list_files("PXD000001")
... except pymzlib.ServiceUnavailableError:
...     files = None          # PRIDE's problem; worth a retry
... except pymzlib.BridgeError as e:
...     raise SystemExit(f"{e.error_type}: {e}")
>>> len(files)
8

```

## Every error, verb by verb

Each verb's spec lists the errors it can raise and when. They are collected here, so you can see in
one place what a call can do wrong. The same tables appear on each verb's reference page. "Kind" is
the bridge's language-neutral name for the error; "pyMzLib raises" is the Python exception.

--8<-- "docs/reference/_generated/index.errors.md"
