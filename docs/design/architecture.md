# Architecture

You don't need any of this to use pyMzLib. It's here because the design has consequences worth
knowing about, and because the reasoning should survive the people who made it.

## The shape of it

```
your code
    │   ordinary Python objects
    ▼
pymzlib.pride                    public API — knows nothing about how the work gets done
    │
pymzlib._bridge                  the ONLY module aware a subprocess exists
    │   JSON over stdout:  {"ok": true, "data": …}  |  {"ok": false, "error": {…}}
    ▼
mzlib-bridge(.exe)               self-contained .NET executable, shipped inside the wheel
    │
mzLib                            the actual proteomics library, in C#
```

When you call `pymzlib.pride.list_files("PXD000001")`, pyMzLib runs the bundled executable with
arguments, that executable calls mzLib's `PrideArchiveClient`, serializes the result as JSON to
standard output, and exits. Python parses it into dataclasses and hands them back.

That's it. It is not clever, and the lack of cleverness is the point.

## Why a subprocess

The obvious alternative is [pythonnet](https://github.com/pythonnet/pythonnet), which loads a
real CLR into the Python process and makes every C# class directly available with no wrapper code
at all. It is genuinely impressive, and it was the wrong choice here.

pythonnet needs a .NET runtime installed on the machine, and since version 3.0 it needs
`Runtime.PythonDLL` or `PYTHONNET_PYDLL` pointed at the exact `python3XX.dll` in use. That is
precisely the class of problem this project exists to eliminate. A wrapper whose install
instructions begin "first, install the .NET runtime, then set an environment variable to the path
of your Python shared library" has already lost, no matter how complete its API surface.

A self-contained single-file executable has none of that. It carries its own runtime, consults no
system installation, and needs no configuration. Costs it does have:

- **Startup.** About 120 ms per call warm, and roughly a second on the very first call of a session while the bundle extracts (measured, not estimated). Irrelevant next to an HTTP request; unacceptable in
  a per-spectrum loop.
- **Serialization.** Everything crossing the boundary must be JSON-representable. Fine for
  manifests and file paths; wrong for a million-point spectrum.
- **Size.** The wheel is about 115 MB because it contains a .NET runtime and mzLib's
  dependencies. See [design decisions](decisions.md#d8-payload-size-is-not-a-design-constraint).

These costs bound where the design applies, and that boundary is honest: pyMzLib is for
coarse-grained operations. When a fine-grained numeric surface is added, the transport underneath
will change — which is why it's confined to one module.

## The two load-bearing properties

### The transport is hidden

`_bridge.py` is the only file that knows a subprocess is involved. Everything above it deals in
Python objects. A future in-process binding, or a long-lived local server for stateful work, can
replace it with no change to any public API. Nothing outside `_bridge.py` may import
`subprocess`, and no public function may leak a process detail — a `CalledProcessError` reaching
a caller would be a bug.

### The wire format is language-neutral

Nothing in the JSON envelope assumes Python. It has no Python-shaped names, no Python types, no
assumptions about the caller. A Rust crate, an R package, or a shell script would use the
identical contract with **no change on the .NET side**:

```bash
$ mzlib-bridge pride files --accession PXD000001
{"ok":true,"data":{"accession":"PXD000001","file_count":8,…}}
```

This costs nothing to maintain and preserves an option worth having. See
[D6](decisions.md#d6-the-wire-contract-stays-language-neutral).

## The envelope

Every invocation writes exactly one JSON object to standard output and nothing else — diagnostics
go to stderr, so a caller can parse stdout without stripping log lines.

```json
{"ok": true,  "data": { … }, "error": null}
{"ok": false, "data": null,  "error": {"type": "HttpRequestException", "message": "…"}}
```

Exit codes: `0` success, `1` handled failure, `2` malformed command.

Failures cross the boundary as **structured data, never a stack trace**. A .NET exception dump is
unparseable noise to a caller in another language; `{"type": "HttpRequestException"}` lets that
caller distinguish a network blip from a bad accession programmatically. The `type` field carries
through to Python as `BridgeError.error_type`.

The envelope is versioned. `mzlib-bridge version` reports a `protocol` number, and pyMzLib checks
it, so a Python package and an executable built from different sources fail loudly rather than
producing subtly wrong results.

## Dependencies, or the lack of them

pyMzLib declares **zero** third-party runtime dependencies. Not few — zero.

Python permits exactly one version of any package per environment. There is no side-by-side
loading and no binding redirect, so two packages with incompatible requirements simply cannot
coexist. Every dependency a library declares is a constraint it imposes on everyone who installs
it. A scientific package that requires `numpy<2` can lock a user out of half the ecosystem.

By requiring nothing, pyMzLib cannot participate in that failure mode at all. This is worth
protecting: a convenience that would cost a runtime dependency is not worth it, and if a future
capability genuinely needs one, it belongs behind an optional extra
(`pip install pymzlib[something]`), never in the base requirements.

## Where the .NET payload lives

Inside the installed package:

```
site-packages/pymzlib/
    __init__.py
    pride.py
    _bridge.py
    _dotnet/
        win-x64/
            mzlib-bridge.exe        ← the whole .NET runtime + mzLib, in one file
            …
```

`bridge_path()` resolves it by platform, and `PYMZLIB_BRIDGE` overrides it for development. Since
the payload is a binary, each wheel is specific to an operating system — but because there is no
compiled *Python* code, one wheel serves **every Python version**. pyMzLib ships one wheel per
OS, not one per OS × Python version, which is why the supported Python range can be generous
without maintenance cost.
