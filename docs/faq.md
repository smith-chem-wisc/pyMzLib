# FAQ

## Do I need .NET installed?

No. The wheel contains a complete .NET runtime. That is most of why it's 115 MB, and it's the
entire point — you should never have to know pyMzLib is C# underneath.

## Do I need a specific Python version?

Python 3.9 or newer, and nothing more specific than that. pyMzLib has no compiled Python code, so
one wheel works on every Python version — unlike most scientific packages, which need a separate
wheel per version. See [D7](design/decisions.md#d7-python-39-and-newer).

## Will it conflict with my other packages?

It cannot. pyMzLib declares **zero** third-party runtime dependencies, so there is nothing for
pip to resolve and nothing to disagree with anything else you have installed. This was a design
requirement, not a happy accident.

## Does it work with conda?

Yes — `pip install pymzlib` inside a conda environment works, and because pyMzLib has no
dependencies, it can't disturb conda's own resolution. A native bioconda package is
[planned](contributing/releasing.md#bioconda-secondary-automatic-after-one-time-setup).

## Why is the wheel so large?

It contains a .NET runtime plus mzLib and its dependencies. The wheel is about 115 MB and
unpacks to roughly 133 MB on disk. For context: mzLib's own NuGet package
is 31 MB, but a C# developer separately downloads TorchSharp and about a gigabyte of libtorch
components — so this single file is *smaller* than what using mzLib from C# costs. On PyPI it's
ordinary: pyOpenMS is 63 MB, torch is 502 MB.

Most of it genuinely is dead weight for the current feature set, for a reason that's about mzLib's
dependency structure rather than about pyMzLib —
[the details](design/decisions.md#d8-payload-size-is-not-a-design-constraint).

## Is it fast?

Each call costs tens of milliseconds of process startup, so pyMzLib is built for coarse-grained
operations — fetch a manifest, download files, run an analysis — not for calling inside a tight
loop. The computation itself runs at C# speed once started.

If you need fine-grained numeric work from Python, please open an issue. That's a real gap and it
needs a different transport underneath, which the
[architecture](design/architecture.md#the-two-load-bearing-properties) deliberately leaves room
for.

## Can I use mzLib feature X?

Only the areas on the [home page](index.md#whats-covered) are exposed so far, deliberately —
coverage grows by demand rather than by guessing. Opening an issue is the fastest way to change
that; the [extension recipe](contributing/adding-a-capability.md) is short.

## Is there an R / Rust / Julia version?

Not today, but the executable pyMzLib drives has a
[language-neutral contract](design/decisions.md#d6-the-wire-contract-stays-language-neutral) and
would work unchanged from another language. Rust is the most likely next binding. If you'd use
one, say so on the issue tracker — that's the evidence that decides it.

## How does this relate to pyOpenMS / pyteomics?

Complementary, not competing. They wrap different libraries: pyOpenMS wraps OpenMS (C++),
pyteomics is pure Python. pyMzLib brings mzLib's particular strengths — the machinery behind
MetaMorpheus, especially top-down and proteoform work — to the same environment. Use whichever
has what you need; nothing stops you using all three, since pyMzLib can't conflict with them.

---

## Troubleshooting

### `BridgeNotFoundError: No mzLib bridge for this platform`

The installed package has no payload for your platform. Almost always one of:

- You installed from an **sdist** instead of a wheel. An sdist can't carry a platform binary.
  Check `pip install --only-binary=:all: pymzlib`.
- You're running from a **source checkout** where the bridge hasn't been built. See
  [building from source](contributing/building.md).
- Your platform isn't among the published wheels (linux-x64, win-x64, osx-x64, osx-arm64).
  Open an issue — adding one is a line of CI configuration.

### `BridgeError: HttpRequestException … status 503`

PRIDE was unavailable. Not a pyMzLib problem; retry later.

### `BridgeError: NotSupportedException … no HTTPS-reachable location`

The file is published only over Aspera. Filter it out before downloading:

```python
downloadable = [f for f in files if f.downloadable]
```

### `PyMzLibError: bridge speaks protocol N, but this pyMzLib expects M`

The Python package and the executable came from different builds. Reinstall the wheel; if you're
working from source, re-run `publish-bridge.ps1` and check `PYMZLIB_BRIDGE` isn't pointing at a
stale executable.

### It hangs

`download()` has no timeout by default, because real transfers can take hours. Pass one if you'd
rather it gave up: `pymzlib.pride.download(..., timeout=3600)`.
