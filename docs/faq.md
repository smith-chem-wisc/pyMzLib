# FAQ

## Do I need .NET installed?

No. The wheel contains a complete .NET runtime — about half of its ~60 MB, and the entire
point: you should never have to know pyMzLib is C# underneath.

## Do I need a specific Python version?

Python 3.9 or newer, and nothing more specific than that. pyMzLib has no compiled Python code, so
one wheel works on every Python version — unlike most scientific packages, which need a separate
wheel per version. See [D7](design/decisions.md#d7-python-39-and-newer).

## Will it conflict with my other packages?

It cannot. pyMzLib declares **zero** third-party runtime dependencies, so there is nothing for
pip to resolve and nothing to disagree with anything else you have installed. This was a design
requirement, not a happy accident.

## Does it work with conda?

Yes — `pip install mzlib` inside a conda environment works, and because pyMzLib has no
dependencies, it can't disturb conda's own resolution. A native bioconda package is
[planned](contributing/releasing.md#bioconda-secondary-automatic-after-one-time-setup).

## Why is the wheel so large?

It contains a .NET runtime plus mzLib and its dependencies. The wheel is about 60 MB. For
context: mzLib's own NuGet package is 31 MB, but a C# developer separately downloads TorchSharp
and about a gigabyte of libtorch components — so this single file is far *smaller* than what using
mzLib from C# costs. On PyPI it's ordinary: pyOpenMS is 63 MB, torch is 502 MB.

It used to be up to 166 MB. Nearly all of that was libtorch, a machine-learning library mzLib
pulls in for one retention-time predictor that pyMzLib doesn't expose, so it's no longer shipped —
[the details](design/decisions.md#d8-payload-size-is-not-a-design-constraint).

## Is it fast?

Each call costs about 120 ms of process startup (roughly a second on the first call of a session), so pyMzLib is built for coarse-grained
operations — fetch a manifest, download files, run an analysis — not for calling inside a tight
loop. The computation itself runs at C# speed once started.

If you need fine-grained numeric work from Python, please open an issue. That's a real gap and it
needs a different transport underneath, which the
[architecture](design/architecture.md#the-two-load-bearing-properties) deliberately leaves room
for.

## How do I read a hundred files? Should I use a thread pool?

Use the function's `_many` twin - `read_spectra_many(paths)`, `read_records_many(paths)`, and so
on - and let it do the parallel work:

```python
batch = pymzlib.readers.read_spectra_many(paths, threads=4)
```

Do not wrap the single-file function in a `ThreadPoolExecutor` or a `multiprocessing.Pool`. It
works, but each call starts its own bridge process, so a hundred files pay the ~120 ms start-up a
hundred times, and nobody is counting threads: eight Python workers each driving a bridge that
already reads on every core is eight times more threads than cores. A `_many` call is one process
that reads `threads` files at a time and hands back one table. See
[Many files at once](guides/readers.md#many-files-at-once).

## Why is `threads` 1 by default?

Because it costs memory, not correctness. Every mzLib reader holds a whole file in memory while it
works, so `threads=8` can mean eight whole files at once - fine for search results, a lot for
multi-gigabyte `.raw` files. The answer never depends on it: the table is byte-identical at any
thread count, because files are always returned in the order you listed them. So the default is
the setting that cannot run you out of memory, and raising it is only a speed decision. For many
small files use your core count (or `-1`); for large spectra files start at 2, since mzLib's mzML
and Thermo readers already parallelise inside each file.

That is different from FlashLFQ's `max_threads`, where the thread count can change the numbers
themselves; see the [FlashLFQ guide](guides/flashlfq.md).

## Why is a whole column `None`?

Look at the result's `absent_fields`. A column named there is one the function defines but **this
file's format has no column for** - MBR Score in a current FlashLFQ peaks table, `q_value` in an
MSPathFinder targets file, apex intensity in a FLASHDeconv feature file - so every value is `None`
rather than the default mzLib would have filled in (often a zero that looks like a measurement).
`failed_fields` names columns whose read threw on some rows, and `excluded_fields` names fields
that have no column shape at all. The
[readers guide](guides/readers.md#four-ways-a-field-can-have-no-value) sets the four cases side
by side.

## Can I use mzLib feature X?

Only the areas on the [home page](index.md#whats-covered) are exposed so far, deliberately —
coverage grows by demand rather than by guessing. Opening an issue is the fastest way to change
that; the [extension recipe](contributing/adding-a-capability.md) is short.

## Why are my GO terms empty, or every gene `not_in_source`?

You read a FASTA. A FASTA header carries organism (`OS=`), taxonomy id (`OX=`) and gene name
(`GN=`) and nothing else, so it has no GO terms and no Ensembl links to report. Each result says so
in `files[i].absent_fields` rather than letting an empty table pass for "no annotation". Read the
UniProt XML of the same proteome instead: the accessions are the same. See
[Protein databases](guides/proteins.md#before-you-start-xml-or-fasta).

## Is peptide uniqueness decided with I and L as different residues?

No. `pymzlib.proteins.classify_peptides` treats I and L as one residue, because they have the same
mass. A peptide written with L finds a protein that has I at that position. See
[the rules](guides/proteins.md#the-rules-exactly).

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

## My SDRF passes `validate()`, but `assess()` calls it a Skeleton. Which is right?

Both. [`validate()`](guides/sdrf.md#validate-a-deposit) checks **structure**, and a file whose
sample columns all say `"not available"` is structurally perfect: reserved words are the
specification's correct way to say there is no value. [`assess()`](guides/sdrf.md#decide-whether-an-sdrf-is-worth-using)
asks whether the file **says anything** about its samples, and that file does not. Use `validate()`
before you deposit and `assess()` before you rely on a file to group results by biology.

## Why is the age `None` when the cell says `63`?

Because `63` has no unit, and 63 years and 63 days are both plausible in one study. mzLib refuses
to guess, and `refusal` says `"no_unit"` so you can tell this from a curator writing
`"not available"`. The same goes for free text. Ask the data's author, or fix the cell to `63Y`;
do not assume years. See [Parse ages safely](guides/sdrf.md#parse-ages-safely).

## How do I check hundreds of SDRF files without waiting on hundreds of processes?

Use the `*_many` forms — `validate_many()`, `assess_many()`, `samples_many()` — which read the
whole list in one bridge call, with `threads=-1` for every core and `on_error="skip"` to carry on
past a bad file. Do not wrap the single-file functions in a thread pool: each call re-pays about
120 ms of start-up, and the thread count belongs inside the bridge, where it cannot change the
answer. See [Many files in one call](guides/sdrf.md#many-files-in-one-call).

---

## Troubleshooting

### `BridgeNotFoundError: No mzLib bridge for this platform`

The installed package has no payload for your platform. Almost always one of:

- You installed from an **sdist** instead of a wheel. An sdist can't carry a platform binary.
  Check `pip install --only-binary=:all: mzlib`.
- You're running from a **source checkout** where the bridge hasn't been built. See
  [building from source](contributing/building.md).
- Your platform isn't among the published wheels (linux-x64, win-x64, osx-x64, osx-arm64).
  Open an issue — adding one is a line of CI configuration.

### `ServiceUnavailableError: … status 503`

PRIDE was unavailable. Not a pyMzLib problem; retry later.

This arrives as `ServiceUnavailableError`, **not** a plain `BridgeError` — the bridge classifies
408, 429 and 5xx as availability failures so you can retry them and report everything else. Since
`ServiceUnavailableError` subclasses `BridgeError`, `except BridgeError` still catches it; catch the
narrower one when you want to retry. See [the PRIDE guide](guides/pride.md#errors-you-might-hit).

### `ServiceUnavailableError: Received an unexpected EOF or 0 bytes from the transport stream`

The connection dropped **part-way through a download**. The server accepted the request and started
sending, then went away — so it is an outage like any other, and retrying is the right response.

Large downloads are the usual place to meet it, simply because they are exposed for longer. If it
repeats at roughly the same point every time, suspect something between you and EBI (a proxy or
scanner cutting long transfers) rather than EBI itself.

!!! note "Before pyMzLib 0.1.0.dev4 this surfaced as `BridgeError`"
    Older versions reported it under the raw .NET type name, which made an outage look like a
    contract break — retry loops written around `ServiceUnavailableError` did not catch it. If you
    wrote a workaround that also catches `BridgeError` for this message, you can drop it.

### `BridgeError: NotSupportedException … no HTTPS-reachable location`

The file is published only over Aspera. Filter it out before downloading:

```python
downloadable = [f for f in files if f.downloadable]
```

### `PyMzLibError: bridge speaks protocol N, but this pyMzLib expects M`

The Python package and the executable came from different builds. Reinstall the wheel; if you're
working from source, re-run `publish-bridge.ps1` and check `PYMZLIB_BRIDGE` isn't pointing at a
stale executable.

### `UsageError: 'readers read-…' needs the bridge from pyMzLib 0.2.0 or later`

`PYMZLIB_BRIDGE` points at a bridge built from an older pyMzLib, which does not have the function
you called. pyMzLib asks the bridge which commands it has before calling a new one, so this is
raised before anything runs. Rebuild the bridge from this source tree (`publish-bridge.ps1`), or
unset `PYMZLIB_BRIDGE` to use the one in the wheel. `pymzlib.bridge_version()["verbs"]` lists what
the bridge in use can do.

### It hangs

`download()` has no timeout by default, because real transfers can take hours. Pass one if you'd
rather it gave up: `pymzlib.pride.download(..., timeout=3600)`.
