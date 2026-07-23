# Design decisions

Decisions that shaped pyMzLib, with the reasoning attached. They're recorded because the
reasoning is the useful part: a decision without its rationale gets re-litigated or, worse,
quietly reversed by someone who assumes it was arbitrary.

The authoritative copy lives in `.project/state.yaml`.

---

## D1 · Start with one narrow vignette, not a general wrapper

**Decision.** The first release covers PRIDE Archive access only.

**Why.** The interesting risk in a project like this is never "can C# be called from Python" —
it can, several ways. The risk is whether the result can be *delivered* without friction. A tiny
vignette exercises the entire delivery path — build, bundle, wheel, install, run on a foreign
machine — for a fraction of the effort of a broad API. If delivery works, breadth is ordinary
incremental labor. If it doesn't, no amount of API coverage would have saved the project.

---

## D2 · It must install with one command and nothing else

**Decision.** A pip-installable wheel carrying its own self-contained .NET payload, with zero
third-party Python runtime dependencies. Any approach failing this is disqualified regardless of
its other merits.

**Why.** This came directly from the person the project is for: *"I hate having to worry about
which Python I have installed and which libraries are compatible."* Treating that as a
requirement rather than a preference is what killed the otherwise-attractive options — see
[D5](#d5-a-subprocess-not-an-in-process-binding). The consequence is inherited by every future
contribution: a dependency added here is a dependency imposed on every user.

---

## D5 · A subprocess, not an in-process binding

**Decision.** mzLib runs as a self-contained executable invoked per call, speaking a versioned
JSON envelope.

**Why.** [pythonnet](https://github.com/pythonnet/pythonnet) would have exposed all of mzLib for
almost no wrapper code — an enormous temptation. It requires a .NET runtime on the machine and,
since 3.0, an environment variable pointing at the exact `python3XX.dll`. That violates
[D2](#d2-it-must-install-with-one-command-and-nothing-else) outright.

NativeAOT with a C ABI and `cffi` satisfies D2 and is faster, but requires hand-written flat C
signatures and manual marshalling for every function, and NativeAOT is hostile to the reflection
that mzLib's dependencies rely on. It is likely the right answer eventually for fine-grained
numeric work; it is the wrong place to start.

**Cost accepted.** Tens of milliseconds per call, and JSON-representable data only. Both are
irrelevant for repository access and disqualifying for per-spectrum loops — a real boundary on
where the current design applies. See [architecture](architecture.md#why-a-subprocess).

---

## D6 · The wire contract stays language-neutral

**Decision.** Nothing in the JSON envelope may assume a Python caller.

**Why.** The same executable a Python package drives could be driven by a Rust crate, an R
package, or a shell script with no change on the .NET side. Preserving that costs nothing today —
it is a naming and framing discipline, not extra code — and preserves a genuinely valuable
option. A Rust binding is a plausible future direction, and this is what keeps it cheap.

**In practice.** Document each verb as a contract, not as "what `pymzlib` calls". Keep all
transport knowledge inside `_bridge.py`.

---

## D7 · Python 3.9 and newer

**Decision.** `requires-python = ">=3.9"`, deliberately permissive.

**Why.** The prevailing advice is 3.11 or newer, and for most scientific packages that's correct
— but the reason is compiled extensions, which must be rebuilt for every supported Python
version, so each version is real recurring cost. pyMzLib contains no compiled Python code, so one
wheel serves every version and an old floor costs essentially nothing.

What it buys is the audience this project is for: someone on a locked-down cluster with a Python
they can't change. Being the package that refuses to install would be an odd outcome for a
project whose premise is removing friction.

**Cost accepted.** Library code here must avoid post-3.9 syntax — no `match`, no `tomllib`, some
typing features out of reach. That's a discipline on the authors, not the users. CI tests the
floor and the ceiling; with no compiled extension, versions in between can't differ in ways these
tests would catch.

---

## D8 · Payload size is not a design constraint

**Decision.** Don't trade architecture, correctness, or clarity to shrink the wheel.

**Why.** The wheel is about 115 MB, which sounds alarming until it's placed in context. mzLib's
own NuGet package is 31 MB, but NuGet resolves dependencies separately — a C# consumer also pulls
TorchSharp and roughly a gigabyte of libtorch native components. Our single file, which contains
mzLib *and* every dependency *and* the .NET runtime, is smaller than what a C# developer
downloads. On PyPI it is unremarkable: pyOpenMS ships 63 MB, torch ships 502 MB.

**What remains true and worth saying:** most of that payload is dead weight. `torch_cpu.dll`
(238 MB uncompressed) arrives because `UsefulProteomicsDatabases → Proteomics → Omics →
TorchSharp`, in order to make one HTTPS request. A bridge built against only `MzLibUtil` and the
three PRIDE source files is 34 MB with byte-identical behavior. That's a finding about mzLib's
dependency structure — every C# consumer wanting one corner of mzLib pays the same tax — and it
belongs upstream, not worked around here.

**One practical consequence.** PyPI rejects individual files over 100 MB by default. Requesting an
increase is routine and routinely granted; it's a step on the release checklist, not a constraint
on the design.

---

## Deliberately not decided yet

| Question | Why it's open |
|---|---|
| How wide the API surface should get | Needs real demand data. Guessing produces a large surface nobody asked for and can't be removed. |
| Whether to bind Rust as well | Kept cheap by [D6](#d6-the-wire-contract-stays-language-neutral); revisit after the first release. |
| Whether pyMzLib eventually lives inside the mzLib repository | The pyOpenMS precedent (in-tree) is strong, and gets stronger as the surface grows. At three methods, out-of-tree costs less and imposes nothing on mzLib contributors. |
