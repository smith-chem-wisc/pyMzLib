# Contributing to pyMzLib

We want contributing here to be as easy and transparent as possible — whether that's reporting a
bug, asking for a piece of mzLib to be exposed, or sending a pull request.

Adapted from [mzLib's contributing guide](https://github.com/smith-chem-wisc/mzLib/blob/master/CONTRIBUTING.md),
with the parts that differ for a two-language project spelled out.

## The most useful thing you can do

**Tell us which part of mzLib you want from Python.** Coverage is deliberately partial and grows
by demand rather than by guessing, so a request with a concrete use case genuinely sets
priorities. [Open a capability request](https://github.com/smith-chem-wisc/pyMzLib/issues/new/choose) —
you don't need to know any C#.

## How we work

We use GitHub to host code, track issues and feature requests, and accept pull requests, with
[GitHub Flow](https://guides.github.com/introduction/flow/index.html) — all changes go through a
pull request:

1. Fork the repo and create a branch from `main`.
2. Add tests for anything you add. See [testing](#testing) — coverage is enforced.
3. Update the documentation if you changed behavior or an API.
4. Make sure the test suite passes and the docs build.
5. Open the pull request.

Setup instructions are in
[building from source](https://smith-chem-wisc.github.io/pyMzLib/contributing/building/), and
[adding a capability](https://smith-chem-wisc.github.io/pyMzLib/contributing/adding-a-capability/)
walks the whole path from an mzLib method to a documented, tested Python function.

## Rules that aren't negotiable

Three constraints define this project. A pull request that breaks one won't be merged, however
good it is otherwise — so it's worth knowing them before you write code. Each links to the
reasoning.

1. **No third-party runtime dependency in `pymzlib`.** Not one. Optional extras are fine; base
   requirements are not. Python permits only one version of a package per environment, so every
   dependency we declare is a constraint imposed on everyone who installs us.
   ([D2](https://smith-chem-wisc.github.io/pyMzLib/design/decisions/#d2-it-must-install-with-one-command-and-nothing-else))
2. **The transport stays hidden.** Only `_bridge.py` knows a subprocess exists. Nothing above it
   may import `subprocess`, see an exit code, or leak raw JSON.
   ([architecture](https://smith-chem-wisc.github.io/pyMzLib/design/architecture/))
3. **The wire format stays language-neutral.** Nothing in the JSON envelope may assume a Python
   caller — a Rust or R binding should be able to use the same contract unchanged.
   ([D6](https://smith-chem-wisc.github.io/pyMzLib/design/decisions/#d6-the-wire-contract-stays-language-neutral))

## Testing

Both halves are tested, and both are gated in CI.

| | Command | Gate |
|---|---|---|
| Python | `pytest -m "not network"` | 90% line + branch coverage |
| Python (live) | `pytest -m network` | hits the real PRIDE Archive |
| C# | `dotnet test pkg/bridge.tests/` | 85% over hand-written bridge code |

**Offline tests must stay offline.** They run against recorded fixtures, so they work on a plane
and can't be broken by an EBI outage. Anything touching the network gets `@pytest.mark.network`.

**Both tiers matter, for different reasons.** Offline tests cover the logic — how a call is
assembled, how an error is classified — and they're what runs on every push. Live tests are the
only thing that would notice PRIDE or mzLib changing underneath us. Neither substitutes for the
other.

Some notes on the coverage numbers, because they're easy to misread:

- The C# threshold is lower than the Python one on purpose. The verb handlers construct a
  `PrideArchiveClient` internally and talk to EBI, so they can't be unit-tested without a
  refactor to inject an `HttpClient`. They're covered end-to-end by the Python live tests
  instead. That refactor is planned; the threshold rises when it lands.
- The raw `dotnet test` coverage number is meaningless here — it measures all of mzLib, which has
  its own test suite in its own repository, and reports about 0.6%. Use
  `pkg/build/check-bridge-coverage.ps1`, which scopes to code we actually wrote.
- Chasing 100% isn't the goal. In a wrapper, high line coverage is the easy half: the real risk
  is that the JSON the C# side emits stops matching what Python expects. That's what the wire
  format tests in `pkg/bridge.tests/WireFormatTests.cs` and the live tests exist for.

## Reporting bugs

Report bugs through [GitHub issues](https://github.com/smith-chem-wisc/pyMzLib/issues). Good bug
reports tend to include a summary, steps to reproduce, what you expected, what actually happened,
and anything you already tried. Please include the output of:

```python
import pymzlib; print(pymzlib.__version__, pymzlib.bridge_version())
```

That one line identifies your platform, the bundled .NET runtime, and which mzLib build you have
— which usually narrows the problem immediately.

## Coding style

**Python** — `snake_case`, Google-style docstrings (the API reference is generated from them),
`ruff` for linting, `from __future__ import annotations` in every module so modern type syntax
stays legal on Python 3.9.

**C#** — match mzLib's conventions: XML doc comments on public members that document the failure
contract, explicit argument validation, `ConfigureAwait(false)` on library awaits.

## Licensing of contributions

By contributing, you agree that your contributions will be licensed under this project's
[LGPL-3.0-or-later](LICENSE) license, the same license that covers mzLib. If that's a concern,
please contact the maintainers before opening a pull request.

## Code of Conduct

By participating you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).
