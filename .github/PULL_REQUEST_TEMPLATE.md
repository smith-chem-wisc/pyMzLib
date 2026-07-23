## What this changes

<!-- One or two sentences. If it fixes an issue, "Fixes #123". -->

## Why

<!-- The reasoning, not the diff. Reviewers can read the diff; they can't read your mind. -->

---

### Checklist

- [ ] Tests added or updated, and both suites pass:
      `pytest -m "not network"` and
      `dotnet test pkg/bridge.tests/MzLibBridge.Tests.csproj --filter "TestCategory!=ExternalService"`
- [ ] `./pkg/build/check-bridge-coverage.ps1` passes (C# coverage is gated separately from Python)
- [ ] Coverage gates still met (Python 90%, bridge 85% — see CONTRIBUTING.md)
- [ ] Documentation updated if behavior or an API changed
- [ ] `mkdocs build --strict` passes if docs were touched

### The three rules

Confirm none of these are broken — a PR that breaks one can't be merged, however good it is
otherwise. Reasoning is in [design decisions](https://smith-chem-wisc.github.io/pyMzLib/design/decisions/).

- [ ] **No new third-party runtime dependency** in `pymzlib` (optional extras are fine)
- [ ] **The transport stays hidden** — nothing outside `_bridge.py` knows a subprocess exists
- [ ] **The wire format stays language-neutral** — nothing in the JSON envelope assumes Python

### If this adds a capability

- [ ] It's coarse-grained enough that per-call process startup is irrelevant
- [ ] No mzLib logic reimplemented on the bridge side
- [ ] Offline test with a recorded fixture, plus a live test for the real path
- [ ] Guide page, nav entry, and coverage table updated
