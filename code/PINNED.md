# Code worktrees + pins — pyMzLib

Worktrees under `code/` are gitignored; this file is their record.

## The pin

**The mzLib commit is in [`code/mzlib.pin`](mzlib.pin), and nowhere else.**

That file is the single authoritative value: `.github/workflows/wheels.yml` reads it in every job
that builds the bridge, and `upstream-watch.yml` is the only thing that writes it. This document
carries the *provenance* — which release, which PRs, and why — but deliberately no longer restates
the sha itself.

It used to. The sha was written in three places (this table, `wheels.yml`, and the project
workspace's own copy of this file), and on 2026-08-09 they disagreed: two said `f6b0f0d1`, one still
said `024dc9af`, and only reading the worktree's actual `HEAD` settled it. Nothing was broken by
that — CI reads `wheels.yml`, which was right — but the file whose entire job is to record the pin
was the one that was wrong, which is the failure mode worth designing out rather than re-checking.
A value that must not drift should exist once.

| repo | worktree path | branch | commit | purpose |
|---|---|---|---|---|
| mzLib | `code/mzLib` | *(detached)* | see [`mzlib.pin`](mzlib.pin) | source of `PrideArchiveClient`; referenced by `pkg/bridge/MzLibBridge.csproj` |

Recreate the worktree on another machine (or after deleting it) with:

```powershell
git -C E:\GitClones\mzLib worktree add --detach E:\CodeReview\pyMzLib\code\mzLib (Get-Content code\mzlib.pin).Trim()
```

## Why a worktree and not the mzLib NuGet package

The published `mzLib` package (1.0.584 on 2026-08-09) is built from a hand-authored
`.nuspec` that lists DLLs explicitly and declares heavy dependencies — TorchSharp, libtorch-cpu,
Microsoft.ML, SkiaSharp — plus native vendor DLLs (Thermo RawFileReader, timsdata, baf2sql).
Referencing the package would drag all of that into the wheel. A `ProjectReference` to the
narrowest project that carries `PrideArchiveClient` is both smaller and pinned to an exact
commit, which matters when the bridge and the Python package must be built from the same source.

Even so, the transitive graph is heavy: `UsefulProteomicsDatabases → Proteomics → Chromatography → TorchSharp → libtorch`. See gap G7 — this is the main open engineering question, and
the best fix is upstream in mzLib, not here.

Caveat: `rebuild-mzlib-for-metamorpheus` hardcodes `E:\GitClones\mzLib`; building from this
worktree needs the worktree path passed explicitly.

## Why the pin does not simply track mzLib's master

Because a release must be reproducible, and because mzLib's master moves for reasons that have
nothing to do with the bindings. What the pin should track is mzLib's *releases*, and until
2026-08-09 nothing did: mzLib released 1.0.584 on 2026-08-07 and no job in any of the three binding
repositories noticed, because none of them consumes mzLib as a release at all.

`upstream-watch.yml` closes that. It resolves mzLib's latest release tag to a commit weekly and
opens a pull request bumping `mzlib.pin` — a pull request, not a push, because `wheels.yml` already
runs on `pull_request` and its cross-platform matrix is exactly the gate that should decide whether
a new mzLib is safe to build against. The bump arrives reviewed by CI; a human only merges.

## Pin history

Appended by `upstream-watch.yml` on each bump, newest last. Hand-written entries are fine too — the
workflow only ever appends one row to the end of the file, so this table must stay last in the
document and nothing may be added below it.

| date | commit | mzLib release | why |
|---|---|---|---|
| 2026-07-24 | `525cb7c8` | — | initial pin; source of `PrideArchiveClient` |
| 2026-08-03 | `f6b0f0d1` | — | **#1116** (MSFragger retention time converted from seconds to minutes at the reader, which let the bindings retire the MSFragger seconds caveat) and **#1121** (`PrideArchiveClient.GetProjectFilesFromFtpAsync`, the complete PRIDE FTP file list the REST manifest omits). Belonged to no mzLib release at the time; the commit is *"Remove y ions from the ETD and ECD product sets (#1109) (#1114)"* |
| 2026-08-09 | `e3220d6a` | 1.0.584 | tracked automatically by `upstream-watch.yml`; see the pull request for the commits it brought in |
| 2026-08-15 | `5ba13155` | 1.0.585 | **#1141** moved mzLib to .NET 10, so this bump is not a pin move alone — `MzLibBridge.csproj`, its test project and `DOTNET_VERSION` move with it (an 8.0.x SDK cannot build a net10.0 project). Opened by hand rather than by `upstream-watch.yml`, because the watcher would have opened a red pull request it had no way to fix |
