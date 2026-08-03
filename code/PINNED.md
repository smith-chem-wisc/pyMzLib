# Code worktrees + pins — pyMzLib

Worktrees under `code/` are gitignored; this file is their record.

| repo | worktree path | branch | base commit | purpose |
|---|---|---|---|---|
| mzLib | `code/mzLib` | *(detached)* | `f6b0f0d17f32383918ef895006aaecb71cdb9a7e` | source of `PrideArchiveClient`; referenced by `pkg/bridge/MzLibBridge.csproj` |

Pinned at `origin/master` on 2026-08-03 — *"Remove y ions from the ETD and ECD product sets
(#1109) (#1114)"*. Re-pinned from `525cb7c8` to pick up two merged mzLib changes the bindings
depend on: **#1116** (MSFragger retention time converted from seconds to minutes at the reader —
`MsFraggerPsm.RetentionTime = RetentionTimeInSeconds / 60`, which lets this branch retire the
MSFragger seconds caveat honestly) and **#1121** (`PrideArchiveClient.GetProjectFilesFromFtpAsync`
— the complete PRIDE FTP file list the REST manifest omits, so the bridge can project the true
file list rather than the incomplete one).

Recreate it on another machine (or after deleting it) with:

```powershell
git -C E:\GitClones\mzLib worktree add --detach E:\CodeReview\pyMzLib\code\mzLib f6b0f0d1
```

## Why a worktree and not the mzLib NuGet package

The published `mzLib` package (1.0.583 at time of writing) is built from a hand-authored
`.nuspec` that lists DLLs explicitly and declares heavy dependencies — TorchSharp, libtorch-cpu,
Microsoft.ML, SkiaSharp — plus native vendor DLLs (Thermo RawFileReader, timsdata, baf2sql).
Referencing the package would drag all of that into the wheel. A `ProjectReference` to the
narrowest project that carries `PrideArchiveClient` is both smaller and pinned to an exact
commit, which matters when the bridge and the Python package must be built from the same source.

Even so, the transitive graph is heavy: `UsefulProteomicsDatabases → Proteomics → Chromatography → TorchSharp → libtorch`. See gap G7 — this is the main open engineering question, and
the best fix is upstream in mzLib, not here.

Caveat: `rebuild-mzlib-for-metamorpheus` hardcodes `E:\GitClones\mzLib`; building from this
worktree needs the worktree path passed explicitly.
