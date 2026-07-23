# Code worktrees + pins — pyMzLib

Worktrees under `code/` are gitignored; this file is their record.

| repo | worktree path | branch | base commit | purpose |
|---|---|---|---|---|
| mzLib | `code/mzLib` | *(detached)* | `525cb7c8e3875501868fb961d0bc1bd34ef5174e` | source of `PrideArchiveClient`; referenced by `pkg/bridge/MzLibBridge.csproj` |

Pinned at `origin/master` on 2026-07-23 — *"fix(variants): read UniProt-native deletions (empty
`<variation>`) (#1095)"*.

Recreate it on another machine (or after deleting it) with:

```powershell
git -C E:\GitClones\mzLib worktree add --detach E:\CodeReview\pyMzLib\code\mzLib 525cb7c8
```

## Why a worktree and not the mzLib NuGet package

The published `mzLib` package (1.0.583 at time of writing) is built from a hand-authored
`.nuspec` that lists DLLs explicitly and declares heavy dependencies — TorchSharp, libtorch-cpu,
Microsoft.ML, SkiaSharp — plus native vendor DLLs (Thermo RawFileReader, timsdata, baf2sql).
Referencing the package would drag all of that into the wheel. A `ProjectReference` to the
narrowest project that carries `PrideArchiveClient` is both smaller and pinned to an exact
commit, which matters when the bridge and the Python package must be built from the same source.

Even so, the transitive graph is heavy: `UsefulProteomicsDatabases → Proteomics/Transcriptomics
→ Omics → TorchSharp → libtorch`. See gap G7 — this is the main open engineering question, and
the best fix is upstream in mzLib, not here.

Caveat: `rebuild-mzlib-for-metamorpheus` hardcodes `E:\GitClones\mzLib`; building from this
worktree needs the worktree path passed explicitly.
