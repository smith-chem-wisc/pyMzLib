# PRIDE Archive

The [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
repository — the place most published mass-spectrometry data lives. pyMzLib exposes mzLib's
`PrideArchiveClient`, so listing and retrieving a project's files is a couple of lines.

Everything here goes through the same code MetaMorpheus uses in C#: the same paging, the same
URL resolution, the same safe-download behavior.

## Listing a project's files

```python
import pymzlib

files = pymzlib.pride.list_files("PXD000001")
print(len(files))          # 8
```

You get the **complete** manifest. PRIDE serves file lists in pages; that is handled for you, so
a project with 4,000 files returns 4,000 entries in one list, not the first 100.

An unknown accession returns an empty list rather than raising. That is PRIDE's own behavior and
pyMzLib preserves it rather than inventing an error PRIDE didn't report:

```python
pymzlib.pride.list_files("PXD999999999")   # []
```

## What a file tells you

Each entry is a [`PrideFile`](../reference.md):

```python
f = files[0]

f.file_name         # 'TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw'
f.category          # 'RAW'  — also 'PEAK', 'SEARCH', 'OTHER', ...
f.file_size_bytes   # 253434822
f.size_mb           # 253.4
f.extension         # '.raw'
f.checksum          # '' when PRIDE publishes none
f.https_url         # direct download URL, or None
f.downloadable      # False for Aspera-only files
f.submission_date   # datetime, timezone-aware
```

!!! note "Not every file is reachable over HTTPS"
    PRIDE publishes locations by protocol, and a few files are Aspera-only. Those report
    `downloadable == False` and `https_url is None`. Check before assuming a download will
    succeed:

    ```python
    unreachable = [f for f in files if not f.downloadable]
    ```

## Sizing up a project before you commit

Worth doing before pulling a dataset that turns out to be 400 GB:

```python
raw = [f for f in files if f.category == "RAW"]
gb = pymzlib.pride.total_size_bytes(raw) / 1e9
print(f"{len(raw)} raw files, {gb:.1f} GB")
```

## Downloading

```python
paths = pymzlib.pride.download("PXD000001", "downloads")
```

`paths` is a list of `pathlib.Path` objects for what was written. The destination directory is
created if it doesn't exist.

### Downloading only part of a project

Filters combine with AND:

```python
# Just the raw files
pymzlib.pride.download("PXD000001", "downloads", category="RAW")

# Just search results in two formats
pymzlib.pride.download("PXD000001", "downloads", extensions=[".mzid", ".mztab"])

# Raw files that are also .raw (belt and braces)
pymzlib.pride.download("PXD000001", "downloads", category="RAW", extensions=[".raw"])
```

### Resuming

```python
pymzlib.pride.download("PXD000001", "downloads", overwrite=False)
```

With `overwrite=False`, a file already present at the destination is skipped without a request.
Re-running after an interruption picks up where it stopped.

### Interruptions don't corrupt anything

Each file streams to a temporary `.partial` name and is moved into place only once complete. A
cancelled or crashed download leaves no truncated file behind that a later run would mistake for
a good one. This matters more than it sounds: a silently truncated `.raw` produces a search that
runs to completion and gives wrong answers.

## Timeouts

`list_files()` defaults to a 300-second timeout. `download()` has **no** timeout by default,
because multi-gigabyte transfers legitimately take hours:

```python
pymzlib.pride.list_files("PXD000001", timeout=60)      # give up after a minute
pymzlib.pride.download("PXD000001", "out", timeout=3600)  # cap it at an hour
```

## A worked example

Fetch every raw file from a project, but only if the total is manageable:

```python
import pymzlib

ACCESSION = "PXD000001"
LIMIT_GB = 5

files = pymzlib.pride.list_files(ACCESSION)
raw = [f for f in files if f.category == "RAW" and f.downloadable]
size_gb = pymzlib.pride.total_size_bytes(raw) / 1e9

if size_gb > LIMIT_GB:
    raise SystemExit(f"{size_gb:.1f} GB exceeds the {LIMIT_GB} GB limit; refine the filter first")

print(f"Downloading {len(raw)} files ({size_gb:.1f} GB)…")
paths = pymzlib.pride.download(ACCESSION, f"data/{ACCESSION}", category="RAW", overwrite=False)
print(f"Wrote {len(paths)} files")
```

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | The accession is blank, or `page_size` isn't positive. Raised before any network call. |
| `BridgeError` with `error_type='HttpRequestException'` | PRIDE returned a failure status, or was unreachable. |
| `BridgeError` with `error_type='NotSupportedException'` | A selected file has no HTTPS location (Aspera-only). Filter on `downloadable` first. |

```python
try:
    pymzlib.pride.download("PXD000001", "downloads")
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```
