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

Paging is handled for you: PRIDE serves file lists in pages, so a project with 4,000 files returns
4,000 entries in one list, not the first 100.

!!! warning "This is PRIDE's REST manifest, which is not always the whole project"
    `list_files()` returns what PRIDE's REST API publishes — and that is sometimes incomplete. For
    PXD000001 it returns **8** files while the project's FTP tree holds **13**, omitting the two
    largest (the ~450 MB `.mzML` and `.mzXML` conversions). When completeness matters, use
    [`list_ftp_files()`](#the-complete-file-list-from-the-ftp-tree) instead.

Accessions are case-insensitive and whitespace is trimmed — `"pxd000001"` works — but an accession
that doesn't exist **raises** rather than returning nothing:

```python
pymzlib.pride.list_files("PXD999999999")   # ProjectNotFoundError
pymzlib.pride.list_files("banana")         # UsageError — not a valid accession at all
```

PRIDE itself answers an unknown accession with an empty result rather than a 404, and early
versions of pyMzLib passed that straight through. That was a mistake: an empty list is
indistinguishable from "this project genuinely has no matching files", so a typo produced a script
that reported *0 files, done* and carried on. A wrong answer that looks like a right answer is
worse than an error.

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

!!! warning "`total_size_bytes()` has two blind spots"
    It sums over the REST manifest — which is **incomplete** — and PRIDE frequently reports the
    *decompressed* size for `.gz` files. For PXD000001 it returns 0.51 GB where the project on disk
    is 1.44 GB. The two errors run in opposite directions and do not cancel. For a size that covers
    the whole project, use `approximate_total_size_bytes()` over the FTP listing below.

## The complete file list (from the FTP tree)

When you need to know *everything* a project holds — not just what its REST manifest lists —
`list_ftp_files()` walks the project's FTP directory tree and returns the authoritative list,
subdirectories included:

```python
ftp = pymzlib.pride.list_ftp_files("PXD000001")
print(len(ftp))          # 13, not the 8 the REST manifest reports
```

Each entry is a [`PrideFtpFile`](../reference.md):

```python
f = ftp[0]

f.relative_path          # 'run1.raw', or 'generated/summary.mztab' for a nested file
f.file_name              # 'summary.mztab' — the bare leaf name
f.url                    # the HTTPS URL to download from
f.approximate_size_bytes # PRIDE's rounded index size (see below)
f.approximate_size_mb    # the same, in MB
f.extension              # '.raw'
```

The whole project's size, covering every file this time:

```python
gb = pymzlib.pride.approximate_total_size_bytes(ftp) / 1e9
print(f"{len(ftp)} files, ~{gb:.2f} GB")     # ~1.44 GB
```

!!! note "Why *approximate*"
    PRIDE's FTP directory index rounds each size to about three significant figures (its `429M`
    becomes 449,839,104 bytes), so these are good for a project-size estimate but are **not** exact.
    When you need the precise transfer size of one file — to budget a download against — issue an
    HTTP `HEAD` against its `url` and read the `Content-Length`.

**Which listing should you use?** Reach for `list_files()` when you want the rich per-file metadata
(category, checksum, controlled-vocabulary locations) the REST manifest carries — it is what
`download()` and `download_files()` select from. Reach for `list_ftp_files()` when completeness or a
true project size matters, and the REST manifest might be hiding files. An unknown accession
**raises** here too — mzLib resolves the project before walking, so a typo fails loudly rather than
returning an empty list.

## Downloading what you selected

This is usually the one you want. Filter the manifest with the full expressiveness of Python, then
hand the result straight back:

```python
files = pymzlib.pride.list_files("PXD000001")
small = [f for f in files if f.size_mb < 5 and f.downloadable]

pymzlib.pride.download_files(small, "downloads")
```

`download()`'s `category` and `extensions` filters can only say what they were built to say.
*"Under 5 MB"*, *"the three most recent"*, or *"everything except the MGF"* cannot be expressed in
that vocabulary at all — and they are all one list comprehension away.

This matters more than it sounds on real projects: PXD000001 publishes `.mztab.gz`, `.mgf.gz` and
`.xml.gz`, which all have `extension == ".gz"`. There is no extension filter that selects the
mzTab without also taking the 16 MB MGF. There is an obvious list comprehension that does.

## Downloading a whole project

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
| `ServiceUnavailableError` | PRIDE was down, rate-limiting, or timed out (HTTP 408/429/5xx). **Not your bug** — retry later. |
| `BridgeError` with `error_type='HttpRequestException'` | PRIDE answered with a client error such as 404. Something about the request is wrong. |
| `BridgeError` with `error_type='NotSupportedException'` | A selected file has no HTTPS location (Aspera-only). Filter on `downloadable` first. |

`ServiceUnavailableError` is a subclass of `BridgeError`, so catching `BridgeError` still catches
everything. Separating them lets you retry the failures worth retrying and report the rest:

```python
import time

for attempt in range(3):
    try:
        files = pymzlib.pride.list_files("PXD000001")
        break
    except pymzlib.ServiceUnavailableError:
        time.sleep(30)      # EBI's problem; it may well pass
else:
    raise SystemExit("PRIDE unavailable after three attempts")
```

```python
try:
    pymzlib.pride.download("PXD000001", "downloads")
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```
