# PRIDE Archive

The [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/) is EBI's public proteomics data
repository, where most published mass-spectrometry data lives. pyMzLib exposes mzLib's
`PrideArchiveClient`, so finding a dataset, listing its files and downloading them takes a few
lines. Everything goes through the code MetaMorpheus uses: the same paging, the same URL
resolution, the same safe-download behaviour.

| You want to | Call | mzLib does it with |
|---|---|---|
| Find datasets about a subject | [`search()`](#finding-projects-in-the-first-place) | `PrideArchiveClient.SearchProjectsAsync` |
| See a project's files, with categories and checksums | [`list_files()`](#listing-a-projects-files) | `PrideArchiveClient.GetProjectFilesAsync` |
| See **every** file, including ones the manifest omits | [`list_ftp_files()`](#the-complete-file-list-from-the-ftp-tree) | `PrideArchiveClient.GetProjectFilesFromFtpAsync` |
| Size a project before downloading it | [`approximate_total_size_bytes()`](#sizing-up-a-project-before-you-commit) | `TotalApproximateSizeBytes` |
| Download the files you picked | [`download_files()`](#downloading-what-you-selected) | `PrideArchiveClient.DownloadFilesAsync` |
| Download by category or extension | [`download()`](#downloading-a-whole-project) | `WhereCategory`, `WhereExtension` |

Every `>>>` example on this page runs in CI against answers recorded from the live PRIDE Archive
through the real bridge. Blocks titled **Not run** download files or show retry logic, and say so.

## Finding projects in the first place

Everything else on this page takes an accession you already have. `search()` produces them, so you
can go from a subject to a dataset without leaving Python:

```pycon
>>> import pymzlib
>>> hits = pymzlib.pride.search("plasmodium falciparum schizont")
>>> len(hits)
6
>>> hits[0].accession
'PXD070842'
>>> hits[0].organisms
['Homo sapiens (human)', 'Plasmodium falciparum (isolate 3d7)']

```

Paging is handled for you, and no accession is repeated. From a hit, the rest of this page follows:
`pymzlib.pride.list_files(hits[0].accession)`.

### Why it matched

`highlights` is the one thing search returns that the project metadata cannot: the snippets that
matched, keyed by the field each was found in.

```pycon
>>> hits[0].matched_fields
['references', 'title']
>>> hits[0].highlights["title"][0][:60]
'High-resolution spatial proteomics of <em>Plasmodium</em> <e'

```

The keys vary per hit and per query, and the `<em>` markup is PRIDE's. Strip it if the value is
going anywhere other than a highlighted view.

### A search hit is not a project's metadata

!!! warning "Controlled-vocabulary fields arrive flattened to display strings"
    PRIDE serves search from a **separate Elasticsearch projection**. The same project reports its
    instruments as `["Q Exactive"]` here and as structured terms with accessions from the metadata
    endpoint. Contacts collapse from ten-field objects to a display name, and publications to one
    pre-formatted citation string.

    That is a property of PRIDE's wire, not a simplification pyMzLib chose: **the accessions are
    not sent.** Resolving a display name against a vocabulary to manufacture one would hand you an
    identifier PRIDE never asserted. Follow the hit's `accession` when you need the vocabulary.

Two more consequences to know before you index anything:

- **`project_file_names` is not the manifest.** It carries names only, with no sizes, categories or
  download locations. Use `list_files()` or `list_ftp_files()` to act on files.
- **`sdrf` is not a file.** It is the project's SDRF metadata flattened by the search index into one
  space-joined bag of term values, with the rows and columns gone. Nothing can be fetched with it.
  For a real SDRF, see the [SDRF guide](sdrf.md).

### Zero means "not reported"

PRIDE omits nothing as `null`, so an absent value arrives as `0`, `""` or `[]`. Several fields are
sparse: across a sample of hits, `project_tags`, `sdrf` and the bot/hub/organic download split
were each empty on most. **A `download_count` of 0 does not mean nobody downloaded it.**

The same applies to `keywords`. PRIDE ships empty and whitespace-only strings inside it, and pyMzLib
passes them through rather than filtering, so that it agrees with mzLib and with the Rust and R
bindings about what a project's keywords are. This hit's only keyword is an empty string. Filter
before you join them:

```pycon
>>> hits[0].keywords
['']
>>> [k for k in hits[0].keywords if k.strip()]
[]

```

### Dates here are calendar dates, not timestamps

A search hit's `submission_date`, `publication_date` and `updated_date` are `datetime.date`, where a
`PrideFile`'s are `datetime`. That follows the wire: the search endpoint sends a date with no time
and no offset, so presenting it as a timestamp would attach a midnight PRIDE never reported.

```pycon
>>> hits[0].submission_date
datetime.date(2025, 11, 17)

```

!!! note "A live index has no stable cursor"
    PRIDE pages search results from a live index. A result set that changes *during* a multi-page
    fetch shifts its own paging: a project published mid-fetch is served on two pages and
    deduplicated, so it comes back once, but a project *removed* mid-fetch can fall between two
    pages and be missed. A search whose hits fit on one page cannot be affected.

## Listing a project's files

```pycon
>>> files = pymzlib.pride.list_files("PXD000001")
>>> len(files)
8

```

Paging is handled for you: PRIDE serves file lists in pages, and you get one list however many
pages the project spans.

!!! warning "This is PRIDE's REST manifest, which is not always the whole project"
    `list_files()` returns what PRIDE's REST API publishes, and that is sometimes incomplete. For
    PXD000001 the FTP tree holds files the manifest omits, including the two largest (the
    ~450 MB `.mzML` and `.mzXML` conversions). When completeness matters, use
    [`list_ftp_files()`](#the-complete-file-list-from-the-ftp-tree).

Accessions are case-insensitive and whitespace is trimmed, so `"pxd000001"` works. But an accession
that does not exist **raises** rather than returning nothing:

```pycon
>>> pymzlib.pride.list_files("PXD999999999")               # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib.pride.ProjectNotFoundError: ...
>>> pymzlib.pride.list_files("banana")                     # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
pymzlib._bridge.UsageError: ...

```

PRIDE itself answers an unknown accession with an empty result rather than a 404, and early
versions of pyMzLib passed that straight through. That was a mistake: an empty list cannot be told
apart from "this project has no matching files", so a typo produced a script that reported *0
files, done* and carried on. A wrong answer that looks like a right answer is worse than an error.

## What a file tells you

Each entry is a [`PrideFile`][pymzlib.pride.PrideFile]. The manifest is in repository order, so pick
files by what they are, not by position:

```pycon
>>> raw = next(f for f in files if f.category == "RAW")
>>> raw.file_name
'TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw'
>>> raw.file_size_bytes, round(raw.size_mb, 1), raw.extension
(220475548, 220.5, '.raw')
>>> raw.checksum                       # '' when PRIDE publishes none
''
>>> raw.downloadable
True
>>> sorted({f.category for f in files})
['OTHER', 'PEAK', 'RAW', 'RESULT', 'SEARCH']

```

!!! note "Not every file is reachable over HTTPS"
    PRIDE publishes locations by protocol, and a few files are Aspera-only. Those report
    `downloadable == False` and `https_url is None`. Every PXD000001 file is reachable:

    ```pycon
    >>> [f.file_name for f in files if not f.downloadable]
    []

    ```

## Sizing up a project before you commit

Worth doing before you pull a dataset that turns out to be 400 GB:

```pycon
>>> raws = [f for f in files if f.category == "RAW"]
>>> f"{len(raws)} raw file, {pymzlib.pride.total_size_bytes(raws) / 1e9:.2f} GB"
'1 raw file, 0.22 GB'

```

!!! warning "`total_size_bytes()` has two blind spots"
    It sums over the REST manifest, which is **incomplete**, and PRIDE often reports the
    *decompressed* size for `.gz` files. The two errors run in opposite directions and do not
    cancel. Compare the manifest's total with the FTP listing's, below:

    ```pycon
    >>> round(pymzlib.pride.total_size_bytes(files) / 1e9, 2)
    0.51

    ```

## The complete file list (from the FTP tree)

When you need *everything* a project holds, not only what its REST manifest lists,
`list_ftp_files()` walks the project's FTP directory tree, subdirectories included:

```pycon
>>> ftp = pymzlib.pride.list_ftp_files("PXD000001")
>>> len(ftp), len(files)
(14, 8)
>>> round(pymzlib.pride.approximate_total_size_bytes(ftp) / 1e9, 2)
1.44

```

The FTP tree holds almost three times the bytes the manifest admits to. The files only it lists
are the ones most people want, the modern open-format conversions:

```pycon
>>> manifest = {f.file_name for f in files}
>>> [f.relative_path for f in ftp if f.file_name not in manifest and f.extension == ".mzml"]
['TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML']

```

Each entry is a [`PrideFtpFile`][pymzlib.pride.PrideFtpFile]. A file in a subdirectory keeps its
path, and `file_name` is the last segment:

```pycon
>>> nested = next(f for f in ftp if "/" in f.relative_path)
>>> nested.relative_path, nested.file_name
('generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz', 'PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz')
>>> nested.url.startswith("https://ftp.pride.ebi.ac.uk/")
True

```

!!! note "Why *approximate*"
    PRIDE's FTP directory index rounds each size to about three significant figures (its `429M`
    becomes 449,839,104 bytes), so these sizes are good for a project estimate but are **not**
    exact. For the precise transfer size of one file, issue an HTTP `HEAD` against its `url` and
    read the `Content-Length`.

**Which listing should you use?** Use `list_files()` when you want the per-file metadata the REST
manifest carries (category, checksum, controlled-vocabulary locations); it is what `download()` and
`download_files()` select from. Use `list_ftp_files()` when completeness or a true project size
matters. An unknown accession **raises** `ProjectNotFoundError` there too: mzLib resolves the
project before walking it, so a typo fails loudly rather than returning an empty list.

!!! warning "Downloading a file that only the FTP listing found"
    `download()` and `download_files()` work on the **REST manifest**, so a file that appears *only*
    in `list_ftp_files()` cannot be passed to them. Fetch it from its `url` with a plain HTTPS GET:

    ```python title="Not run: downloads a 450 MB file from PRIDE"
    import urllib.request

    hidden = next(f for f in ftp if f.extension == ".mzml")
    urllib.request.urlretrieve(hidden.url, hidden.file_name)
    ```

    The field names differ by type: a `PrideFtpFile` has `url`, always an HTTPS location, while a
    `PrideFile` has `https_url`, which is `None` for Aspera-only files.

## Downloading what you selected

This is usually the one you want. Filter the manifest with all of Python, then hand the result
back:

```pycon
>>> small = [f for f in files if f.size_mb < 5 and f.downloadable]
>>> sorted(f.file_name for f in small)
['F063721.dat-mztab.txt', 'PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz', 'erwinia_carotovora.fasta']

```

```python title="Not run: downloads files from PRIDE"
pymzlib.pride.download_files(small, "downloads")
```

`download()`'s `category` and `extensions` filters can say only what they were built to say.
*"Under 5 MB"*, *"the three most recent"* or *"everything except the MGF"* cannot be expressed in
that vocabulary, and each is one list comprehension away.

That matters on real projects. PXD000001 publishes an `.mztab.gz`, an `.mgf.gz` and an `.xml.gz`,
which all have the extension `.gz`, so no extension filter selects the mzTab without also taking
the MGF:

```pycon
>>> sorted(f.file_name for f in files if f.extension == ".gz")
['PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz', 'PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz', 'PRIDE_Exp_Complete_Ac_22134.xml.gz']

```

## Downloading a whole project

```python title="Not run: downloads files from PRIDE"
paths = pymzlib.pride.download("PXD000001", "downloads")
```

`paths` is a list of `pathlib.Path` objects, one per file written. The destination directory is
created if it does not exist.

### Downloading only part of a project

Filters combine with AND:

```python title="Not run: downloads files from PRIDE"
# Just the raw files
pymzlib.pride.download("PXD000001", "downloads", category="RAW")

# Just search results in two formats
pymzlib.pride.download("PXD000001", "downloads", extensions=[".mzid", ".mztab"])

# Raw files that are also .raw
pymzlib.pride.download("PXD000001", "downloads", category="RAW", extensions=[".raw"])
```

A filter that matches nothing raises `UsageError` rather than returning an empty list, for the same
reason an unknown accession does.

### Resuming

```python title="Not run: downloads files from PRIDE"
pymzlib.pride.download("PXD000001", "downloads", overwrite=False)
```

With `overwrite=False`, a file already present at the destination is skipped without a request, so
re-running after an interruption picks up where it stopped.

### Interruptions don't corrupt anything

Each file streams to a temporary `.partial` name and moves into place only once complete. A
cancelled or crashed download leaves no truncated file for a later run to mistake for a good one.
That matters: a silently truncated `.raw` produces a search that runs to completion and gives wrong
answers.

## Timeouts

`list_files()` and `search()` default to a 300-second timeout. `download()` has **no** timeout by
default, because multi-gigabyte transfers can take hours:

```python title="Not run: contacts PRIDE"
pymzlib.pride.list_files("PXD000001", timeout=60)          # give up after a minute
pymzlib.pride.download("PXD000001", "out", timeout=3600)   # cap it at an hour
```

## A worked example

Fetch every raw file from a project, but only if the total is manageable:

```python title="Not run: downloads files from PRIDE"
import pymzlib

ACCESSION = "PXD000001"
LIMIT_GB = 5

files = pymzlib.pride.list_files(ACCESSION)
raw = [f for f in files if f.category == "RAW" and f.downloadable]
size_gb = pymzlib.pride.total_size_bytes(raw) / 1e9

if size_gb > LIMIT_GB:
    raise SystemExit(f"{size_gb:.1f} GB exceeds the {LIMIT_GB} GB limit; refine the filter first")

print(f"Downloading {len(raw)} files ({size_gb:.1f} GB)")
paths = pymzlib.pride.download(ACCESSION, f"data/{ACCESSION}", category="RAW", overwrite=False)
print(f"Wrote {len(paths)} files")
```

## Errors you might hit

| Exception | Means |
|---|---|
| `UsageError` | The accession is blank or malformed, `page_size` is not positive, or a filter matched nothing. Raised before any network call where it can be. |
| `ProjectNotFoundError` | PRIDE has no project with that accession, or it has no files. |
| `ServiceUnavailableError` | PRIDE was down, rate-limiting or timed out, **or the connection dropped part-way through a transfer**. **Not your bug**: retry later. |
| `BridgeError` with `error_type='HttpRequestException'` | PRIDE answered with a client error such as 404. Something about the request is wrong. |
| `BridgeError` with `error_type='NotSupportedException'` | A selected file has no HTTPS location (Aspera-only). Filter on `downloadable` first. |

Since mzLib#1350, every PRIDE transport failure reaches the bridge as an `HttpRequestException`,
carrying its HTTP status where there was one, and no message carries a URL. A connection that
drops mid-download is still an outage, so it still raises `ServiceUnavailableError`.

!!! note "A download that dies in the middle is an outage, not a bug"
    A request that fails outright reports a status code. One that fails **after** the response has
    begun has no status left to report, because the server already said 200, so it surfaces as a
    truncated stream. Both are the same thing from your side, and both raise
    `ServiceUnavailableError`, so the retry loop below covers `download()` as well as `list_files()`.

`ServiceUnavailableError` is a subclass of `BridgeError`, so catching `BridgeError` catches
everything. Separating them lets you retry the failures worth retrying and report the rest:

```python title="Not run: shows how to handle a PRIDE outage, which a recording cannot produce"
import time

for attempt in range(3):
    try:
        files = pymzlib.pride.list_files("PXD000001")
        break
    except pymzlib.ServiceUnavailableError:
        time.sleep(30)      # EBI's problem; it may well pass
else:
    raise SystemExit("PRIDE unavailable after three attempts")

try:
    pymzlib.pride.download("PXD000001", "downloads")
except pymzlib.BridgeError as e:
    print(f"{e.error_type}: {e}")
```

[Every error each verb can raise](../errors.md) is listed on one page.

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
repository the data came from:

--8<-- "docs/reference/_generated/cite.pride.md"
