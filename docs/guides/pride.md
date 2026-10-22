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
| Check a download against the MD5 and exact size PRIDE recorded | [`file_checksums()`](#checking-a-download-against-prides-checksums) | `PrideArchiveClient.GetFileChecksumsAsync` |

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
project before walking it, so a typo fails loudly rather than returning an empty list. To check a
file you downloaded, use the third list,
[`file_checksums()`](#checking-a-download-against-prides-checksums).

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
re-running after an interruption picks up where it stopped. The file it skips is not checked, so a
damaged copy stays damaged. [Check it](#checking-a-download-against-prides-checksums) first.

### Interruptions don't corrupt anything

Each file streams to a temporary `.partial` name and moves into place only once complete. A
cancelled or crashed download leaves no truncated file under its real name for a later run to
mistake for a good one. That matters: a silently truncated `.raw` produces a search that runs to
completion and gives wrong answers.

A transfer that fails after mzLib's retries leaves `<name>.partial` and `<name>.partial.validator`
in the folder on purpose: the next call for that file resumes from them instead of starting again
(mzLib#1412). Leave them in place, or delete both to start over.

### Check large files yourself

mzLib can check a download against PRIDE's checksums, and since 1.0.595 it also catches an EBI
error page served at the right size (mzLib#1427, #1449). That check is not yet on the path
`download()` uses, so a file can be reported as downloaded when its bytes are wrong. For files
that matter, check each one against `file_checksums()` before you search it, as the next section
shows.

## Checking a download against PRIDE's checksums

**Is the file on your disk the file PRIDE holds?** `download()` cannot tell you yet. PRIDE records
an MD5 and the exact size of every file a submitter deposits, and `file_checksums()` returns that
record:

```pycon
>>> checksums = pymzlib.pride.file_checksums("PXD000001")
>>> len(checksums)
12
>>> raw = checksums["TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw"]
>>> raw.md5, raw.size_bytes
('c567cecbb9866796e690a4958746f7b6', 220475548)

```

It is a dict keyed by file name, so look a file up by the name it was downloaded under. Each value
is a [`PrideFileChecksum`][pymzlib.pride.PrideFileChecksum]:

- `md5` is 32 hexadecimal characters, as PRIDE wrote them. Compare it without regard to case.
- `size_bytes` is the exact length a download of the file comes to. Neither other listing gives
  you that. `list_ftp_files()` rounds its sizes, and `list_files()` can be wrong:

```pycon
>>> rest = next(f for f in files if f.file_name == "PRIDE_Exp_Complete_Ac_22134.xml.gz")
>>> rest.file_size_bytes, checksums[rest.file_name].size_bytes
(10677205, 10668000)

```

The download is 10,668,000 bytes, the checksum list's number. Because each size is exact, so is
their sum: what the deposited files weigh, which the other two listings can only estimate. It
leaves out the `generated/` files, which have no row (below).

```pycon
>>> round(sum(c.size_bytes for c in checksums.values()) / 1e9, 2)
1.43

```

### The check

A file is PRIDE's when its length and its MD5 both match. The length costs nothing to check. The
MD5 reads the whole file, which takes a few seconds per gigabyte. Check both, because the length
alone does not prove the bytes: on 7 and 8 October 2026, EBI served about 23 of about 70 large
transfers at the right size with its own error text written over part of them (mzLib#1449).

```pycon
>>> import hashlib
>>> from pathlib import Path
>>> def check(path, expected):
...     """None when the file is the one PRIDE recorded, otherwise what is wrong with it."""
...     size = path.stat().st_size
...     if size != expected.size_bytes:
...         return f"{size} bytes where PRIDE records {expected.size_bytes}"
...     md5 = hashlib.md5()
...     with path.open("rb") as handle:
...         for block in iter(lambda: handle.read(1 << 20), b""):
...             md5.update(block)
...     if md5.hexdigest() != expected.md5.lower():
...         return f"MD5 {md5.hexdigest()} where PRIDE records {expected.md5}"
...     return None

```

Here it is on two damaged copies of the project's 1,645-byte `README.txt`. The first was cut short.
The second is the right size with EBI's error text inside, and only the MD5 catches it:

```pycon
>>> readme = checksums["README.txt"]
>>> damaged = Path("README.txt")
>>> _ = damaged.write_bytes(b"x" * 1000)
>>> check(damaged, readme)
'1000 bytes where PRIDE records 1645'
>>> _ = damaged.write_bytes(b"<Error><Code>ConnectionClosedException</Code></Error>".ljust(1645, b"."))
>>> check(damaged, readme)
'MD5 6a8f1a7d39089d34bf9bc9b19e553dd4 where PRIDE records 492753137013825025e6729c8798664f'

```

After a real download, check every file it wrote:

```python title="Not run: downloads files from PRIDE"
paths = pymzlib.pride.download("PXD000001", "downloads", category="RAW")
checksums = pymzlib.pride.file_checksums("PXD000001")
for path in paths:
    expected = checksums.get(path.name)
    if expected is None:
        print(f"{path.name}: PRIDE has no checksum for it, so it cannot be checked")
        continue
    problem = check(path, expected)
    if problem:
        print(f"{path.name}: {problem}. Download it again.")
```

To download a damaged file again, use `download_files()` with that file, or `download()` with
`overwrite=True`, the default. `overwrite=False` would skip it, because it is already there.

### What the list covers

It covers the files the submitter deposited at the project's root. It is a third list, not a view
of either of the others, and three things follow from that.

**PRIDE's own `generated/` files have no row.** A file without a row cannot be checked; that does
not make it invalid. In PXD000001 those are the two files under `generated/`:

```pycon
>>> [f.relative_path for f in ftp if f.relative_path not in checksums]   # doctest: +NORMALIZE_WHITESPACE
['generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz',
 'generated/PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz']

```

`list_files()` names those two without their folder, so the same two are the manifest's files
without a row:

```pycon
>>> [f.file_name for f in files if f.file_name not in checksums]
['PRIDE_Exp_Complete_Ac_22134.pride.mztab.gz', 'PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz']

```

**It names files the REST manifest omits,** such as the 450 MB `.mzML` that only the FTP listing
showed above:

```pycon
>>> "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML" in checksums
True

```

**It can name files that are not deposits at all.** PXD015239's list includes a 2.3 GB NFS temp
file, `.nfs80360b91005fae5300002478`. Look up the files you have by name; do not read the list as
the project's contents.

### Two traps

**`PrideFile.checksum` is not this MD5.** Where PRIDE fills it in, it is a 40-character SHA-1, a
different digest of the same bytes. For every PXD000001 file it is empty:

```pycon
>>> {f.checksum for f in files}
{''}

```

**An empty dict is an answer.** It means the project exists but PRIDE holds no checksums for it,
so nothing can be checked. An accession PRIDE does not know raises `ProjectNotFoundError`, and so
does a project not yet published: PRIDE holds no checksums before publication.

## Retries and timeouts

**mzLib retries for you.** A request that fails in a way that can pass (no answer, an HTTP 408,
429 or 5xx, or a 403 from EBI's FTP host, its rate limit) is retried three times, after 5, 20 and
60 seconds, and that request alone (mzLib#1412). A paged manifest retries only the page that
failed, and a download resumes where it stopped.

**pyMzLib adds no timeout of its own by default**, on any PRIDE function. Each attempt is already
limited to 100 seconds, so an outage that hangs is reported as `ServiceUnavailableError` after
about 8 minutes per request. Pass `timeout=` to cut a call off sooner; it then raises
`BridgeTimeoutError`, which is not a `ServiceUnavailableError`:

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
| `ProjectNotFoundError` | PRIDE has no project with that accession, or it has no files. From `file_checksums()`, also a project not yet published. |
| `ServiceUnavailableError` | PRIDE was down, rate-limiting or timed out, **or the connection dropped part-way through a transfer**. **Not your bug**: retry later. |
| `BridgeError` with `error_type='HttpRequestException'` | PRIDE answered with a client error such as 404, so something about the request is wrong. **Two exceptions** look the same but are EBI failures: an FTP listing that came back empty, and a 403 from EBI's FTP host, each still failing after mzLib's retries. Treat those as an outage. |
| `BridgeError` with `error_type='NotSupportedException'` | A selected file has no HTTPS location (Aspera-only). Filter on `downloadable` first. |
| `BridgeError` with `error_type='MzLibException'` from `file_checksums()` | PRIDE sent a checksum list mzLib cannot read: a new header, a malformed row, or a name listed twice. Not your bug; mzLib does not retry it, so report it. |

Since mzLib#1350, every PRIDE transport failure reaches the bridge as an `HttpRequestException`,
carrying its HTTP status where there was one, and no message carries a URL. A connection that
drops mid-download is still an outage, so it still raises `ServiceUnavailableError`.

!!! note "A download that dies in the middle is an outage, not a bug"
    A request that fails outright reports a status code. One that fails **after** the response has
    begun has no status left to report, because the server already said 200, so it surfaces as a
    truncated stream. Both are the same thing from your side, and both raise
    `ServiceUnavailableError`, so the handling below covers `download()` as well as `list_files()`.

`ServiceUnavailableError` is a subclass of `BridgeError`, so catching `BridgeError` catches
everything. By the time it reaches you, mzLib has already retried, so retrying again at once
rarely helps. Wait minutes, not seconds:

```python title="Not run: shows how to handle a PRIDE outage, which a recording cannot produce"
import time

try:
    files = pymzlib.pride.list_files("PXD000001")
except pymzlib.ServiceUnavailableError:
    time.sleep(600)         # mzLib retried already; give EBI ten minutes
    files = pymzlib.pride.list_files("PXD000001")

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
