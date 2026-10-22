using System.Diagnostics.CodeAnalysis;
using System.Net;
using System.Text;
using System.Text.Json;
using UsefulProteomicsDatabases;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the verb handlers, driven through a stub HTTP handler rather than the live service.
/// </summary>
/// <remarks>
/// The stub-handler approach is borrowed directly from mzLib's own <c>PrideArchiveClientTests</c>:
/// <see cref="PrideArchiveClient"/> accepts an injected <see cref="HttpClient"/> precisely so its
/// behavior can be pinned without a network. The bridge now obtains its client through
/// <see cref="Program.PrideClientFactory"/> for the same reason.
/// <para>
/// These cover what the live canary cannot: paging, empty results, and failure classification are
/// all conditions EBI will not produce on demand.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class VerbHandlerTests
{
    /// <summary>An HttpMessageHandler that returns a caller-supplied response and records requests.</summary>
    private sealed class StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) : HttpMessageHandler
    {
        public List<string> RequestedUris { get; } = [];

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            RequestedUris.Add(request.RequestUri!.ToString());
            return Task.FromResult(responder(request));
        }
    }

    private static HttpResponseMessage Json(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };

    private static string FileJson(string fileName, string category = "RAW") =>
        $$"""
        {
          "fileName": "{{fileName}}",
          "fileSizeBytes": 1024,
          "fileCategory": { "cvLabel": "PRIDE", "accession": "PRIDE:0000404", "name": "category", "value": "{{category}}" },
          "publicFileLocations": [
            { "cvLabel": "PRIDE", "accession": "PRIDE:0000469", "name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/x/{{fileName}}" }
          ],
          "submissionDate": "2019-01-15T09:42:57.000+00:00"
        }
        """;

    private static void UseStub(Func<HttpRequestMessage, HttpResponseMessage> responder) =>
        Program.PrideClientFactory = () => new PrideArchiveClient(new HttpClient(new StubHandler(responder))
        {
            BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress),
        });

    [TearDown]
    public void RestoreTheRealClient() =>
        Program.PrideClientFactory = () => new PrideArchiveClient();

    // ---- pride files ---------------------------------------------------------

    [Test]
    public async Task PrideFiles_ReportsCountAndTotalSize()
    {
        UseStub(_ => Json($"[{FileJson("a.raw")},{FileJson("b.raw")}]"));

        JsonElement data = await InvokeAsync("pride", "files", "--accession", "PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("accession").GetString(), Is.EqualTo("PXD012345"));
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("total_size_bytes").GetInt64(), Is.EqualTo(2048));
            Assert.That(data.GetProperty("files").GetArrayLength(), Is.EqualTo(2));
        });
    }

    [Test]
    public async Task PrideFiles_UnknownAccession_ReportsAnEmptyManifestNotAnError()
    {
        // PRIDE answers an unknown accession with an empty result rather than a 404. Inventing an
        // error here would misrepresent what the repository actually said.
        UseStub(_ => Json("[]"));

        JsonElement data = await InvokeAsync("pride", "files", "--accession", "PXDBOGUS");

        Assert.That(data.GetProperty("file_count").GetInt32(), Is.Zero);
    }

    [Test]
    public async Task PrideFiles_PageSizeReachesTheRequestUri()
    {
        string? seen = null;
        Program.PrideClientFactory = () =>
        {
            var handler = new StubHandler(request => { seen = request.RequestUri!.ToString(); return Json("[]"); });
            return new PrideArchiveClient(new HttpClient(handler) { BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress) });
        };

        await InvokeAsync("pride", "files", "--accession", "PXD012345", "--page-size", "7");

        Assert.That(seen, Does.Contain("pageSize=7"));
    }

    [TestCase("GetProjectFilesAsync", "files", "--accession", "PXD012345")]
    [TestCase("SearchProjectsAsync", "search", "--keyword", "phospho")]
    public async Task Pride_WithoutAPageSize_SendsMzLibsDefault(string method, string verb, string option, string value)
    {
        // The default is mzLib's, read from its signature here, never written down: the bridge once
        // passed a 100 of its own to both.
        object mzLibDefault = typeof(PrideArchiveClient).GetMethod(method)!.GetParameters()
            .Single(p => p.Name == "pageSize").DefaultValue!;
        string? seen = null;
        Program.PrideClientFactory = () =>
        {
            var handler = new StubHandler(request => { seen = request.RequestUri!.AbsoluteUri; return Json("[]"); });
            return new PrideArchiveClient(new HttpClient(handler) { BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress) });
        };

        await InvokeAsync("pride", verb, option, value);

        Assert.That(seen, Does.Contain($"pageSize={mzLibDefault}"));
    }

    [TestCase("0")]
    [TestCase("-5")]
    public void PrideFiles_APageSizeBelowOne_IsAUsageError(string pageSize)
    {
        // Left to mzLib, it crossed as a correctness failure (its ArgumentOutOfRangeException).
        UseStub(_ => Json("[]"));

        var ex = Assert.ThrowsAsync<Program.UsageException>(
            async () => await InvokeAsync("pride", "files", "--accession", "PXD012345", "--page-size", pageSize));
        Assert.That(ex!.Message, Does.Contain("--page-size must be greater than zero"));
    }

    [TestCase("files", "--accession", "PXD012345")]
    [TestCase("search", "--keyword", "phospho")]
    public void Pride_APageSizeWithoutAValue_IsAUsageError(string verb, string option, string value)
    {
        UseStub(_ => Json("[]"));

        var ex = Assert.ThrowsAsync<Program.UsageException>(
            async () => await InvokeAsync("pride", verb, option, value, "--page-size"));
        Assert.That(ex!.Message, Is.EqualTo("Option --page-size was given without a value."));
    }

    /// <summary>
    /// A PRIDE manifest whose <c>total_records</c> UNDERSTATES the truth must still arrive whole.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is a regression test for mzLib #1173, and it is here rather than only upstream because
    /// the damage is visible in this verb's own answer: <c>file_count</c> and <c>total_size_bytes</c>
    /// are computed from whatever the pager returned, so a truncated manifest reports a smaller
    /// project rather than an error. A caller downloading "all" the files of a project would have
    /// silently skipped the tail, and nothing in the envelope would have said so.
    /// </para>
    /// <para>
    /// The shape is the one the pager's own comments name as the ambiguous case: page 0 comes back
    /// FULL (two records against a requested page size of two) while the header claims two records
    /// total. Full-and-satisfied is indistinguishable from full-and-one-page-short, so the fix
    /// spends one speculative request before believing the header. Stopping at
    /// <c>Count &gt;= total</c> - what the code did before #1173 - drops <c>c.raw</c> here.
    /// </para>
    /// </remarks>
    [Test]
    public async Task PrideFiles_UnderstatedTotalRecords_StillReturnsTheTail()
    {
        UseStub(request =>
        {
            // The pager asks for page 0, then (because page 0 looked full) probes page 1.
            string uri = request.RequestUri!.ToString();
            string body = uri.Contains("page=1", StringComparison.Ordinal)
                ? $"[{FileJson("c.raw")}]"
                : $"[{FileJson("a.raw")},{FileJson("b.raw")}]";

            HttpResponseMessage response = Json(body);
            // Understated on EVERY page, which is what a server with a stale count actually does.
            response.Headers.Add("total_records", "2");
            return response;
        });

        JsonElement data = await InvokeAsync("pride", "files", "--accession", "PXD012345", "--page-size", "2");

        string[] names = data.GetProperty("files").EnumerateArray()
            .Select(f => f.GetProperty("file_name").GetString()!)
            .ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(names, Is.EqualTo(new[] { "a.raw", "b.raw", "c.raw" }),
                "the record past the understated total was dropped - mzLib #1173 has regressed");
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("total_size_bytes").GetInt64(), Is.EqualTo(3072));
        });
    }

    // ---- pride search --------------------------------------------------------
    //
    // The discovery verb (mzLib #1187). Its payload is NOT a project's metadata: PRIDE serves
    // search from an Elasticsearch projection with every controlled-vocabulary field flattened to
    // a display string, so these fixtures are string collections on purpose.

    private static string SearchHitJson(string accession, string title = "A project") =>
        $$"""
        {
          "accession": "{{accession}}",
          "title": "{{title}}",
          "submissionDate": "2025-11-17",
          "publicationDate": "2026-06-01",
          "keywords": ["phospho", "  ", ""],
          "instruments": ["Q Exactive"],
          "organismsPart": ["liver"],
          "highlights": { "title": ["<em>phospho</em>proteome"] },
          "downloadCount": 0
        }
        """;

    [Test]
    public async Task PrideSearch_ReturnsOneResultPerHit()
    {
        UseStub(_ => Json($"[{SearchHitJson("PXD000111")},{SearchHitJson("PXD000222")}]"));

        JsonElement data = await InvokeAsync("pride", "search", "--keyword", "phospho");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("keyword").GetString(), Is.EqualTo("phospho"));
            Assert.That(data.GetProperty("result_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("results")[0].GetProperty("accession").GetString(),
                Is.EqualTo("PXD000111"));
        });
    }

    [Test]
    public async Task PrideSearch_DatesCrossAsBareCalendarDates()
    {
        // The endpoint sends "2025-11-17" with no time and no offset, which is why mzLib types
        // these DateTime rather than DateTimeOffset. Serialising the DateTime directly would emit
        // "2025-11-17T00:00:00" -- a midnight PRIDE never sent, and exactly the fabrication that
        // type split exists to prevent. Every binding would then inherit the invented time.
        UseStub(_ => Json($"[{SearchHitJson("PXD000111")}]"));

        JsonElement hit = (await InvokeAsync("pride", "search", "--keyword", "x"))
            .GetProperty("results")[0];

        Assert.Multiple(() =>
        {
            Assert.That(hit.GetProperty("submission_date").GetString(), Is.EqualTo("2025-11-17"));
            Assert.That(hit.GetProperty("publication_date").GetString(), Is.EqualTo("2026-06-01"));
        });
    }

    [Test]
    public async Task PrideSearch_ADateTheEndpointOmittedIsNullNotYearOne()
    {
        // The DTO's date members are non-nullable, so an omitted date deserialises to
        // default(DateTime) -- 0001-01-01. Crossing that as a real date hands a caller a year-1
        // value to do arithmetic on, indistinguishable from a genuine (if absurd) date.
        UseStub(_ => Json("""[{ "accession": "PXD000111", "title": "No dates" }]"""));

        JsonElement hit = (await InvokeAsync("pride", "search", "--keyword", "x"))
            .GetProperty("results")[0];

        Assert.Multiple(() =>
        {
            Assert.That(hit.GetProperty("submission_date").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(hit.GetProperty("updated_date").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public async Task PrideSearch_HighlightKeysAreNotRewrittenIntoSnakeCase()
    {
        // The keys are whichever PRIDE fields the query hit -- PRIDE's names, not ours. The
        // serializer's SnakeCaseLower policy renames PROPERTIES; DictionaryKeyPolicy is unset, so
        // these pass through. Pinned because an options change could quietly start rewriting them,
        // and the keys are the only thing that says why a project matched.
        UseStub(_ => Json("""
            [{ "accession": "PXD000111",
               "highlights": { "projectDescription": ["a <em>hit</em>"], "organismsPart": ["liver"] } }]
            """));

        JsonElement highlights = (await InvokeAsync("pride", "search", "--keyword", "x"))
            .GetProperty("results")[0].GetProperty("highlights");

        Assert.Multiple(() =>
        {
            Assert.That(highlights.TryGetProperty("projectDescription", out _), Is.True,
                "the key was rewritten; it is PRIDE's field name, not a property of ours");
            Assert.That(highlights.TryGetProperty("organismsPart", out _), Is.True);
        });
    }

    [Test]
    public async Task PrideSearch_BlankKeywordsAreCarriedRatherThanFiltered()
    {
        // PRIDE ships empty and whitespace-only strings inside keywords on roughly 9% of hits.
        // Dropping them in the bridge would be a repair: every binding would then report a
        // different keyword list than mzLib does, and none of them could tell it had happened.
        UseStub(_ => Json($"[{SearchHitJson("PXD000111")}]"));

        JsonElement keywords = (await InvokeAsync("pride", "search", "--keyword", "x"))
            .GetProperty("results")[0].GetProperty("keywords");

        Assert.That(keywords.GetArrayLength(), Is.EqualTo(3));
    }

    [Test]
    public async Task PrideSearch_FlattenedCvFieldsKeepTheMetadataSpellingOfTheirName()
    {
        // PRIDE spells it "organismsPart" on search and "organismParts" on projects/{accession}.
        // mzLib binds the search spelling onto the metadata name so the two DTOs agree about the
        // same concept; the wire publishes the metadata name in snake_case.
        UseStub(_ => Json($"[{SearchHitJson("PXD000111")}]"));

        JsonElement hit = (await InvokeAsync("pride", "search", "--keyword", "x"))
            .GetProperty("results")[0];

        Assert.That(hit.GetProperty("organism_parts")[0].GetString(), Is.EqualTo("liver"));
    }

    [Test]
    public async Task PrideSearch_NoHits_IsAnEmptyResultNotAnError()
    {
        // A keyword that matches nothing is a real answer, unlike an accession that matches
        // nothing -- there is no typo to protect the caller from.
        UseStub(_ => Json("[]"));

        JsonElement data = await InvokeAsync("pride", "search", "--keyword", "qwertyuiop");

        Assert.That(data.GetProperty("result_count").GetInt32(), Is.Zero);
    }

    [Test]
    public void PrideSearch_AnOverLongKeyword_IsAUsageErrorNotAnOutage()
    {
        // PRIDE answers a very long keyword with HTTP 500, and a 5xx is the one signature the
        // classifier reads as ServiceUnavailable -- so the caller's bug would be reported as EBI
        // being down, and the live-test helpers would SKIP rather than fail. Refused here first.
        string tooLong = new('x', PrideArchiveClient.MaxKeywordLength + 1);

        var usage = Assert.ThrowsAsync<Program.UsageException>(
            async () => await InvokeAsync("pride", "search", "--keyword", tooLong));

        Assert.That(usage!.Message, Does.Contain("at most 1000 characters"));
    }

    [Test]
    public void PrideSearch_WithoutAKeyword_IsAUsageError() =>
        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync("pride", "search"));

    [Test]
    public void PrideSearch_APageSizeOfZero_IsAUsageError() =>
        Assert.ThrowsAsync<Program.UsageException>(
            async () => await InvokeAsync("pride", "search", "--keyword", "x", "--page-size", "0"));

    [Test]
    public async Task PrideSearch_PageSizeAndKeywordReachTheRequestUri()
    {
        string? seen = null;
        Program.PrideClientFactory = () =>
        {
            // AbsoluteUri, not ToString(): Uri.ToString() renders a canonical UNESCAPED form for
            // display, so a correctly escaped query reads as though nothing was escaped at all.
            var handler = new StubHandler(request => { seen = request.RequestUri!.AbsoluteUri; return Json("[]"); });
            return new PrideArchiveClient(new HttpClient(handler) { BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress) });
        };

        await InvokeAsync("pride", "search", "--keyword", "phospho proteome", "--page-size", "7");

        Assert.Multiple(() =>
        {
            Assert.That(seen, Does.Contain("search/projects"));
            Assert.That(seen, Does.Contain("pageSize=7"));
            Assert.That(seen, Does.Contain("keyword=phospho%20proteome"));
        });
    }

    [Test]
    public async Task PrideSearch_AKeywordCannotSmuggleInAQueryParameter()
    {
        // The one that matters. A keyword is caller-supplied text interpolated into a query string,
        // so an unescaped '&' would let "x&pageSize=9999" silently override the paging the caller
        // asked for -- or add a parameter PRIDE honours and nobody intended. mzLib escapes with
        // Uri.EscapeDataString; this pins that it still does, from the far side of the wire.
        string? seen = null;
        Program.PrideClientFactory = () =>
        {
            var handler = new StubHandler(request => { seen = request.RequestUri!.AbsoluteUri; return Json("[]"); });
            return new PrideArchiveClient(new HttpClient(handler) { BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress) });
        };

        await InvokeAsync("pride", "search", "--keyword", "liver&pageSize=9999", "--page-size", "7");

        Assert.Multiple(() =>
        {
            Assert.That(seen, Does.Contain("keyword=liver%26pageSize%3D9999"));
            Assert.That(seen, Does.Contain("pageSize=7"));
            Assert.That(seen, Does.Not.Contain("pageSize=9999"),
                "the keyword's own '&' created a second pageSize parameter");
        });
    }

    // ---- pride ftp-files -----------------------------------------------------
    //
    // The complete listing, walked from the FTP directory tree (mzLib #1121). The stub serves the
    // project metadata (for the publication date that locates the FTP root) and Apache autoindex
    // pages, dispatched by URI — the same fixture shape mzLib's own PrideFtpListingTests uses.

    private static HttpResponseMessage Html(string body) =>
        new(HttpStatusCode.OK) { Content = new StringContent(body, Encoding.UTF8, "text/html") };

    private const string FtpProjectJson =
        """{ "accession": "PXD000001", "publicationDate": "2012-03-13" }""";

    private const string FtpRootUrl = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/";
    private const string FtpGeneratedUrl = FtpRootUrl + "generated/";

    // Three top-level files (one of which the REST manifest omits) plus a subdirectory to descend.
    private const string FtpRootIndexHtml = """
        <html><body><table>
        <tr><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th></tr>
        <tr><td><a href="/pride/data/archive/2012/03/">Parent Directory</a></td><td>&nbsp;</td><td align="right">  - </td></tr>
        <tr><td><a href="README.txt">README.txt</a></td><td align="right">2021-10-20 04:56  </td><td align="right">1.6K</td></tr>
        <tr><td><a href="run1.raw">run1.raw</a></td><td align="right">2021-10-20 04:56  </td><td align="right">210M</td></tr>
        <tr><td><a href="hidden_from_rest.mzML">hidden_from_rest.mzML</a></td><td align="right">2021-10-20 04:56  </td><td align="right">429M</td></tr>
        <tr><td><a href="generated/">generated/</a></td><td align="right">2021-10-20 04:56  </td><td align="right">  - </td></tr>
        </table></body></html>
        """;

    private const string FtpGeneratedIndexHtml = """
        <html><body><table>
        <tr><th><a href="?C=N;O=D">Name</a></th></tr>
        <tr><td><a href="/pride/data/archive/2012/03/PXD000001/">Parent Directory</a></td><td align="right">  - </td></tr>
        <tr><td><a href="summary.mztab">summary.mztab</a></td><td align="right">2021-10-20 04:56  </td><td align="right">844K</td></tr>
        </table></body></html>
        """;

    private static void UseFtpStub() => UseStub(request =>
    {
        string uri = request.RequestUri!.ToString();
        if (uri.EndsWith("/projects/PXD000001", StringComparison.Ordinal)) return Json(FtpProjectJson);
        if (uri == FtpRootUrl) return Html(FtpRootIndexHtml);
        if (uri == FtpGeneratedUrl) return Html(FtpGeneratedIndexHtml);
        return new HttpResponseMessage(HttpStatusCode.NotFound);
    });

    [Test]
    public async Task PrideFtpFiles_WalksTheTreeAndSurfacesFilesTheRestManifestOmits()
    {
        UseFtpStub();

        JsonElement data = await InvokeAsync("pride", "ftp-files", "--accession", "PXD000001");
        string[] paths = data.GetProperty("files").EnumerateArray()
            .Select(f => f.GetProperty("relative_path").GetString()!).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("accession").GetString(), Is.EqualTo("PXD000001"));
            // The whole point of the verb: the complete set, including a subdirectory file and the
            // one a REST manifest would hide. "Parent Directory" and the column-sort links are not files.
            Assert.That(paths, Is.EquivalentTo(new[]
            {
                "README.txt", "run1.raw", "hidden_from_rest.mzML", "generated/summary.mztab",
            }));
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(4));
        });
    }

    [Test]
    public async Task PrideFtpFiles_ProjectsRelativePathLeafNameUrlAndApproximateSize()
    {
        UseFtpStub();

        JsonElement data = await InvokeAsync("pride", "ftp-files", "--accession", "PXD000001");
        JsonElement nested = data.GetProperty("files").EnumerateArray()
            .Single(f => f.GetProperty("relative_path").GetString() == "generated/summary.mztab");

        Assert.Multiple(() =>
        {
            // A nested file keeps its subdirectory in relative_path but exposes a bare leaf file_name,
            // and its url is the HTTPS location a caller can download from.
            Assert.That(nested.GetProperty("file_name").GetString(), Is.EqualTo("summary.mztab"));
            Assert.That(nested.GetProperty("url").GetString(), Is.EqualTo(FtpGeneratedUrl + "summary.mztab"));
            Assert.That(nested.GetProperty("approximate_size_bytes").GetInt64(),
                Is.EqualTo((long)Math.Round(844 * 1024.0)));

            // The reported total is the sum of the (approximate) index sizes, named to say so.
            long expectedTotal = (long)Math.Round(1.6 * 1024)
                + 210L * 1024 * 1024 + 429L * 1024 * 1024 + (long)Math.Round(844 * 1024.0);
            Assert.That(data.GetProperty("approximate_total_size_bytes").GetInt64(), Is.EqualTo(expectedTotal));
        });
    }

    // ---- ftp-files failure classification ------------------------------------
    //
    // The bridge half of the contract the Python layer's ProjectNotFoundError re-map depends on.
    // The FTP walk resolves the project over the REST API first, so an unknown accession is a 404
    // there (not on the autoindex): mzLib turns that into an MzLibException, which must cross as the
    // type string "MzLibException" that list_ftp_files() keys on. Previously covered only in Python.

    [Test]
    public void PrideFtpFiles_UnknownProject_SurfacesMzLibExceptionTheBindingsRemap()
    {
        // Every request 404s, so the REST project lookup returns not-found and mzLib throws
        // MzLibException before any FTP walk begins.
        UseStub(_ => new HttpResponseMessage(HttpStatusCode.NotFound));

        Exception? ex = Assert.CatchAsync(async () =>
            await InvokeAsync("pride", "ftp-files", "--accession", "PXD999999"));

        Assert.That(Program.ClassifyError(ex!), Is.EqualTo("MzLibException"),
            "an unknown accession must cross as MzLibException — the type list_ftp_files re-maps to ProjectNotFoundError");
    }

    [Test]
    public void PrideFtpFiles_ProjectLookupUnavailable_ClassifiesAsServiceUnavailableNotNotFound()
    {
        // A 503 on the project lookup is an outage, not an absent project. mzLib throws
        // HttpRequestException (not MzLibException), so the not-found re-map must NOT swallow it.
        UseStub(_ => new HttpResponseMessage(HttpStatusCode.ServiceUnavailable));

        Exception? ex = Assert.CatchAsync(async () =>
            await InvokeAsync("pride", "ftp-files", "--accession", "PXD999999"));

        Assert.That(Program.ClassifyError(ex!), Is.EqualTo(Program.ServiceUnavailableType));
    }

    [Test]
    public void PrideFtpFiles_NoAccession_IsAUsageError() =>
        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync("pride", "ftp-files"));

    // ---- pride checksums -----------------------------------------------------
    //
    // PRIDE's checksum list (mzLib #1427): text/plain, a header and one tab-separated row per file.
    // The stub answers the checksum route from the test's body and the project route from a JSON
    // project or a 404, the two routes GetFileChecksumsAsync can touch.

    private const string ChecksumHeader = "File-Name\tFile-MD5Checksum\tFile-Size\n";

    // Not in name order, with a size past int.MaxValue and PXD015239's NFS debris row, all of which
    // must cross unchanged.
    private const string ChecksumList = ChecksumHeader
        + "run2.raw\tc567cecbb9866796e690a4958746f7b6\t2285905384\n"
        + "README.txt\t492753137013825025e6729c8798664f\t1645\n"
        + ".nfs80360b91005fae5300002478\te6392618a1e53a280b6385158a91bfeb\t10\n";

    private static HttpResponseMessage Text(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(body, Encoding.UTF8, "text/plain") };

    private static StubHandler UseChecksumStub(string body, bool projectExists)
    {
        var handler = new StubHandler(request =>
        {
            string path = request.RequestUri!.AbsolutePath;
            if (path.Contains("/files/checksum/", StringComparison.Ordinal)) return Text(body);
            if (path.Contains("/projects/", StringComparison.Ordinal) && projectExists)
                return Json("""{ "accession": "PXD000001", "title": "t" }""");
            return new HttpResponseMessage(HttpStatusCode.NotFound);
        });
        Program.PrideClientFactory = () => new PrideArchiveClient(new HttpClient(handler)
        {
            BaseAddress = new Uri(PrideArchiveClient.DefaultBaseAddress),
        });
        return handler;
    }

    [Test]
    public async Task PrideChecksums_ProjectsEveryRowSortedByFileName()
    {
        StubHandler handler = UseChecksumStub(ChecksumList, projectExists: true);

        JsonElement data = await InvokeAsync("pride", "checksums", "--accession", "PXD000001");
        JsonElement[] rows = data.GetProperty("files").EnumerateArray().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("accession").GetString(), Is.EqualTo("PXD000001"));
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(3));
            // mzLib returns a map; the wire sorts by name (ordinal: '.' < 'R' < 'r').
            Assert.That(rows.Select(r => r.GetProperty("file_name").GetString()), Is.EqualTo(new[]
            {
                ".nfs80360b91005fae5300002478", "README.txt", "run2.raw",
            }));
            Assert.That(rows[1].GetProperty("md5").GetString(), Is.EqualTo("492753137013825025e6729c8798664f"));
            Assert.That(rows[1].GetProperty("size_bytes").GetInt64(), Is.EqualTo(1645));
            // A raw file past 2 GiB: the size is a 64-bit count, not an int.
            Assert.That(rows[2].GetProperty("size_bytes").GetInt64(), Is.EqualTo(2_285_905_384L));
            Assert.That(rows[0].EnumerateObject().Select(p => p.Name),
                Is.EqualTo(new[] { "file_name", "md5", "size_bytes" }));
            // No total: mzLib has none over this list, so the wire sends none.
            Assert.That(data.TryGetProperty("total_size_bytes", out _), Is.False);
            // One request: a non-empty list needs no existence check.
            Assert.That(handler.RequestedUris, Is.EqualTo(new[]
            {
                "https://www.ebi.ac.uk/pride/ws/archive/v3/files/checksum/PXD000001",
            }));
        });
    }

    [Test]
    public async Task PrideChecksums_TheMd5CrossesAsPrideWroteIt()
    {
        UseChecksumStub(ChecksumHeader + "run1.raw\t0123456789ABCDEF0123456789ABCDEF\t10\n", projectExists: true);

        JsonElement data = await InvokeAsync("pride", "checksums", "--accession", "PXD000001");

        // mzLib keeps PRIDE's text, upper case included; the bridge does not rewrite it either.
        Assert.That(data.GetProperty("files")[0].GetProperty("md5").GetString(),
            Is.EqualTo("0123456789ABCDEF0123456789ABCDEF"));
    }

    [Test]
    public async Task PrideChecksums_AnExistingProjectWithNoChecksums_IsAnEmptyListNotAnError()
    {
        StubHandler handler = UseChecksumStub("", projectExists: true);

        JsonElement data = await InvokeAsync("pride", "checksums", "--accession", "PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("files").GetArrayLength(), Is.Zero);
            // The empty body made mzLib ask whether the project exists.
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(2));
            Assert.That(handler.RequestedUris[1], Does.EndWith("/projects/PXD000001"));
        });
    }

    [Test]
    public void PrideChecksums_UnknownAccession_SurfacesMzLibExceptionTheBindingsRemap()
    {
        // PRIDE answers an accession it does not know with 200 and an empty body, then 404 on the project.
        UseChecksumStub("", projectExists: false);

        Exception? ex = Assert.CatchAsync(async () =>
            await InvokeAsync("pride", "checksums", "--accession", "PXD999999999"));

        Assert.Multiple(() =>
        {
            Assert.That(Program.ClassifyError(ex!), Is.EqualTo("MzLibException"));
            // The wording pyMzLib matches to raise ProjectNotFoundError.
            Assert.That(ex!.Message, Does.Contain("has no project with accession 'PXD999999999'"));
        });
    }

    [TestCase("File-Name\tMD5\tFile-Size\nrun1.raw\t0123456789abcdef0123456789abcdef\t10\n", "starts with")]
    [TestCase(ChecksumHeader + "run1.raw\t0123456789abcdef\t10\n", "not 32 hexadecimal characters")]
    [TestCase(ChecksumHeader + "run1.raw\t0123456789abcdef0123456789abcdef\t-10\n", "not a whole number of bytes")]
    [TestCase(ChecksumHeader + "a.raw\t0123456789abcdef0123456789abcdef\t10\na.raw\t0123456789abcdef0123456789abcdef\t10\n", "twice")]
    public void PrideChecksums_AMalformedList_IsACorrectnessFailureAndIsNotRetried(string body, string expected)
    {
        StubHandler handler = UseChecksumStub(body, projectExists: true);

        Exception? ex = Assert.CatchAsync(async () =>
            await InvokeAsync("pride", "checksums", "--accession", "PXD000001"));

        Assert.Multiple(() =>
        {
            Assert.That(Program.ClassifyError(ex!), Is.EqualTo("MzLibException"));
            Assert.That(ex!.Message, Does.Contain(expected));
            // Not the not-found wording, so no binding reports a broken list as a typo.
            Assert.That(ex.Message, Does.Not.Contain("has no project with accession"));
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(1), "a broken list is not an outage");
        });
    }

    [Test]
    public void PrideChecksums_A404OnTheList_IsACorrectnessFailure()
    {
        // 404 is not transient, so mzLib does not retry it and the test does not wait out a backoff.
        UseStub(_ => new HttpResponseMessage(HttpStatusCode.NotFound));

        Exception? ex = Assert.CatchAsync(async () =>
            await InvokeAsync("pride", "checksums", "--accession", "PXD000001"));

        Assert.That(Program.ClassifyError(ex!), Is.EqualTo(nameof(HttpRequestException)));
    }

    // The options after the verb, '|'-separated (an array argument would be spread by TestCase).
    [TestCase("", TestName = "PrideChecksums_NoAccession_IsAUsageError")]
    [TestCase("--accession", TestName = "PrideChecksums_AnAccessionWithoutAValue_IsAUsageError")]
    [TestCase("--accession|   ", TestName = "PrideChecksums_ABlankAccession_IsAUsageError")]
    public void PrideChecksums_AMissingAccession_IsAUsageError(string options)
    {
        StubHandler handler = UseChecksumStub(ChecksumList, projectExists: true);
        string[] args = new[] { "pride", "checksums" }
            .Concat(options.Split('|', StringSplitOptions.RemoveEmptyEntries)).ToArray();

        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(args));
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    // ---- failure classification ---------------------------------------------
    //
    // This is the behavior the whole external-service convention rests on: a caller must be able to
    // tell "EBI is down" from "this is broken" without reading prose.

    [TestCase(500, Description = "server error")]
    [TestCase(503, Description = "service unavailable")]
    [TestCase(429, Description = "rate limited")]
    [TestCase(408, Description = "request timeout")]
    public void AvailabilityFailures_AreLabelledServiceUnavailable(int status)
    {
        var exception = new HttpRequestException($"PRIDE Archive request failed with status {status} Whatever for 'projects/x/files'.");

        Assert.That(Program.ClassifyError(exception), Is.EqualTo(Program.ServiceUnavailableType));
    }

    [TestCase(404, Description = "wrong URL is our problem, not EBI's")]
    [TestCase(400, Description = "malformed request is our problem")]
    public void ClientErrors_KeepTheirOwnTypeSoTheyAreNotExcusedAsAnOutage(int status)
    {
        var exception = new HttpRequestException($"PRIDE Archive request failed with status {status} Whatever for 'projects/x/files'.");

        Assert.That(Program.ClassifyError(exception), Is.EqualTo(nameof(HttpRequestException)));
    }

    [Test]
    public void TimeoutsCountAsUnavailable()
    {
        Assert.That(Program.ClassifyError(new TaskCanceledException()), Is.EqualTo(Program.ServiceUnavailableType));
    }

    [Test]
    public void UnreachableHostCountsAsUnavailable()
    {
        Assert.That(Program.ClassifyError(new System.Net.Sockets.SocketException(10061)),
            Is.EqualTo(Program.ServiceUnavailableType));
    }

    [Test]
    public void ATransportFailureWithNoStatusCodeCountsAsUnavailable()
    {
        // The regression this guards was live and inverted the whole convention. Connection
        // refused, DNS failure and TLS failure are the commonest outage modes there are, and all
        // three produce an HttpRequestException with no status because no response ever arrived.
        // HttpClient wraps the SocketException, so the socket arm never sees it, and these used to
        // be classified as correctness failures - turning every real outage into a red build.
        // Constructed the way HttpClient actually constructs them: the real cause is the inner
        // exception, which is also what distinguishes these from mzLib's hand-thrown ones.
        Assert.Multiple(() =>
        {
            Assert.That(Program.ClassifyError(new HttpRequestException(
                    "Connection refused (ftp.pride.ebi.ac.uk:443)",
                    new System.Net.Sockets.SocketException(10061))),
                Is.EqualTo(Program.ServiceUnavailableType), "connection refused");
            Assert.That(Program.ClassifyError(new HttpRequestException(
                    "No such host is known.", new System.Net.Sockets.SocketException(11001))),
                Is.EqualTo(Program.ServiceUnavailableType), "DNS failure");
            Assert.That(Program.ClassifyError(new HttpRequestException(
                    "The SSL connection could not be established.",
                    new System.Security.Authentication.AuthenticationException("bad cert"))),
                Is.EqualTo(Program.ServiceUnavailableType), "TLS failure");
        });
    }

    [Test]
    public void ATruncatedResponseBodyCountsAsUnavailable()
    {
        // The live regression, 5 August 2026: the PRIDE download canary went red on
        // "Received an unexpected EOF or 0 bytes from the transport stream." EBI hung up partway
        // through the body. By then HttpClient had already returned the headers, so this never
        // reached the HttpRequestException arm - it surfaced on the content stream as a bare
        // IOException and crossed the wire as the type "IOException", which the Python and R
        // suites read as a contract break rather than an outage.
        Assert.Multiple(() =>
        {
            Assert.That(Program.ClassifyError(new IOException(
                    "Received an unexpected EOF or 0 bytes from the transport stream.")),
                Is.EqualTo(Program.ServiceUnavailableType), "TLS stream truncated mid-body");
            Assert.That(Program.ClassifyError(new HttpIOException(
                    HttpRequestError.ResponseEnded, "The response ended prematurely.")),
                Is.EqualTo(Program.ServiceUnavailableType), "response ended prematurely");
            Assert.That(Program.ClassifyError(new IOException(
                    "Unable to read data from the transport connection.",
                    new System.Net.Sockets.SocketException(10054))),
                Is.EqualTo(Program.ServiceUnavailableType), "connection reset mid-read");
        });
    }

    [Test]
    public void AFullDiskIsNotExcusedAsAnOutage()
    {
        // The counterweight to the test above, and the reason that arm is narrow rather than a
        // blanket `IOException => ServiceUnavailable`. The bridge writes downloads to disk, so the
        // same code path raises IOException for reasons that are entirely ours. Excusing those
        // would skip the canary on a broken runner and report nothing.
        Assert.Multiple(() =>
        {
            Assert.That(Program.ClassifyError(new IOException("There is not enough space on the disk.")),
                Is.EqualTo(nameof(IOException)), "disk full");
            Assert.That(Program.ClassifyError(new FileNotFoundException("Could not find file 'x.raw'.")),
                Is.EqualTo(nameof(FileNotFoundException)), "missing file");
            Assert.That(Program.ClassifyError(new DirectoryNotFoundException("Could not find a part of the path.")),
                Is.EqualTo(nameof(DirectoryNotFoundException)), "missing directory");
        });
    }

    [Test]
    public void TheStatusPatternCannotBeSteeredByCallerSuppliedText()
    {
        // An unanchored "status NNN" match could be driven by any caller string that reaches the
        // message - an accession of "x status 503 x" would classify itself as an outage and be
        // skipped by every test suite.
        var exception = new HttpRequestException(
            "PRIDE Archive paging exceeded 3 pages for accession 'x status 503 x'.");

        Assert.That(Program.ClassifyError(exception), Is.EqualTo(nameof(HttpRequestException)));
    }

    [Test]
    public void AProgrammingErrorIsNeverExcusedAsAnOutage()
    {
        // The failure mode this guards against is the worst one: a real bug quietly reported as
        // "the service is down", skipped by every test suite, and never seen again.
        Assert.That(Program.ClassifyError(new NullReferenceException()), Is.EqualTo(nameof(NullReferenceException)));
        Assert.That(Program.ClassifyError(new JsonException("bad json")), Is.EqualTo(nameof(JsonException)));
    }

    [Test]
    public void StatusCodeIsReadFromThePropertyWhenItIsSet()
    {
        // mzLib composes its message by hand and leaves StatusCode null, but other callers may not.
        var exception = new HttpRequestException("no status in this text", null, HttpStatusCode.ServiceUnavailable);

        Assert.That(Program.ClassifyError(exception), Is.EqualTo(Program.ServiceUnavailableType));
    }

    // ---- the filter must never fail open -------------------------------------
    //
    // A filter that was asked for but degenerates to nothing used to leave the selection null,
    // which downloaded the entire project - hundreds of megabytes where an error was wanted.

    [TestCase("   ", Description = "whitespace category")]
    [TestCase("", Description = "empty category")]
    public void ABlankCategoryIsRejectedRatherThanSelectingEverything(string category)
    {
        UseStub(_ => Json($"[{FileJson("a.raw")}]"));

        Assert.That(async () => await InvokeAsync(
                "pride", "download", "--accession", "PXD012345", "--dest", "out", "--category", category),
            Throws.InstanceOf<Program.UsageException>());
    }

    [Test]
    public void AnExtensionListThatNamesNothingIsRejected()
    {
        UseStub(_ => Json($"[{FileJson("a.raw")}]"));

        Assert.That(async () => await InvokeAsync(
                "pride", "download", "--accession", "PXD012345", "--dest", "out", "--ext", ",,"),
            Throws.InstanceOf<Program.UsageException>());
    }

    [Test]
    public void AnOptionWrittenWithoutAValueStillCountsAsProvided()
    {
        // "--category --ext .raw" parses the category as a valueless flag. If that reads as
        // "category was not requested", the caller's intent is discarded and the download widens
        // to the whole project.
        var arguments = new Program.Arguments(["pride", "download", "--category", "--dest", "out"]);

        Assert.Multiple(() =>
        {
            Assert.That(arguments.WasProvided("category"), Is.True);
            Assert.That(arguments.Optional("category"), Is.Null);
        });
    }

    [Test]
    public async Task AnAbsentFilterSelectsEveryFile()
    {
        // The counterpart: omitting the options entirely is a legitimate "give me all of it".
        // Asserting only Throws.Nothing would keep passing if a regression narrowed the selection
        // to none, which is the exact failure this pair of tests exists to detect.
        UseStub(_ => Json($"[{FileJson("a.raw")},{FileJson("b.raw")}]"));

        JsonElement data = await InvokeAsync("pride", "download", "--accession", "PXD012345", "--dest", "out");

        Assert.That(data.GetProperty("downloaded_count").GetInt32(), Is.EqualTo(2));
    }

    // ---- explicit selection over stdin ---------------------------------------
    //
    // The headline capability of this change, and previously untested on the C# side.

    [Test]
    public async Task NamesFromStdinSelectExactlyThoseFiles()
    {
        UseStub(_ => Json($"[{FileJson("a.raw")},{FileJson("b.raw")},{FileJson("c.raw")}]"));
        Console.SetIn(new StringReader("a.raw\nc.raw\n"));

        JsonElement data = await InvokeAsync(
            "pride", "download", "--accession", "PXD012345", "--dest", "out", "--names-from-stdin");

        Assert.That(data.GetProperty("downloaded_count").GetInt32(), Is.EqualTo(2));
    }

    [Test]
    public async Task AUtf8BomOnTheFirstSelectedNameIsStripped()
    {
        // A BOM (Windows PowerShell 5.1 pipes one) must not turn "a.raw" into "﻿a.raw", which
        // matches nothing and fails as "not in project" after every other file has downloaded.
        UseStub(_ => Json($"[{FileJson("a.raw")},{FileJson("b.raw")}]"));
        Console.SetIn(new StringReader("﻿a.raw\nb.raw\n"));

        JsonElement data = await InvokeAsync(
            "pride", "download", "--accession", "PXD012345", "--dest", "out", "--names-from-stdin");

        Assert.That(data.GetProperty("downloaded_count").GetInt32(), Is.EqualTo(2));
    }

    [Test]
    public void ARequestedNameThatIsNotInTheProjectIsReported()
    {
        // Silently downloading fewer files than asked for is the failure mode this whole change
        // set exists to remove: a typo would otherwise produce a short download and a success.
        UseStub(_ => Json($"[{FileJson("a.raw")}]"));
        Console.SetIn(new StringReader("a.raw\ntypo.raw\n"));

        var ex = Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "pride", "download", "--accession", "PXD012345", "--dest", "out", "--names-from-stdin"));

        Assert.That(ex!.Message, Does.Contain("typo.raw"));
    }

    [Test]
    public void AnExplicitSelectionAndAFilterAreContradictory()
    {
        UseStub(_ => Json($"[{FileJson("a.raw")}]"));
        Console.SetIn(new StringReader("a.raw\n"));

        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "pride", "download", "--accession", "PXD012345", "--dest", "out",
            "--names-from-stdin", "--category", "RAW"));
    }

    [Test]
    public void AnEmptySelectionOnStdinIsRejected()
    {
        UseStub(_ => Json($"[{FileJson("a.raw")}]"));
        Console.SetIn(new StringReader("   \n\n"));

        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "pride", "download", "--accession", "PXD012345", "--dest", "out", "--names-from-stdin"));
    }

    [Test]
    public void AnHttpFailureWithNoStatusAndNoInnerCauseIsOurProblem()
    {
        // The other side of the inner-exception rule. mzLib composes exceptions by hand to signal
        // conditions that are ours — the paging guard, for one — and those must not be excused as
        // outages, or they would be skipped by every suite and never seen again.
        Assert.That(Program.ClassifyError(new HttpRequestException("something we did wrong")),
            Is.EqualTo(nameof(HttpRequestException)));
    }

    // ---- helper --------------------------------------------------------------

    /// <summary>Runs a verb through the real dispatcher and returns its data payload as JSON.</summary>
    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        object data = await Program.DispatchAsync(args);
        return JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(data, Program.JsonOptions));
    }
}
