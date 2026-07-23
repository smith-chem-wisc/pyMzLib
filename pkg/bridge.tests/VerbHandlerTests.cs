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
