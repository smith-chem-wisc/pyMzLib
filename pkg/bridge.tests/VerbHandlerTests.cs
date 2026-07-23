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

    // ---- helper --------------------------------------------------------------

    /// <summary>Runs a verb through the real dispatcher and returns its data payload as JSON.</summary>
    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        object data = await Program.DispatchAsync(args);
        return JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(data, Program.JsonOptions));
    }
}
