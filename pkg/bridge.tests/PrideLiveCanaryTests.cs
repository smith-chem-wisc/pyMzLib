using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Live canary against the real PRIDE Archive API.
/// </summary>
/// <remarks>
/// Everything else in this suite runs against stubs, which means it would keep passing forever
/// even if EBI changed the shape of its responses tomorrow. This is the test that would notice.
/// <para>
/// It carries <c>[Category("ExternalService")]</c> so CI runs it in a dedicated, non-blocking job
/// rather than the required unit-test run — matching mzLib's convention — and every call routes
/// through <see cref="ExternalServiceTestHelper.RunAsync"/> so that a PRIDE outage <b>skips</b>
/// with an explanatory message while a genuine contract break still <b>fails</b>. The point is
/// that a red build should always mean "we broke something", never "EBI is having a bad morning".
/// </para>
/// </remarks>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideLiveCanaryTests
{
    /// <summary>A small, long-stable public project. Changing this changes what the canary proves.</summary>
    private const string CanaryAccession = "PXD000001";

    [Test]
    public Task TheApiStillAnswersAndTheManifestStillParses() =>
        ExternalServiceTestHelper.RunAsync("PRIDE Archive", async () =>
        {
            JsonElement data = await InvokeAsync("pride", "files", "--accession", CanaryAccession);

            Assert.Multiple(() =>
            {
                Assert.That(data.GetProperty("file_count").GetInt32(), Is.GreaterThan(0),
                    "PRIDE answered but reported no files for a project known to have them — the response shape has probably changed.");
                Assert.That(data.GetProperty("total_size_bytes").GetInt64(), Is.GreaterThan(0));
            });
        });

    [Test]
    public Task TheFieldsWeDependOnAreStillPopulated() =>
        ExternalServiceTestHelper.RunAsync("PRIDE Archive", async () =>
        {
            JsonElement data = await InvokeAsync("pride", "files", "--accession", CanaryAccession);
            JsonElement file = data.GetProperty("files")[0];

            // Each of these is something the Python layer reads. If PRIDE stops sending one, we
            // want to hear about it here rather than from a user.
            Assert.Multiple(() =>
            {
                Assert.That(file.GetProperty("file_name").GetString(), Is.Not.Empty);
                Assert.That(file.GetProperty("category").GetString(), Is.Not.Empty);
                Assert.That(file.GetProperty("file_size_bytes").GetInt64(), Is.GreaterThan(0));
                Assert.That(file.GetProperty("submission_date").GetString(), Is.Not.Empty);
            });
        });

    [Test]
    public Task AtLeastOneFileIsStillReachableOverHttps() =>
        ExternalServiceTestHelper.RunAsync("PRIDE Archive", async () =>
        {
            JsonElement data = await InvokeAsync("pride", "files", "--accession", CanaryAccession);

            bool anyHttps = data.GetProperty("files").EnumerateArray()
                .Any(f => f.GetProperty("https_url").ValueKind is not JsonValueKind.Null);

            // The FTP-to-HTTPS upgrade is an assumption about how EBI publishes locations, not a
            // guarantee. If it ever stops holding, downloads break and this is the early warning.
            Assert.That(anyHttps, Is.True,
                "No file exposed an HTTPS location — the FTP-to-HTTPS upgrade assumption may no longer hold.");
        });

    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        object data = await Program.DispatchAsync(args);
        return JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(data, Program.JsonOptions));
    }
}
