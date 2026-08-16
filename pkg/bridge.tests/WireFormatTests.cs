using System.Text.Json;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the JSON contract itself.
/// </summary>
/// <remarks>
/// These matter more than they look. The wire format is a published interface — pyMzLib parses
/// it, and by design any other language could too — so a change to a key name or a null is a
/// breaking change to consumers that no compiler will catch. Pinning the shape here means such
/// a change has to be deliberate.
/// </remarks>
[TestFixture]
public class WireFormatTests
{
    private static PrideArchiveFile FileWith(params CvParam[] locations) => new()
    {
        FileName = "run1.raw",
        FileSizeBytes = 1234,
        Checksum = "abc",
        FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "Raw", "RAW"),
        PublicFileLocations = [.. locations],
        SubmissionDate = new DateTimeOffset(2012, 3, 13, 0, 0, 0, TimeSpan.Zero),
    };

    private static CvParam Ftp(string url) => new("PRIDE", "PRIDE:0000469", "FTP Protocol", url);

    private static CvParam Aspera(string url) => new("PRIDE", "PRIDE:0000468", "Aspera Protocol", url);

    private static JsonElement Serialize(object value) =>
        JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(value, Program.JsonOptions));

    [Test]
    public void KeysAreSnakeCaseSoTheyReadAsIdiomaticPython()
    {
        JsonElement json = Serialize(Program.ToWireFile(FileWith(Ftp("ftp://x/run1.raw"))));

        Assert.Multiple(() =>
        {
            Assert.That(json.TryGetProperty("file_name", out _), Is.True);
            Assert.That(json.TryGetProperty("file_size_bytes", out _), Is.True);
            Assert.That(json.TryGetProperty("category_accession", out _), Is.True);
            Assert.That(json.TryGetProperty("fileName", out _), Is.False, "camelCase would break callers");
        });
    }

    [Test]
    public void CategoryIsFlattenedToItsValueAndAccession()
    {
        // A caller should not have to understand controlled-vocabulary structure to read a
        // manifest, but the accession stays available because names are not stable identifiers.
        JsonElement json = Serialize(Program.ToWireFile(FileWith(Ftp("ftp://x/run1.raw"))));

        Assert.Multiple(() =>
        {
            Assert.That(json.GetProperty("category").GetString(), Is.EqualTo("RAW"));
            Assert.That(json.GetProperty("category_accession").GetString(), Is.EqualTo("PRIDE:0000404"));
        });
    }

    [Test]
    public void FtpLocationsAreUpgradedToHttps()
    {
        JsonElement json = Serialize(Program.ToWireFile(FileWith(Ftp("ftp://ftp.pride.ebi.ac.uk/x/run1.raw"))));

        Assert.That(json.GetProperty("https_url").GetString(),
            Is.EqualTo("https://ftp.pride.ebi.ac.uk/x/run1.raw"));
    }

    [Test]
    public void AsperaOnlyFilesReportNullRatherThanAnUnusableUrl()
    {
        // pyMzLib's `downloadable` property is derived from this being null. Returning the
        // Aspera URL here would make a file look fetchable and fail later, at download time.
        JsonElement json = Serialize(Program.ToWireFile(FileWith(Aspera("prd_ascp@fasp.ebi.ac.uk:x/run1.raw"))));

        Assert.That(json.GetProperty("https_url").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void EveryPublishedLocationIsPreserved()
    {
        JsonElement json = Serialize(Program.ToWireFile(
            FileWith(Ftp("ftp://x/run1.raw"), Aspera("prd_ascp@fasp.ebi.ac.uk:x/run1.raw"))));

        Assert.That(json.GetProperty("locations").GetArrayLength(), Is.EqualTo(2));
    }

    [Test]
    public void DatesAreIso8601WithAnOffset()
    {
        JsonElement json = Serialize(Program.ToWireFile(FileWith(Ftp("ftp://x/run1.raw"))));

        // Python parses these with datetime.fromisoformat; a naive timestamp would silently
        // produce timezone-unaware datetimes on the other side.
        Assert.That(json.GetProperty("submission_date").GetString(), Does.StartWith("2012-03-13T00:00:00"));
        Assert.That(json.GetProperty("submission_date").GetString(), Does.Contain("+00:00"));
    }

    [Test]
    public async Task VersionReportsWhichMzLibItWasBuiltAgainst()
    {
        // A native C# consumer never needs this: FlashLFQ references mzLib and can read its version
        // directly. A binding holds a prebuilt binary fetched from a release and cannot, so without
        // this field "which mzLib produced this result?" has no runtime answer in any of the three
        // languages — only a pin file in a repository the user may not have.
        JsonElement json = Serialize(new Program.Envelope
        {
            Ok = true,
            Data = await Program.DispatchAsync(["version"]),
        });

        JsonElement data = json.GetProperty("data");

        Assert.Multiple(() =>
        {
            // The compatibility contract stays `protocol`. Adding a field must not move it, which is
            // exactly why the bindings could take this addition without a protocol bump.
            Assert.That(data.GetProperty("protocol").GetInt32(), Is.EqualTo(1));

            // Read out of the linked assembly's informational version, so it describes the mzLib
            // that is genuinely in this binary rather than one passed in alongside it.
            string? mzlib = data.GetProperty("mzlib").GetString();
            Assert.That(mzlib, Is.Not.Null.And.Not.Empty,
                "the bridge should report the mzLib it was built against; a null here means the "
                + "build recorded no source commit, which makes a release untraceable to its source");
            Assert.That(mzlib, Does.Contain("+"),
                "the value is `<version>+<commit>`; without the commit it names no particular mzLib");
        });
    }

    [Test]
    public void SuccessEnvelopeCarriesDataAndANullError()
    {
        JsonElement json = Serialize(new Program.Envelope { Ok = true, Data = new { a = 1 } });

        Assert.Multiple(() =>
        {
            Assert.That(json.GetProperty("ok").GetBoolean(), Is.True);
            Assert.That(json.GetProperty("data").GetProperty("a").GetInt32(), Is.EqualTo(1));
            Assert.That(json.GetProperty("error").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void FailureEnvelopeNamesTheErrorType()
    {
        // The type is the machine-readable half: it is what lets a caller distinguish a network
        // failure from a bad accession without matching on prose.
        JsonElement json = Serialize(new Program.Envelope
        {
            Ok = false,
            Error = new Program.ErrorInfo { Type = "HttpRequestException", Message = "status 503" },
        });

        Assert.Multiple(() =>
        {
            Assert.That(json.GetProperty("ok").GetBoolean(), Is.False);
            Assert.That(json.GetProperty("data").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(json.GetProperty("error").GetProperty("type").GetString(),
                Is.EqualTo("HttpRequestException"));
        });
    }

    [Test]
    public void EnvelopeAlwaysHasAllThreeKeys()
    {
        // Nulls are written rather than omitted so a consumer in a statically typed language can
        // deserialize into one fixed shape regardless of outcome.
        JsonElement json = Serialize(new Program.Envelope { Ok = true, Data = null });

        Assert.Multiple(() =>
        {
            Assert.That(json.TryGetProperty("ok", out _), Is.True);
            Assert.That(json.TryGetProperty("data", out _), Is.True);
            Assert.That(json.TryGetProperty("error", out _), Is.True);
        });
    }
}
