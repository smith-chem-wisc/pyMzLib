using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the <c>readers</c> verbs — format dispatch and the view reporting that keeps a caller
/// from assuming a uniform shape that does not exist.
/// </summary>
/// <remarks>
/// <para>
/// The parsing is mzLib's and mzLib tests it. What lives only in the bridge is the projection:
/// which views a reader is reported to support, and the translation of mzLib's two exception types
/// into a usage error a caller can act on. Those need no real data — mzLib dispatches on the
/// filename for every type exercised here, so empty files with the right names suffice, and the
/// suite stays fast and offline.
/// </para>
/// <para>
/// The exception is <c>.mztab</c>, which mzLib disambiguates by reading the first five lines. That
/// is used deliberately below to cover the content-sniffing branch.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReadingTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-readers-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void Cleanup()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    /// <summary>Writes a file with the given name and returns its full path.</summary>
    private string Touch(string name, string contents = "")
    {
        string path = Path.Combine(_tempDirectory, name);
        File.WriteAllText(path, contents);
        return path;
    }

    /// <summary>Runs a verb and returns its data as JSON, the shape a caller actually receives.</summary>
    private static JsonElement Run(params string[] args)
    {
        var arguments = new Program.Arguments(args);
        object data = arguments.Verb switch
        {
            "readers identify" => Reading.Identify(arguments),
            "readers formats" => Reading.Formats(arguments),
            _ => throw new ArgumentException($"Not a readers verb: {arguments.Verb}"),
        };

        return JsonSerializer.SerializeToElement(data, Program.JsonOptions);
    }

    private static string[] ViewsOf(JsonElement element) =>
        element.GetProperty("views").EnumerateArray().Select(v => v.GetString()!).ToArray();

    // ---- identify: the view families -------------------------------------------------------

    [Test]
    public void Identify_MetaMorpheusPsmtsv_ReportsQuantifiableView()
    {
        JsonElement result = Run("readers", "identify", "--path", Touch("AllPSMs.psmtsv"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("psmtsv"));
        Assert.That(result.GetProperty("reader").GetString(), Is.EqualTo("PsmFromTsvFile"));
        Assert.That(ViewsOf(result), Does.Contain("quantifiable"));
    }

    [Test]
    public void Identify_MsFraggerPsm_ReportsQuantifiableView()
    {
        JsonElement result = Run("readers", "identify", "--path", Touch("psm.tsv"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("MsFraggerPsm"));
        Assert.That(ViewsOf(result), Does.Contain("quantifiable"));
    }

    [Test]
    public void Identify_ToppicPrsm_ReportsNoUniformView()
    {
        // The honesty this verb exists for. TopPIC is read perfectly well by mzLib but implements
        // no cross-format interface, so there is nothing to project it onto. A caller must be able
        // to discover that here rather than as a cast failure inside a later call — and the
        // Readers plan's original acceptance criterion assumed the opposite.
        JsonElement result = Run("readers", "identify", "--path", Touch("run_prsm.tsv"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("ToppicPrsm"));
        Assert.That(ViewsOf(result), Is.Empty, "TopPIC implements no uniform view interface");
    }

    [Test]
    public void Identify_Ms1Feature_ReportsMs1FeaturesViewNotQuantifiable()
    {
        JsonElement result = Run("readers", "identify", "--path", Touch("run_ms1.feature"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("Ms1Feature"));
        Assert.That(ViewsOf(result), Is.EqualTo(new[] { "ms1_features" }));
    }

    [Test]
    public void Identify_SpectraFile_ReportsSpectraView()
    {
        JsonElement result = Run("readers", "identify", "--path", Touch("run.mzML"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("MzML"));
        Assert.That(ViewsOf(result), Is.EqualTo(new[] { "spectra" }),
            "a spectra file is not a result file and must not claim a record view");
    }

    [Test]
    public void Identify_CasanovoMzTab_ReportsSpectralMatchView_FromFileContents()
    {
        // .mztab is one of the three types mzLib disambiguates by reading the file rather than the
        // name, so this also covers the content-sniffing branch.
        JsonElement result = Run("readers", "identify", "--path",
            Touch("denovo.mztab", "MTD\tmzTab-version\t1.0.0\nMTD\tsoftware[1]\t[MS, MS:1003281, casanovo, 5.0.0]\n"));

        Assert.That(result.GetProperty("file_type").GetString(), Is.EqualTo("CasanovoMzTab"));
        Assert.That(ViewsOf(result), Is.EqualTo(new[] { "spectral_match" }));
    }

    // ---- identify: the wire shape and its deliberate omission -------------------------------

    [Test]
    public void Identify_ReportsAbsolutePath_AndNoSoftwareField()
    {
        string path = Touch("AllPSMs.psmtsv");
        JsonElement result = Run("readers", "identify", "--path", path);

        Assert.That(result.GetProperty("path").GetString(), Is.EqualTo(Path.GetFullPath(path)));
        Assert.That(result.GetProperty("extension").GetString(), Is.EqualTo(".psmtsv"));

        // Pinned deliberately: mzLib's IResultFile.Software is Unspecified for everything its own
        // factory returns, so a software field here could only be reconstructed — i.e. invented.
        // If someone adds one, this test should make them justify where the value came from.
        Assert.That(result.TryGetProperty("software", out _), Is.False,
            "software cannot be answered reliably by mzLib; see Reading.Identify remarks");
    }

    // ---- identify: failures -----------------------------------------------------------------

    [Test]
    public void Identify_MissingFile_IsAUsageErrorNamingThePath()
    {
        string missing = Path.Combine(_tempDirectory, "absent.psmtsv");

        var exception = Assert.Throws<Program.UsageException>(
            () => Run("readers", "identify", "--path", missing));

        // mzLib throws a bare FileNotFoundException carrying neither message nor file name, which
        // would reach the caller as an empty error.
        Assert.That(exception!.Message, Does.Contain("absent.psmtsv"));
    }

    [Test]
    public void Identify_UnsupportedExtension_IsAUsageErrorPointingAtTheFormatsVerb()
    {
        var exception = Assert.Throws<Program.UsageException>(
            () => Run("readers", "identify", "--path", Touch("notes.docx")));

        Assert.That(exception!.Message, Does.Contain("readers formats"),
            "an unsupported file is a bad argument, and the caller needs to know what IS supported");
    }

    [Test]
    public void Identify_MissingPathOption_IsAUsageError() =>
        Assert.Throws<Program.UsageException>(() => Run("readers", "identify"));

    // ---- formats -----------------------------------------------------------------------------

    [Test]
    public void Formats_EnumeratesEverySupportedTypeWithAnExtensionAndReader()
    {
        JsonElement result = Run("readers", "formats");
        JsonElement formats = result.GetProperty("formats");

        Assert.That(result.GetProperty("format_count").GetInt32(), Is.EqualTo(formats.GetArrayLength()));
        Assert.That(formats.GetArrayLength(), Is.EqualTo(Enum.GetValues<Readers.SupportedFileType>().Length),
            "the table is enumerated from mzLib, so it cannot drift from what mzLib dispatches");

        foreach (JsonElement format in formats.EnumerateArray())
        {
            Assert.That(format.GetProperty("file_type").GetString(), Is.Not.Empty);
            Assert.That(format.GetProperty("extension").ValueKind, Is.Not.EqualTo(JsonValueKind.Null),
                $"every supported type should map to an extension: {format.GetProperty("file_type")}");
            Assert.That(format.GetProperty("reader").ValueKind, Is.Not.EqualTo(JsonValueKind.Null),
                $"every supported type should map to a reader: {format.GetProperty("file_type")}");
        }
    }

    [Test]
    public void Formats_ExactlyThreeTypesOfferTheQuantifiableView()
    {
        // The headline fact about this tranche, pinned so a change in mzLib is DETECTED rather than
        // silently widening or narrowing what pyMzLib claims. mzLib reads 29 formats; only these
        // three implement IQuantifiableResultFile and can therefore feed flashlfq.quantify().
        // If mzLib adds one, this test fails and the docs get updated — which is the point.
        string[] quantifiable = Run("readers", "formats").GetProperty("formats").EnumerateArray()
            .Where(f => ViewsOf(f).Contains("quantifiable"))
            .Select(f => f.GetProperty("file_type").GetString()!)
            .ToArray();

        Assert.That(quantifiable, Is.EquivalentTo(new[] { "psmtsv", "osmtsv", "MsFraggerPsm" }));
    }

    [Test]
    public void Formats_MostTypesHaveNoUniformViewAtAll()
    {
        int viewless = Run("readers", "formats").GetProperty("formats").EnumerateArray()
            .Count(f => ViewsOf(f).Length == 0);

        Assert.That(viewless, Is.GreaterThan(10),
            "an empty view list is the common case, not an edge case");
    }
}
