using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
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
            "readers read-results" => Reading.ReadResults(arguments),
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
    public void Identify_UnsupportedExtension_IsAUsageErrorPointingAtTheFormatsListing()
    {
        var exception = Assert.Throws<Program.UsageException>(
            () => Run("readers", "identify", "--path", Touch("notes.docx")));

        // Points at the capability, not at a command: the wire contract is language-neutral, and
        // telling a Python caller to "run readers formats" names a CLI they do not have.
        Assert.That(exception!.Message, Does.Contain("formats listing"),
            "an unsupported file is a bad argument, and the caller needs to know what IS supported");
        Assert.That(exception.Message, Does.Not.Contain("Run '"),
            "the bridge must not prescribe a command; its callers are not all shells");
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

    // ---- read-results ------------------------------------------------------------------------

    /// <summary>
    /// The reader fixtures in the pinned mzLib worktree, located relative to this source file so
    /// they resolve both locally and in CI.
    /// </summary>
    private static string FixtureDirectory([CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "code", "mzLib", "mzLib", "Test", "FileReadingTests");
    }

    /// <summary>A named mzLib fixture, or an ignored test when the worktree is absent.</summary>
    private static string Fixture(params string[] parts)
    {
        string path = Path.Combine(new[] { FixtureDirectory() }.Concat(parts).ToArray());
        if (!File.Exists(path))
            Assert.Ignore($"mzLib reader fixture not present in the worktree: {path}");
        return path;
    }

    private static string Psmtsv() => Fixture("SearchResults", "BottomUpExample.psmtsv");

    private static string MsFragger() =>
        Fixture("ExternalFileTypes", "FraggerPsm_FragPipev21.1_psm.tsv");

    [Test]
    public void ReadResults_ReturnsEveryRecordColumnarByDefault()
    {
        JsonElement result = Run("readers", "read-results", "--path", Psmtsv());

        int count = result.GetProperty("record_count").GetInt32();
        Assert.That(count, Is.GreaterThan(0));
        Assert.That(result.GetProperty("returned_count").GetInt32(), Is.EqualTo(count),
            "there is no default row cap — the ordinary call returns the whole file");
        Assert.That(result.GetProperty("truncated").GetBoolean(), Is.False);

        // Columnar: one array per field, each as long as the record count.
        JsonElement columns = result.GetProperty("columns");
        foreach (JsonElement name in result.GetProperty("column_names").EnumerateArray())
            Assert.That(columns.GetProperty(name.GetString()!).GetArrayLength(), Is.EqualTo(count));
    }

    [Test]
    public void ReadResults_Limit_TruncatesAndDisclosesIt()
    {
        JsonElement result = Run("readers", "read-results", "--path", Psmtsv(), "--limit", "2");

        Assert.That(result.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
        Assert.That(result.GetProperty("record_count").GetInt32(), Is.GreaterThan(2));
        Assert.That(result.GetProperty("truncated").GetBoolean(), Is.True,
            "a short answer and a complete one must not look alike");
    }

    [Test]
    public void ReadResults_OffsetPastTheEnd_ReturnsNothingAndStillReportsTheTotal()
    {
        JsonElement result = Run("readers", "read-results", "--path", Psmtsv(), "--offset", "10000");

        Assert.That(result.GetProperty("returned_count").GetInt32(), Is.Zero);
        Assert.That(result.GetProperty("record_count").GetInt32(), Is.GreaterThan(0));
        Assert.That(result.GetProperty("truncated").GetBoolean(), Is.True,
            "records were skipped, so the answer is incomplete even though the limit never bit");
    }

    [Test]
    public void ReadResults_MsFragger_DisclosesThatRetentionTimeIsSeconds()
    {
        // The single most consequential caveat: these values are seconds while MetaMorpheus's are
        // minutes, and nothing in mzLib converts them. A caller comparing the two silently gets a
        // 60x error, so the wire must say so.
        string[] caveats = Run("readers", "read-results", "--path", MsFragger())
            .GetProperty("caveats").EnumerateArray().Select(c => c.GetString()!).ToArray();

        Assert.That(caveats.Any(c => c.Contains("SECONDS")), Is.True);
        Assert.That(caveats.Any(c => c.Contains("is_decoy")), Is.True);
        Assert.That(caveats.Any(c => c.Contains("THEORETICAL")), Is.True);
    }

    [Test]
    public void ReadResults_Psmtsv_DisclosesAmbiguityAndSilentRowLoss()
    {
        string[] caveats = Run("readers", "read-results", "--path", Psmtsv())
            .GetProperty("caveats").EnumerateArray().Select(c => c.GetString()!).ToArray();

        Assert.That(caveats, Is.Not.Empty);
        Assert.That(caveats.Any(c => c.Contains("FIRST candidate")), Is.True);
    }

    [Test]
    public void ReadResults_ReportsThatEveryRowWasRead()
    {
        // mzLib drops a malformed row silently; this counts the difference so a partial read is
        // visible rather than passing for a complete one.
        JsonElement result = Run("readers", "read-results", "--path", Psmtsv());

        Assert.That(result.GetProperty("rows_not_read").GetInt32(), Is.Zero);
    }

    [Test]
    public void ReadResults_Out_WritesATabSeparatedTableAndOmitsTheInlineRecords()
    {
        string destination = Path.Combine(_tempDirectory, "records.tsv");

        JsonElement result = Run("readers", "read-results", "--path", Psmtsv(), "--out", destination);

        Assert.That(result.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null),
            "writing to disk exists to keep the table OUT of the envelope");
        Assert.That(result.GetProperty("returned_count").GetInt32(), Is.Zero);

        JsonElement output = result.GetProperty("output");
        Assert.That(output.GetProperty("format").GetString(), Is.EqualTo("tsv"));
        Assert.That(File.Exists(destination), Is.True);

        string[] lines = File.ReadAllLines(destination);
        int columnCount = result.GetProperty("column_names").GetArrayLength();
        Assert.That(lines[0].Split('\t'), Has.Length.EqualTo(columnCount));
        Assert.That(lines[0], Does.Not.Contain(","),
            "tab-separated, not comma-separated: these fields contain commas");
        Assert.That(lines, Has.Length.EqualTo(output.GetProperty("row_count").GetInt32() + 1));
    }

    [Test]
    public void ReadResults_Out_CreatesMissingDirectories()
    {
        string destination = Path.Combine(_tempDirectory, "nested", "deeper", "records.tsv");

        Run("readers", "read-results", "--path", Psmtsv(), "--out", destination);

        Assert.That(File.Exists(destination), Is.True);
    }

    [Test]
    public void ReadResults_FormatWithoutTheView_IsAUsageErrorNamingWhatItDoesHave()
    {
        string toppic = Fixture("ExternalFileTypes", "ToppicPrsm_TopPICv1.6.2_prsm.tsv");

        var exception = Assert.Throws<Program.UsageException>(
            () => Run("readers", "read-results", "--path", toppic));

        // mzLib throws the same exception for "unknown extension" and "wrong interface"; the
        // caller deserves to know which, and what the file can do instead.
        Assert.That(exception!.Message, Does.Contain("ToppicPrsm"));
        Assert.That(exception.Message, Does.Contain("no cross-format record view"));
        Assert.That(exception.Message, Does.Not.Contain("Run '"),
            "the bridge must not prescribe a command; its callers are not all shells");
    }

    [Test]
    public void ReadResults_SpectraFile_IsAUsageErrorNamingTheSpectraView()
    {
        var exception = Assert.Throws<Program.UsageException>(
            () => Run("readers", "read-results", "--path", Touch("run.mzML")));

        Assert.That(exception!.Message, Does.Contain("spectra"));
    }

    [Test]
    public void ReadResults_MissingFile_IsAUsageError() =>
        Assert.Throws<Program.UsageException>(() => Run(
            "readers", "read-results", "--path", Path.Combine(_tempDirectory, "absent.psmtsv")));

    [Test]
    public void ReadResults_NegativeLimit_IsAUsageError() =>
        Assert.Throws<Program.UsageException>(() => Run(
            "readers", "read-results", "--path", Psmtsv(), "--limit", "-1"));

    [Test]
    public void ReadResults_NegativeOffset_IsAUsageError() =>
        Assert.Throws<Program.UsageException>(() => Run(
            "readers", "read-results", "--path", Psmtsv(), "--offset", "-5"));

    [Test]
    public void ReadResults_ReportsRetentionTimeUnitPerFormat()
    {
        // The units differ by format and mzLib normalises nothing, so the unit must cross as a
        // VALUE. Reported in prose only, a caller has to grep a sentence for "SECONDS" - which is
        // what one did before this field existed.
        Assert.That(Run("readers", "read-results", "--path", Psmtsv())
            .GetProperty("retention_time_unit").GetString(), Is.EqualTo("minutes"));

        Assert.That(Run("readers", "read-results", "--path", MsFragger())
            .GetProperty("retention_time_unit").GetString(), Is.EqualTo("seconds"));
    }

    [Test]
    public void ReadResults_MassCaveat_DoesNotClaimPsmtsvReportsObservedMass()
    {
        // Regression for a caveat that was WRONG and actively misled a reader: it said the psmtsv
        // formats report the observed mass. They do not - they report "Peptide Monoisotopic Mass",
        // the theoretical value, exactly as MSFragger does. Verified on BottomUpExample.psmtsv,
        // where record 1 is 1959.90366 (theoretical) and NOT 1959.9122 (the precursor mass).
        string[] caveats = Run("readers", "read-results", "--path", MsFragger())
            .GetProperty("caveats").EnumerateArray().Select(c => c.GetString()!).ToArray();

        Assert.That(caveats.Any(c => c.Contains("psmtsv formats report the observed mass")), Is.False,
            "both formats report the THEORETICAL mass; claiming otherwise invents a discrepancy");

        double firstMass = Run("readers", "read-results", "--path", Psmtsv(), "--limit", "1")
            .GetProperty("columns").GetProperty("monoisotopic_mass")[0].GetDouble();
        Assert.That(firstMass, Is.EqualTo(1959.90366).Within(1e-5),
            "psmtsv monoisotopic_mass is the theoretical Peptide Monoisotopic Mass column");
    }

    [Test]
    public void ReadResults_ProteinColumnsUseMzLibsOwnTupleNames()
    {
        // mzLib's tuple is (proteinAccessions, geneName, organism). Pluralising the last two here
        // would force every reader to hold a translation table against the mzLib source.
        string[] names = Run("readers", "read-results", "--path", Psmtsv())
            .GetProperty("column_names").EnumerateArray().Select(n => n.GetString()!).ToArray();

        Assert.That(names, Does.Contain("protein_accessions"));
        Assert.That(names, Does.Contain("gene_name"));
        Assert.That(names, Does.Contain("organism"));
    }

    [Test]
    public void ReadResults_IsDecoyIsNullWhereTheFormatCannotReportDecoys()
    {
        // MSFragger psm.tsv has no target/decoy column; mzLib hardcodes false. A boolean that
        // silently means "unknown" for one format and "target" for another is the worst thing this
        // view could hand back, and both bake-off rounds named it as such.
        JsonElement fragger = Run("readers", "read-results", "--path", MsFragger());
        foreach (JsonElement value in fragger.GetProperty("columns").GetProperty("is_decoy").EnumerateArray())
            Assert.That(value.ValueKind, Is.EqualTo(JsonValueKind.Null));

        // The psmtsv family genuinely reads decoys, so it must still report real booleans.
        JsonElement psmtsv = Run("readers", "read-results", "--path", Psmtsv());
        foreach (JsonElement value in psmtsv.GetProperty("columns").GetProperty("is_decoy").EnumerateArray())
            Assert.That(value.ValueKind, Is.EqualTo(JsonValueKind.False).Or.EqualTo(JsonValueKind.True));
    }

    [Test]
    public void Caveats_AreAsciiSoTheyReadCorrectlyInAConsole()
    {
        // These strings are printed by humans on Windows lab machines, where a cp1252 console turns
        // an em-dash into a replacement character and makes every caveat look corrupted.
        foreach (string path in new[] { Psmtsv(), MsFragger() })
        {
            foreach (JsonElement caveat in Run("readers", "read-results", "--path", path)
                         .GetProperty("caveats").EnumerateArray())
            {
                string text = caveat.GetString()!;
                Assert.That(text.All(c => c < 128), Is.True, $"non-ASCII in caveat: {text}");
            }
        }
    }

    [Test]
    public void ReadResults_ZeroLimit_ReturnsNoRecordsButStillReportsTheTotal()
    {
        JsonElement result = Run("readers", "read-results", "--path", Psmtsv(), "--limit", "0");

        Assert.That(result.GetProperty("returned_count").GetInt32(), Is.Zero);
        Assert.That(result.GetProperty("record_count").GetInt32(), Is.GreaterThan(0));
        Assert.That(result.GetProperty("truncated").GetBoolean(), Is.True);
    }
}
