using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>library read</c>: an mzLib <c>.msl</c> spectral library, one row per precursor or per
/// fragment ion, through the shared reader engine.
/// </summary>
/// <remarks>
/// <para>
/// mzLib ships no <c>.msl</c> in its test data (its own tests write them), so the libraries here
/// are written the same way, with <see cref="MslLibrary.Save"/>. Two of them are also committed to
/// <c>pkg/python/tests/fixtures/library/</c> as the inputs the recorded Python fixtures were made
/// from; <see cref="CommittedLibrariesAreWhatTheGeneratorWrites"/> holds those bytes to
/// <see cref="LibraryFixtures"/>, so the committed binaries cannot drift from the code that
/// documents what is in them.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class LibraryTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-library-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void RemoveTempDirectory()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    private string Synthetic(string name = "synthetic_library.msl")
    {
        string path = Path.Combine(_tempDirectory, name);
        LibraryFixtures.WriteSynthetic(path);
        return path;
    }

    private string FromProsit()
    {
        string path = Path.Combine(_tempDirectory, "myPrositLib.msl");
        LibraryFixtures.WriteFromMsp(LibraryFixtures.PrositMsp(), path);
        return path;
    }

    // ---- the precursor table -------------------------------------------------------------------

    [Test]
    public void RowsAreEveryPrecursorInStoredOrderNotMzOrder()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic()]);
        MslLibraryEntry[] written = LibraryFixtures.SyntheticEntries();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(written.Length));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(written.Length));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(written.Length));
            Assert.That(Ints(data, "precursor_index"), Is.EqualTo(Enumerable.Range(0, written.Length)));
            Assert.That(Column(data, "full_sequence").Select(v => v.GetString()),
                Is.EqualTo(written.Select(e => e.FullSequence)));
            // The fixture is deliberately not written in m/z order, so this order is the file's.
            double[] mz = Column(data, "precursor_mz").Select(v => v.GetDouble()).ToArray();
            Assert.That(mz, Is.Not.Ordered);
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.False);
        });
    }

    [Test]
    public void DocumentedSentinelsCrossAsNull()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic()]);
        int predicted = LibraryFixtures.PredictedIndex;

        Assert.Multiple(() =>
        {
            // The predicted entry: IonMobility 0.0, QValue NaN, Nce 0 and no protein, all mzLib's
            // documented "not available".
            foreach (string name in new[] { "ion_mobility", "q_value", "nce", "protein_accession", "protein_name", "gene_name" })
                Assert.That(Column(data, name)[predicted].ValueKind, Is.EqualTo(JsonValueKind.Null), name);
            Assert.That(Column(data, "dissociation_type")[predicted].GetString(), Is.EqualTo("Unknown"));
            Assert.That(Column(data, "source_type")[predicted].GetString(), Is.EqualTo("Predicted"));

            // The empirical target carries every one of them.
            Assert.That(Column(data, "ion_mobility")[0].GetDouble(), Is.EqualTo(0.95).Within(1e-6));
            Assert.That(Column(data, "q_value")[0].GetDouble(), Is.EqualTo(0.001).Within(1e-6));
            Assert.That(Column(data, "nce")[0].GetInt32(), Is.EqualTo(28));
            Assert.That(Column(data, "protein_accession")[0].GetString(), Is.EqualTo("P02768"));
            Assert.That(Column(data, "gene_name")[0].GetString(), Is.EqualTo("ALB"));
            Assert.That(Column(data, "source_type")[0].GetString(), Is.EqualTo("Empirical"));
            Assert.That(Column(data, "dissociation_type")[0].GetString(), Is.EqualTo("HCD"));
            Assert.That(Column(data, "is_proteotypic")[0].GetBoolean(), Is.True);
            Assert.That(Strings(data, "absent_fields"), Is.Empty);
        });
    }

    [Test]
    public void RetentionTimeUnitIsPerPrecursor()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic()]);
        string?[] units = Column(data, "retention_time_unit").Select(v => v.GetString()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(units[LibraryFixtures.CalibratedIndex], Is.EqualTo("minutes"));
            Assert.That(units.Where((_, i) => i != LibraryFixtures.CalibratedIndex), Is.All.EqualTo("iRT"));
            Assert.That(Column(data, "retention_time")[LibraryFixtures.CalibratedIndex].GetDouble(),
                Is.EqualTo(21.4).Within(1e-4));
            // The unit is not a per-file fact here, so the block's is null.
            Assert.That(data.GetProperty("retention_time_unit").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void LabelsAreMzLibsOwnTargetDecoyEntrapmentLabels()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic()]);
        string?[] labels = Column(data, "decoy_contaminant_target").Select(v => v.GetString()).ToArray();
        bool[] decoy = Column(data, "is_decoy").Select(v => v.GetBoolean()).ToArray();
        bool[] entrapment = Column(data, "is_entrapment").Select(v => v.GetBoolean()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(labels.Take(4), Is.EqualTo(new[] { "T", "D", "ET", "ED" }));
            Assert.That(labels, Has.None.EqualTo("C"));
            Assert.That(decoy.Take(4), Is.EqualTo(new[] { false, true, false, true }));
            Assert.That(entrapment.Take(4), Is.EqualTo(new[] { false, false, true, true }));
            // MslLibrary.TargetCount counts every non-decoy, entrapment targets included.
            Assert.That(data.GetProperty("decoy_count").GetInt32(), Is.EqualTo(decoy.Count(d => d)));
            Assert.That(data.GetProperty("target_count").GetInt32(), Is.EqualTo(decoy.Count(d => !d)));
        });
    }

    [Test]
    public void EnvelopeCarriesTheLibrarysOwnFacts()
    {
        string path = Synthetic();
        JsonElement data = Invoke(null, ["library", "read", "--path", path]);
        using MslLibrary library = MslLibrary.Load(path);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("path").GetString(), Is.EqualTo(Path.GetFullPath(path)));
            Assert.That(data.GetProperty("file_type").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("reader").GetString(), Is.EqualTo("MslLibrary"));
            Assert.That(data.GetProperty("format_version").GetInt32(), Is.EqualTo(MslFormat.CurrentVersion));
            Assert.That(data.GetProperty("is_compressed").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("is_index_only").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("proteoform_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("min_precursor_mz").GetDouble(), Is.EqualTo(library.MinPrecursorMz).Within(1e-3));
            Assert.That(data.GetProperty("max_precursor_mz").GetDouble(), Is.EqualTo(library.MaxPrecursorMz).Within(1e-3));
            Assert.That(data.GetProperty("rows_not_read").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("fragments_included").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("index_only_threshold").GetInt64(),
                Is.EqualTo(MslFileTypeHandler.DefaultIndexOnlyThresholdBytes));
            Assert.That(data.GetProperty("excluded_fields").GetArrayLength(), Is.Zero);
            Assert.That(Strings(data, "caveats"), Has.Some.Contains("iRT"));
            Assert.That(data.GetProperty("output").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void ALibraryWithNoProteinsOrIonMobilityNamesThoseFieldsAbsent()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", FromProsit()]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data, "absent_fields"),
                Is.EquivalentTo(new[] { "ion_mobility", "protein_accession", "protein_name", "gene_name" }));
            Assert.That(Column(data, "protein_accession"), Has.All.Matches<JsonElement>(v => v.ValueKind == JsonValueKind.Null));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(5));
            Assert.That(Column(data, "source_type").Select(v => v.GetString()), Is.All.EqualTo("Predicted"));
        });
    }

    // ---- load mode -----------------------------------------------------------------------------

    [Test]
    public void IndexOnlyAndFullLoadGiveTheSameRows([Values] bool fragments)
    {
        string path = Synthetic();
        string[] extra = fragments ? ["--fragments"] : [];
        JsonElement full = Invoke(null, ["library", "read", "--path", path, .. extra]);
        JsonElement indexOnly = Invoke(null, ["library", "read", "--path", path, "--index-only-threshold", "0", .. extra]);

        Assert.Multiple(() =>
        {
            Assert.That(full.GetProperty("is_index_only").GetBoolean(), Is.False);
            Assert.That(indexOnly.GetProperty("is_index_only").GetBoolean(), Is.True);
            Assert.That(indexOnly.GetProperty("index_only_threshold").GetInt64(), Is.Zero);
            Assert.That(indexOnly.GetProperty("columns").GetRawText(), Is.EqualTo(full.GetProperty("columns").GetRawText()));
        });
    }

    [Test]
    public void ACompressedLibraryLoadsWholeEvenWhenIndexOnlyIsAskedFor()
    {
        string path = Path.Combine(_tempDirectory, "compressed.msl");
        MslWriter.Write(path, LibraryFixtures.SyntheticEntries(), compressionLevel: 3);

        JsonElement data = Invoke(null, ["library", "read", "--path", path, "--index-only-threshold", "0", "--fragments"]);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("is_compressed").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("is_index_only").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("columns").GetRawText(),
                Is.EqualTo(Invoke(null, ["library", "read", "--path", Synthetic(), "--fragments"]).GetProperty("columns").GetRawText()));
        });
    }

    // ---- fragments -----------------------------------------------------------------------------

    [Test]
    public void FragmentsAreOneRowEachWithThePrecursorRepeated()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic(), "--fragments"]);
        MslLibraryEntry[] written = LibraryFixtures.SyntheticEntries();
        int fragmentRows = written.Sum(e => Math.Max(1, e.MatchedFragmentIons.Count));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("fragments_included").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(written.Length));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(fragmentRows));
            Assert.That(Column(data, "fragment_mz"), Has.Length.EqualTo(fragmentRows));
            Assert.That(Ints(data, "precursor_index").Distinct(), Is.EqualTo(Enumerable.Range(0, written.Length)));
            Assert.That(Strings(data, "column_names").Take(LibraryPrecursorColumnCount(data)),
                Is.EqualTo(Strings(Invoke(null, ["library", "read", "--path", Synthetic()]), "column_names")));
            Assert.That(Strings(data, "caveats"), Has.Some.Contains("custom neutral loss"));
        });
    }

    [Test]
    public void FragmentSentinelsCrossAsNullAndANoLossIsZero()
    {
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic(), "--fragments", "--limit", "1"]);
        var rows = Rows(data);

        JsonElement custom = rows.Single(r => r["neutral_loss"].ValueKind == JsonValueKind.Number
                                              && Math.Abs(r["neutral_loss"].GetDouble() - LibraryFixtures.CustomLoss) < 1e-6)["residue_position"];
        Dictionary<string, JsonElement> internalIon = rows.Single(r => r["is_internal"].GetBoolean());
        Dictionary<string, JsonElement> water = rows.Single(r => Math.Abs(r["neutral_loss"].GetDouble() + 18.010565) < 1e-6);
        var terminal = rows.Where(r => !r["is_internal"].GetBoolean()).ToList();

        Assert.Multiple(() =>
        {
            // The custom loss survives; its residue position does not (MslReader.ConvertFragment).
            Assert.That(custom.ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(internalIon["secondary_product_type"].GetString(), Is.EqualTo("y"));
            Assert.That(internalIon["secondary_fragment_number"].GetInt32(), Is.EqualTo(5));
            Assert.That(internalIon["annotation"].GetString(), Is.EqualTo("bIy[2-5]"));
            Assert.That(terminal.Select(r => r["secondary_fragment_number"].ValueKind), Is.All.EqualTo(JsonValueKind.Null));
            Assert.That(terminal.Select(r => r["secondary_product_type"].ValueKind), Is.All.EqualTo(JsonValueKind.Null));
            Assert.That(water["annotation"].GetString(), Does.StartWith("y3-18.0106"));
            Assert.That(rows.Count(r => r["neutral_loss"].GetDouble() == 0.0), Is.GreaterThan(0),
                "a fragment with no loss carries 0, not null");
            // mzLib's writer scales each precursor's most intense fragment to 1.
            Assert.That(rows.Max(r => r["fragment_intensity"].GetDouble()), Is.EqualTo(1.0).Within(1e-6));
            Assert.That(rows.Select(r => r["fragment_count"].GetInt32()), Is.All.EqualTo(rows.Count));
            // Fragments are as stored: the writer sorts them by m/z.
            Assert.That(rows.Select(r => r["fragment_mz"].GetDouble()), Is.Ordered);
        });
    }

    [Test]
    public void APrecursorWithNoFragmentsIsOneRowWithNullFragmentCells()
    {
        string path = Path.Combine(_tempDirectory, "bare.msl");
        MslLibraryEntry bare = LibraryFixtures.SyntheticEntries()[0];
        bare.MatchedFragmentIons = [];
        MslLibrary.Save(path, [bare, LibraryFixtures.SyntheticEntries()[1]]);

        JsonElement data = Invoke(null, ["library", "read", "--path", path, "--fragments"]);
        var rows = Rows(data);

        Assert.Multiple(() =>
        {
            Assert.That(rows[0]["precursor_index"].GetInt32(), Is.Zero);
            Assert.That(rows[0]["fragment_count"].GetInt32(), Is.Zero);
            Assert.That(rows[0]["fragment_mz"].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(rows[0]["product_type"].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(rows.Count(r => r["precursor_index"].GetInt32() == 0), Is.EqualTo(1));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
        });
    }

    // ---- the window and --out ------------------------------------------------------------------

    [Test]
    public void OffsetAndLimitCountPrecursors([Values] bool fragments)
    {
        string[] extra = fragments ? ["--fragments"] : [];
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic(), "--offset", "1", "--limit", "2", .. extra]);

        Assert.Multiple(() =>
        {
            Assert.That(Ints(data, "precursor_index").Distinct(), Is.EqualTo(new[] { 1, 2 }));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("offset").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(LibraryFixtures.SyntheticEntries().Length));
        });
    }

    [Test]
    public void OutWritesTheTableAndReturnsNoColumns([Values] bool fragments)
    {
        string out_ = Path.Combine(_tempDirectory, "nested", "library.tsv");
        string[] extra = fragments ? ["--fragments"] : [];
        JsonElement inline = Invoke(null, ["library", "read", "--path", Synthetic(), .. extra]);
        JsonElement data = Invoke(null, ["library", "read", "--path", Synthetic(), "--out", out_, .. extra]);
        string[] lines = File.ReadAllLines(out_);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("output").GetProperty("row_count").GetInt32(),
                Is.EqualTo(inline.GetProperty("row_count").GetInt32()));
            Assert.That(lines[0].Split('\t'), Is.EqualTo(Strings(data, "column_names")));
            Assert.That(lines, Has.Length.EqualTo(inline.GetProperty("row_count").GetInt32() + 1));
        });
    }

    // ---- many libraries ------------------------------------------------------------------------

    [Test]
    public void ManyLibrariesAreByteIdenticalAtOneAndFourThreads([Values] bool fragments)
    {
        string[] paths = [Synthetic(), FromProsit(), Synthetic("another.msl")];
        string[] extra = fragments ? ["--fragments"] : [];

        string one = Raw(Stdin(paths), ["library", "read", "--paths-stdin", "--threads", "1", .. extra]);
        string four = Raw(Stdin(paths), ["library", "read", "--paths-stdin", "--threads", "4", .. extra]);

        Assert.That(four, Is.EqualTo(one));
        JsonElement data = JsonDocument.Parse(one).RootElement.GetProperty("data");
        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("read_count").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(2 * LibraryFixtures.SyntheticEntries().Length + 5));
            Assert.That(Strings(data, "column_names").Take(2), Is.EqualTo(new[] { "source_index", "source_path" }));
            Assert.That(Ints(data, "source_index").Distinct(), Is.EqualTo(new[] { 0, 1, 2 }));
            JsonElement prosit = data.GetProperty("files")[1];
            Assert.That(prosit.GetProperty("format_version").GetInt32(), Is.EqualTo(MslFormat.CurrentVersion));
            Assert.That(prosit.GetProperty("record_count").GetInt32(), Is.EqualTo(5));
            Assert.That(prosit.GetProperty("absent_fields").GetArrayLength(), Is.EqualTo(4));
            Assert.That(data.GetProperty("fragments_included").GetBoolean(), Is.EqualTo(fragments));
        });
    }

    [Test]
    public void ASkippedLibraryCarriesEveryKeyOfARead()
    {
        string missing = Path.Combine(_tempDirectory, "missing.msl");
        JsonElement data = Invoke(Stdin([Synthetic(), missing]), ["library", "read", "--paths-stdin", "--on-error", "skip"]);
        JsonElement read = data.GetProperty("files")[0];
        JsonElement failed = data.GetProperty("files")[1];

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(failed.EnumerateObject().Select(p => p.Name), Is.EqualTo(read.EnumerateObject().Select(p => p.Name)));
            Assert.That(failed.GetProperty("error").GetProperty("kind").GetString(), Is.EqualTo("usage"));
            Assert.That(failed.GetProperty("format_version").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    // ---- refusals ------------------------------------------------------------------------------

    private static IEnumerable<TestCaseData> UsageCases()
    {
        yield return new TestCaseData((object)new[] { "--fragments", "true" }, "takes no value");
        yield return new TestCaseData((object)new[] { "--index-only-threshold", "-1" }, "zero or greater");
        yield return new TestCaseData((object)new[] { "--index-only-threshold", "1GB" }, "whole number of bytes");
        yield return new TestCaseData((object)new[] { "--index-only-threshold" }, "no value");
        yield return new TestCaseData((object)new[] { "--limit", "-1" }, "--limit");
        yield return new TestCaseData((object)new[] { "--out", "library.csv" }, ".tsv");
        yield return new TestCaseData((object)new[] { "--on-error", "skip" }, "skip");
    }

    [Test]
    [TestCaseSource(nameof(UsageCases))]
    public void BadOptionsAreUsageErrors(string[] extra, string expected)
    {
        JsonElement error = InvokeExpectingError(null, ["library", "read", "--path", Synthetic(), .. extra]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    [Test]
    public void AMissingFileADirectoryAndAnotherExtensionAreUsageErrors()
    {
        string msp = Path.Combine(_tempDirectory, "library.msp");
        File.WriteAllText(msp, "Name: PEPTIDE/2\n");

        JsonElement missing = InvokeExpectingError(null, ["library", "read", "--path", Path.Combine(_tempDirectory, "none.msl")]);
        JsonElement directory = InvokeExpectingError(null, ["library", "read", "--path", _tempDirectory]);
        JsonElement other = InvokeExpectingError(null, ["library", "read", "--path", msp]);
        JsonElement noPath = InvokeExpectingError(null, ["library", "read"]);

        Assert.Multiple(() =>
        {
            Assert.That(missing.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(missing.GetProperty("message").GetString(), Does.Contain("not found"));
            Assert.That(directory.GetProperty("message").GetString(), Does.Contain("directory"));
            Assert.That(other.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(other.GetProperty("message").GetString(), Does.Contain(".msl"));
            Assert.That(noPath.GetProperty("type").GetString(), Is.EqualTo("usage"));
        });
    }

    [Test]
    public void AFileThatIsNotAnMslIsMzLibsCorrectnessError()
    {
        string path = Path.Combine(_tempDirectory, "corrupt.msl");
        File.WriteAllBytes(path, Enumerable.Repeat((byte)7, 200).ToArray());

        JsonElement error = InvokeExpectingError(null, ["library", "read", "--path", path]);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo(nameof(FormatException)));
    }

    [Test]
    public void AnEmptyLibraryHasNoMzRangeAndNoRows()
    {
        string path = Path.Combine(_tempDirectory, "empty.msl");
        MslLibrary.Save(path, []);

        JsonElement data = Invoke(null, ["library", "read", "--path", path, "--fragments"]);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("min_precursor_mz").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("max_precursor_mz").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(Column(data, "precursor_index"), Is.Empty);
        });
    }

    // ---- the committed Python inputs -----------------------------------------------------------

    /// <summary>
    /// The two <c>.msl</c> files the recorded Python fixtures were read from are exactly what
    /// <see cref="LibraryFixtures"/> writes. Set <c>PYMZLIB_WRITE_LIBRARY_FIXTURES=1</c> to rewrite
    /// them after changing the generator, then re-record the fixtures.
    /// </summary>
    [Test]
    public void CommittedLibrariesAreWhatTheGeneratorWrites()
    {
        string directory = Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "library");
        (string Name, Action<string> Write)[] libraries =
        [
            ("synthetic_library.msl", LibraryFixtures.WriteSynthetic),
            ("myPrositLib.msl", path => LibraryFixtures.WriteFromMsp(LibraryFixtures.PrositMsp(), path)),
        ];

        foreach ((string name, Action<string> write) in libraries)
        {
            string committed = Path.Combine(directory, name);
            if (Environment.GetEnvironmentVariable("PYMZLIB_WRITE_LIBRARY_FIXTURES") == "1")
            {
                Directory.CreateDirectory(directory);
                write(committed);
            }

            string fresh = Path.Combine(_tempDirectory, name);
            write(fresh);
            Assert.That(File.Exists(committed), Is.True, $"{committed} is missing; see this test's summary");
            Assert.That(File.ReadAllBytes(committed), Is.EqualTo(File.ReadAllBytes(fresh)),
                $"{name} is not what LibraryFixtures writes; regenerate it (see this test's summary)");
        }
    }

    // ---- harness -------------------------------------------------------------------------------

    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    private static int LibraryPrecursorColumnCount(JsonElement fragmentData) =>
        Strings(fragmentData, "column_names").ToList().IndexOf("fragment_count");

    private static string Stdin(IEnumerable<string> paths) => string.Join('\n', paths) + "\n";

    private static string[] Strings(JsonElement data, string name) =>
        data.GetProperty(name).EnumerateArray().Select(v => v.GetString()!).ToArray();

    private static JsonElement[] Column(JsonElement data, string name) =>
        data.GetProperty("columns").GetProperty(name).EnumerateArray().ToArray();

    private static int[] Ints(JsonElement data, string name) => Column(data, name).Select(v => v.GetInt32()).ToArray();

    private static List<Dictionary<string, JsonElement>> Rows(JsonElement data)
    {
        string[] names = Strings(data, "column_names");
        int count = data.GetProperty("row_count").GetInt32();
        return Enumerable.Range(0, count)
            .Select(i => names.ToDictionary(n => n, n => Column(data, n)[i]))
            .ToList();
    }

    private static JsonElement Invoke(string? stdin, string[] args)
    {
        JsonElement envelope = JsonDocument.Parse(Raw(stdin, args)).RootElement;
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.True, $"Expected success, got: {envelope}");
        return envelope.GetProperty("data").Clone();
    }

    private static JsonElement InvokeExpectingError(string? stdin, string[] args)
    {
        JsonElement envelope = JsonDocument.Parse(Raw(stdin, args)).RootElement;
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.False, $"Expected a failure, got: {envelope}");
        return envelope.GetProperty("error").Clone();
    }

    /// <summary>Runs a verb as Main would and returns the envelope text exactly as it would be written.</summary>
    private static string Raw(string? stdin, string[] args)
    {
        TextReader previousIn = Console.In;
        if (stdin is not null)
            Console.SetIn(new StringReader(stdin));
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.Serialize(new Program.Envelope { Ok = true, Data = data }, Program.JsonOptions);
        }
        catch (Exception exception)
        {
            string type = Program.Unwrap(exception is Program.SourceFailure failure ? failure.InnerException! : exception)
                is Program.UsageException ? "usage" : Program.ClassifyError(exception);
            return JsonSerializer.Serialize(
                new Program.Envelope
                {
                    Ok = false,
                    Error = new Program.ErrorInfo { Type = type, Message = Program.Unwrap(exception).Message },
                },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetIn(previousIn);
        }
    }
}

/// <summary>
/// The spectral libraries <see cref="LibraryTests"/> reads, and that the committed Python inputs
/// were written from.
/// </summary>
/// <remarks>
/// <para>
/// <b>synthetic_library.msl</b>: six precursors that exercise every field <c>library read</c>
/// projects, written deliberately out of m/z order. In order: an empirical target with a protein,
/// ion mobility, q-value and NCE, and fragments that include a water loss, a custom loss and an
/// internal ion; its decoy; an entrapment target whose retention time is calibrated minutes; an
/// entrapment decoy; a predicted precursor with every "not available" sentinel mzLib documents; and
/// a proteoform.
/// </para>
/// <para>
/// <b>myPrositLib.msl</b>: mzLib's own test library <c>myPrositLib.msp</c> (five Prosit-predicted
/// precursors) converted by mzLib's <see cref="MslConverter.FromMspFile"/> and saved with
/// <see cref="MslLibrary.Save"/>, precursors in the text file's order.
/// </para>
/// </remarks>
[ExcludeFromCodeCoverage]
internal static class LibraryFixtures
{
    /// <summary>The custom neutral loss on the target's y4 ion, in Da (not one of the format's named losses).</summary>
    public const double CustomLoss = -43.005814;

    /// <summary>The precursor whose retention time is calibrated minutes rather than iRT.</summary>
    public const int CalibratedIndex = 2;

    /// <summary>The predicted precursor that carries every "not available" sentinel.</summary>
    public const int PredictedIndex = 4;

    public static string PrositMsp([CallerFilePath] string thisFile = "")
    {
        string path = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "..", "code", "mzLib", "mzLib",
            "Test", "FileReadingTests", "SpectralLibraryTests", "SpectralLibraryData", "myPrositLib.msp"));
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present in the worktree: {path}");
        return path;
    }

    public static void WriteSynthetic(string path) => MslLibrary.Save(path, SyntheticEntries());

    public static void WriteFromMsp(string mspPath, string path)
    {
        using MslLibrary converted = MslConverter.FromMspFile(mspPath);
        var entries = Enumerable.Range(0, converted.PrecursorCount).Select(i => converted.GetEntry(i)!).ToList();
        MslLibrary.Save(path, entries);
    }

    /// <summary>Fresh entries every call: <see cref="MslLibrary.Save"/> sorts and rescales fragments in place.</summary>
    public static MslLibraryEntry[] SyntheticEntries() =>
    [
        Entry("PEPTIDEK", 465.7400, 35.2, decoy: false, entrapment: false, ionMobility: 0.95, qValue: 0.001f,
            protein: ("P02768", "Albumin", "ALB"), proteotypic: true, nce: 28, source: MslFormat.SourceType.Empirical,
            dissociation: DissociationType.HCD, fragments:
            [
                Ion(147.1128f, 0.40f, ProductType.y, 1, 8),
                Ion(276.1554f, 1.00f, ProductType.y, 2, 7),
                Ion(227.1026f, 0.30f, ProductType.b, 2, 2),
                Ion(373.2082f, 0.20f, ProductType.y, 3, 6, neutralLoss: -18.010565),
                Ion(461.2242f, 0.15f, ProductType.y, 4, 5, neutralLoss: CustomLoss),
                new MslFragmentIon
                {
                    Mz = 326.1710f, Intensity = 0.10f, ProductType = ProductType.b, SecondaryProductType = ProductType.y,
                    FragmentNumber = 2, SecondaryFragmentNumber = 5, ResiduePosition = 2, Charge = 1,
                },
            ]),
        Entry("KEDITPEP", 465.7401, 35.2, decoy: true, entrapment: false, ionMobility: 0.95, qValue: 0.5f,
            protein: ("DECOY_P02768", "Albumin", "ALB"), proteotypic: false, nce: 28, source: MslFormat.SourceType.Empirical,
            dissociation: DissociationType.HCD, fragments: [Ion(147.1128f, 1.0f, ProductType.y, 1, 8), Ion(244.1292f, 0.5f, ProductType.b, 2, 2)]),
        Entry("PEPTLDEK", 458.7462, 21.4, decoy: false, entrapment: true, ionMobility: 0.93, qValue: 0.002f,
            protein: ("ENTRAP_P02768", "Albumin", "ALB"), proteotypic: true, nce: 28, source: MslFormat.SourceType.EmpiricalRefined,
            dissociation: DissociationType.HCD, fragments: [Ion(147.1128f, 1.0f, ProductType.y, 1, 8)], rtIsCalibrated: true),
        Entry("KEDLTPEP", 458.7463, 33.0, decoy: true, entrapment: true, ionMobility: 0.93, qValue: 0.6f,
            protein: ("DECOY_ENTRAP_P02768", "Albumin", "ALB"), proteotypic: false, nce: 28, source: MslFormat.SourceType.Empirical,
            dissociation: DissociationType.HCD, fragments: [Ion(147.1128f, 1.0f, ProductType.y, 1, 8)]),
        Entry("ELVISLIVESK", 406.2420, 61.7, decoy: false, entrapment: false, ionMobility: 0.0, qValue: float.NaN,
            protein: (string.Empty, string.Empty, string.Empty), proteotypic: false, nce: 0, source: MslFormat.SourceType.Predicted,
            dissociation: DissociationType.Unknown, fragments: [Ion(147.1128f, 1.0f, ProductType.y, 1, 11), Ion(243.1339f, 0.6f, ProductType.b, 2, 2)],
            charge: 3),
        Entry("MAHVRSGKHQIPAEALQ", 952.4930, 48.0, decoy: false, entrapment: false, ionMobility: 0.0, qValue: 0.01f,
            protein: ("P04637", "Cellular tumor antigen p53", "TP53"), proteotypic: true, nce: 30, source: MslFormat.SourceType.Predicted,
            dissociation: DissociationType.ETD, fragments: [Ion(132.0478f, 1.0f, ProductType.c, 1, 1), Ion(203.0849f, 0.8f, ProductType.c, 2, 2)],
            charge: 2, molecule: MslFormat.MoleculeType.Proteoform),
    ];

    private static MslFragmentIon Ion(float mz, float intensity, ProductType type, int number, int residue, double neutralLoss = 0.0) =>
        new()
        {
            Mz = mz, Intensity = intensity, ProductType = type, FragmentNumber = number,
            ResiduePosition = residue, Charge = 1, NeutralLoss = neutralLoss,
        };

    private static MslLibraryEntry Entry(
        string sequence, double mz, double rt, bool decoy, bool entrapment, double ionMobility, float qValue,
        (string Accession, string Name, string Gene) protein, bool proteotypic, int nce, MslFormat.SourceType source,
        DissociationType dissociation, List<MslFragmentIon> fragments, bool rtIsCalibrated = false, int charge = 2,
        MslFormat.MoleculeType molecule = MslFormat.MoleculeType.Peptide) =>
        new()
        {
            FullSequence = sequence,
            BaseSequence = sequence,
            PrecursorMz = mz,
            ChargeState = charge,
            RetentionTime = rt,
            RtIsCalibrated = rtIsCalibrated,
            IsDecoy = decoy,
            IsEntrapment = entrapment,
            IonMobility = ionMobility,
            QValue = qValue,
            ProteinAccession = protein.Accession,
            ProteinName = protein.Name,
            GeneName = protein.Gene,
            IsProteotypic = proteotypic,
            Nce = nce,
            Source = source,
            DissociationType = dissociation,
            MoleculeType = molecule,
            MatchedFragmentIons = fragments,
        };
}
