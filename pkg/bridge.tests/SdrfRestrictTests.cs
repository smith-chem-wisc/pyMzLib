using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Readers;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>sdrf restrict</c>, the projection of mzLib's <c>SdrfSearchScope.Restrict</c>.
/// </summary>
/// <remarks>
/// <para>
/// The document is mzLib's <b>PXD008841.subset</b> (40 rows, 4 files, TMT10; see
/// <see cref="SdrfIsobaricDesignTests"/>). <see cref="Messy"/> is a search that exercises every
/// report but one: two calibrated derivatives, the raw file of one of them as well (ambiguous), an
/// mzML whose acquisition the caller names (sample_A.mzML from pool 5's fraction 1), and a file the
/// SDRF does not describe. <see cref="Clean"/> is the four files, all calibrated. shared_stem needs a
/// deposit that lists one run twice, which no mzLib fixture does, so it is built here.
/// </para>
/// <para>
/// <b>PXD008841_searched.sdrf.tsv</b> is the clean restriction, written by this verb's <c>--out</c>
/// (mzLib's <c>WriteResults</c>) and committed as the pyMzLib guide's input for the restrict-then-design
/// example; <see cref="Restrict_TheCommittedSearchedDocumentIsWhatTheVerbWrites"/> holds it to the
/// verb.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SdrfRestrictTests
{
    private const string Pxd008841 = "PXD008841.subset.sdrf.tsv";
    private const string Pool1 = "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr0";
    private const string Pool5 = "HJOSLO2U_QEHF_20150329_TMT_pool5_300ugIPG37-49_7of15ul_fr0";

    private static readonly string Messy = string.Join("\n",
        $@"D:\search\{Pool1}1-calib.mzML",
        $"{Pool1}2-calib.mzML",
        $"{Pool1}2.raw",
        $"sample_A.mzML\t{Pool5}1.raw",
        "not_in_sdrf.mzML");

    private static readonly string[] CleanFiles =
    {
        $@"D:\search\{Pool1}1-calib.mzML", $@"D:\search\{Pool1}2-calib.mzML",
        $@"D:\search\{Pool5}1-calib.mzML", $@"D:\search\{Pool5}2-calib.mzML",
    };

    private static readonly string Clean = string.Join("\n", CleanFiles);

    // The same four files under this platform's own separator. restrict reduces a path written with
    // either separator on any platform, but the design readers use Path.GetFileName, which splits
    // only this platform's, so a design test given D:\search\... on Linux keys on the whole string
    // (SDRF-B5 asks mzLib for one rule).
    private static readonly string[] NativeCleanFiles = CleanFiles
        .Select(f => Path.Combine(Path.GetTempPath(), "search", f[(f.LastIndexOf('\\') + 1)..]))
        .ToArray();

    private static readonly string NativeClean = string.Join("\n", NativeCleanFiles);

    [Test]
    public void Restrict_JoinsDerivativesAndReportsWhatItCouldNotPlace()
    {
        JsonElement data = InvokeWithStdin(Messy, "sdrf", "restrict", "--path", Fixture(Pxd008841));

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data.GetProperty("searched_files")), Has.Length.EqualTo(5));
            Assert.That(data.GetProperty("acquired_name_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("input_row_count").GetInt32(), Is.EqualTo(40));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(30), "three of the four files, ten channels each");
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(30));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.False);
            Assert.That(Strings(data.GetProperty("searched_without_row")), Is.EqualTo(new[] { "not_in_sdrf.mzML" }));
            Assert.That(Strings(data.GetProperty("dropped_data_files")), Is.EqualTo(new[] { $"{Pool5}2.raw" }));
            Assert.That(Strings(data.GetProperty("ambiguous")),
                Is.EqualTo(new[] { $"{Pool1}2: {Pool1}2-calib.mzML, {Pool1}2.raw" }),
                "mzLib's 'stem: searched, searched'; the first in name order kept the rows");
            Assert.That(Strings(data.GetProperty("shared_stem")), Is.Empty);
            Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });

        string[] header = Strings(data.GetProperty("column_names"));
        int dataFile = Array.IndexOf(header, "comment[data file]");
        Assert.That(header[dataFile + 1], Is.EqualTo("comment[searched data file]"), "directly after comment[data file]");

        var pairs = data.GetProperty("rows").EnumerateArray()
            .Select(r => (Acquired: r[dataFile].GetString(), Searched: r[dataFile + 1].GetString()))
            .Distinct().ToArray();
        Assert.That(pairs, Is.EqualTo(new[]
        {
            ($"{Pool1}1.raw", $"{Pool1}1-calib.mzML"),
            ($"{Pool1}2.raw", $"{Pool1}2-calib.mzML"),
            ($"{Pool5}1.raw", "sample_A.mzML"),
        }), "the searched file's NAME, no folder; a known acquisition wins over the suffix rule");
    }

    [Test]
    public void Restrict_IsMzLibsRestrictionCellForCell()
    {
        JsonElement data = InvokeWithStdin(Messy, "sdrf", "restrict", "--path", Fixture(Pxd008841));

        SdrfSearchScoping mzLib = SdrfSearchScope.Restrict(new SdrfDocument(Fixture(Pxd008841)),
            new[] { $@"D:\search\{Pool1}1-calib.mzML", $"{Pool1}2-calib.mzML", $"{Pool1}2.raw", "sample_A.mzML", "not_in_sdrf.mzML" },
            new Dictionary<string, string> { ["sample_A.mzML"] = $"{Pool5}1.raw" });

        Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(mzLib.Document.Header.ToArray()));
        string[][] rows = data.GetProperty("rows").EnumerateArray().Select(Strings).ToArray();
        Assert.That(rows, Is.EqualTo(mzLib.Document.Results.Select(r => r.Cells.ToArray()).ToArray()));
        Assert.That(Strings(data.GetProperty("searched_without_row")), Is.EqualTo(mzLib.SearchedWithoutRow));
        Assert.That(Strings(data.GetProperty("dropped_data_files")), Is.EqualTo(mzLib.DroppedDataFiles));
        Assert.That(Strings(data.GetProperty("ambiguous")), Is.EqualTo(mzLib.Ambiguous));
        Assert.That(Strings(data.GetProperty("shared_stem")), Is.EqualTo(mzLib.SharedStem));
    }

    [Test]
    public void Restrict_OneRunListedTwiceIsKeptAndReportedAsASharedStem()
    {
        string path = TempSdrf(
            new[] { "S1", "1", "run 1", "label free sample", "X.raw" },
            new[] { "S1", "1", "run 1 mzML", "label free sample", "X.mzML" },
            new[] { "S2", "2", "run 2", "label free sample", "Y.raw" });
        try
        {
            JsonElement data = InvokeWithStdin("X-calib.mzML", "sdrf", "restrict", "--path", path);

            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(2), "both rows of the run are kept");
            Assert.That(Strings(data.GetProperty("shared_stem")), Is.EqualTo(new[] { "X-calib.mzML: X.raw, X.mzML" }));
            Assert.That(Strings(data.GetProperty("dropped_data_files")), Is.EqualTo(new[] { "Y.raw" }));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void Restrict_AShortRowIsPaddedWithEmptyCells()
    {
        string path = TempSdrf(
            new[] { "S1", "1", "run 1", "label free sample", "X.raw", "1" },
            new[] { "S2", "2", "run 2", "label free sample", "Y.raw" });
        try
        {
            JsonElement read = Invoke("sdrf", "read", "--path", path);
            Assert.That(read.GetProperty("rows")[1].GetArrayLength(), Is.EqualTo(5), "sdrf read keeps the row ragged");

            JsonElement data = InvokeWithStdin("Y.raw", "sdrf", "restrict", "--path", path);
            string[] row = Strings(data.GetProperty("rows")[0]);
            Assert.That(row, Has.Length.EqualTo(7), "six columns padded, plus comment[searched data file]");
            Assert.That(row[^1], Is.EqualTo(""));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void Restrict_WindowsThePayloadButWritesTheWholeDocumentIntoAFolderItCreates()
    {
        string root = Path.Combine(Path.GetTempPath(), $"restrict-out-{Guid.NewGuid():N}");
        string output = Path.Combine(root, "new", "searched.sdrf.tsv");
        try
        {
            JsonElement data = InvokeWithStdin(Clean, "sdrf", "restrict", "--path", Fixture(Pxd008841),
                "--limit", "2", "--offset", "1", "--out", output);

            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("offset").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
            JsonElement written = data.GetProperty("written");
            Assert.That(written.GetProperty("path").GetString(), Is.EqualTo(output));
            Assert.That(written.GetProperty("row_count").GetInt32(), Is.EqualTo(40));

            var back = new SdrfDocument(output);
            Assert.That(back.Results, Has.Count.EqualTo(40), "the whole document, not the window");
            Assert.That(back.Header.ToArray(), Does.Contain("comment[searched data file]"));
            Assert.That(File.ReadAllText(output), Does.Contain("\r\n"), "an in-memory document is written CRLF");
        }
        finally
        {
            if (Directory.Exists(root))
                Directory.Delete(root, recursive: true);
        }
    }

    [Test]
    public void Restrict_WrittenPathIsTheFileMzLibWrote()
    {
        // mzLib's WriteResults adds .sdrf.tsv to a name that does not end in it.
        string root = Path.Combine(Path.GetTempPath(), $"restrict-name-{Guid.NewGuid():N}");
        string output = Path.Combine(root, "searched.tsv");
        try
        {
            JsonElement data = InvokeWithStdin(Clean, "sdrf", "restrict", "--path", Fixture(Pxd008841), "--out", output);

            string written = data.GetProperty("written").GetProperty("path").GetString()!;
            Assert.That(written, Is.EqualTo(output + ".sdrf.tsv"));
            Assert.That(File.Exists(written), Is.True);
            Assert.That(File.Exists(output), Is.False);
        }
        finally
        {
            if (Directory.Exists(root))
                Directory.Delete(root, recursive: true);
        }
    }

    [TestCase("searched.txt")]
    [TestCase("searched")]
    [TestCase("searched.sdrf.csv")]
    public void Restrict_OutMustBeTsvAndIsCheckedBeforeAnyRead(string output)
    {
        // The SDRF does not exist and stdin is empty: the --out check must come first (PYB-4).
        JsonElement error = InvokeExpectingError("sdrf", "restrict", "--path", "absent.sdrf.tsv", "--out", output);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(".tsv"));
    }

    [Test]
    public void Restrict_TheCommittedSearchedDocumentIsWhatTheVerbWrites()
    {
        string output = Path.Combine(Path.GetTempPath(), $"restrict-{Guid.NewGuid():N}.sdrf.tsv");
        try
        {
            InvokeWithStdin(Clean, "sdrf", "restrict", "--path", Fixture(Pxd008841), "--out", output);

            // Compared cell for cell, not byte for byte: git may check the committed copy out with
            // either line ending.
            var fresh = new SdrfDocument(output);
            var committed = new SdrfDocument(Fixture("PXD008841_searched.sdrf.tsv"));
            Assert.That(committed.Header.ToArray(), Is.EqualTo(fresh.Header.ToArray()));
            Assert.That(committed.Results.Select(r => r.Cells.ToArray()).ToArray(),
                Is.EqualTo(fresh.Results.Select(r => r.Cells.ToArray()).ToArray()));
        }
        finally
        {
            File.Delete(output);
        }
    }

    [Test]
    public void Restrict_ThenIsobaricDesignKeysOnTheSearchedName()
    {
        // The design readers match searched files by exact name; the restricted document names them
        // in comment[searched data file], which both readers key on first (MAP-12).
        string output = Path.Combine(Path.GetTempPath(), $"restrict-design-{Guid.NewGuid():N}.sdrf.tsv");
        try
        {
            InvokeWithStdin(NativeClean, "sdrf", "restrict", "--path", Fixture(Pxd008841), "--out", output);
            JsonElement design = InvokeWithStdin(NativeClean, "sdrf", "isobaric-design", "--path", output,
                "--kit", "TMT10", "--plex-file-name-pattern", @"TMT_?pool(?<plex>\d+)", "--searched-files-stdin");

            Assert.That(design.GetProperty("is_valid").GetBoolean(), Is.True, design.GetProperty("report").GetString());
            Assert.That(design.GetProperty("file_key_column").GetString(), Is.EqualTo("comment[searched data file]"));
            Assert.That(Strings(design.GetProperty("columns").GetProperty("full_path")).Distinct(), Is.EqualTo(NativeCleanFiles));

            // Without the restriction, the same search matches no row: the SDRF names X.raw, the
            // search read X-calib.mzML.
            JsonElement unrestricted = InvokeWithStdin(NativeClean, "sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
                "--kit", "TMT10", "--plex-file-name-pattern", @"TMT_?pool(?<plex>\d+)", "--searched-files-stdin");
            Assert.That(unrestricted.GetProperty("is_valid").GetBoolean(), Is.False);
            Assert.That(Strings(unrestricted.GetProperty("refusals")), Has.Some.Contain("has no SDRF row"));
        }
        finally
        {
            File.Delete(output);
        }
    }

    [Test]
    public void Restrict_ADocumentWithNoDataFileColumnIsMzLibsFailureNotAUsageError()
    {
        string path = Path.Combine(Path.GetTempPath(), $"no-data-file-{Guid.NewGuid():N}.sdrf.tsv");
        var header = new SdrfHeader(new[] { "source name", "assay name" });
        new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "S1", "run 1" }) }).WriteResults(path);
        try
        {
            JsonElement error = InvokeExpectingErrorWithStdin("X.raw", "sdrf", "restrict", "--path", path);
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("ArgumentException"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("comment[data file]"));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void Restrict_UsageErrors()
    {
        string path = Fixture(Pxd008841);
        Assert.Multiple(() =>
        {
            Assert.That(Type(InvokeExpectingErrorWithStdin("", "sdrf", "restrict", "--path", path)), Is.EqualTo("usage"));
            Assert.That(InvokeExpectingErrorWithStdin("\n  \n", "sdrf", "restrict", "--path", path)
                .GetProperty("message").GetString(), Does.Contain("No searched files"));
            Assert.That(Type(InvokeExpectingErrorWithStdin("a.mzML\tb.raw\tc", "sdrf", "restrict", "--path", path)),
                Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingErrorWithStdin("\tb.raw", "sdrf", "restrict", "--path", path)),
                Is.EqualTo("usage"));
            Assert.That(InvokeExpectingErrorWithStdin("a.mzML\tb.raw\na.mzML\tc.raw", "sdrf", "restrict", "--path", path)
                .GetProperty("message").GetString(), Does.Contain("two acquired files"));
            Assert.That(InvokeExpectingErrorWithStdin("a.mzML", "sdrf", "restrict", "--path", "absent.sdrf.tsv")
                .GetProperty("message").GetString(), Does.Contain("not found"));
            Assert.That(Type(InvokeExpectingErrorWithStdin("a.mzML", "sdrf", "restrict", "--path", path, "--out")),
                Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingErrorWithStdin("a.mzML", "sdrf", "restrict", "--path", path, "--limit", "-1")),
                Is.EqualTo("usage"));
        });
    }

    [Test]
    public void Restrict_TheSameAcquiredNameTwiceIsNotAConflict()
    {
        JsonElement data = InvokeWithStdin($"sample_A.mzML\t{Pool5}1.raw\nsample_A.mzML\t{Pool5}1.raw",
            "sdrf", "restrict", "--path", Fixture(Pxd008841));

        Assert.That(data.GetProperty("acquired_name_count").GetInt32(), Is.EqualTo(1));
        Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(10));
    }

    [Test]
    public void Restrict_CarriesItsCaveats()
    {
        string[] caveats = Strings(InvokeWithStdin(Clean, "sdrf", "restrict", "--path", Fixture(Pxd008841))
            .GetProperty("caveats"));

        Assert.That(caveats, Has.Length.EqualTo(6));
        Assert.That(caveats[0], Does.Contain("STEM"));
        Assert.That(caveats, Has.Some.Contain(".sdrf.tsv"));
    }

    // ---- helpers ------------------------------------------------------------------------------

    /// <summary>A small SDRF written by mzLib to a temporary file; rows may be shorter than the header.</summary>
    private static string TempSdrf(params string[][] rows)
    {
        string[] columns =
            { "source name", "characteristics[biological replicate]", "assay name", "comment[label]", "comment[data file]", "comment[technical replicate]" };
        var header = new SdrfHeader(columns);
        string path = Path.Combine(Path.GetTempPath(), $"restrict-{Guid.NewGuid():N}.sdrf.tsv");
        new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r))).WriteResults(path);
        return path;
    }

    private static string Fixture(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "pkg", "python", "tests", "fixtures", name);
    }

    private static string Type(JsonElement error) => error.GetProperty("type").GetString()!;

    private static string[] Strings(JsonElement array) =>
        array.EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(null, args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) => Unwrap(Envelope(null, args), expectOk: false);

    private static JsonElement InvokeWithStdin(string stdin, params string[] args) => Unwrap(Envelope(stdin, args), expectOk: true);

    private static JsonElement InvokeExpectingErrorWithStdin(string stdin, params string[] args) =>
        Unwrap(Envelope(stdin, args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk), $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    private static JsonElement Envelope(string? stdin, string[] args)
    {
        TextReader previousIn = Console.In;
        Console.SetIn(new StringReader(stdin ?? ""));
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.SerializeToElement(new { ok = true, data }, Program.JsonOptions);
        }
        catch (Program.UsageException usage)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = "usage", message = usage.Message } }, Program.JsonOptions);
        }
        catch (Exception exception)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = Program.ClassifyError(exception), message = exception.Message } },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetIn(previousIn);
        }
    }
}
