using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Omics.Modifications;
using Readers;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>sdrf isobaric-design</c>, the projection of mzLib's <c>SdrfIsobaricDesign</c>.
/// </summary>
/// <remarks>
/// The documents are mzLib's own isobaric design fixtures
/// (<c>Test/FileReadingTests/ExternalFileTypes/SdrfDesign</c>, cut from bigbio/sdrf-annotated-datasets,
/// Apache-2.0, at 4f823dc), copied byte for byte into pyMzLib's fixture directory so the Python tests
/// share them. <b>PXD008841.subset</b>: 40 rows, 4 files, TMT10, two plexes whose name survives only in
/// a file-name token spelled two ways (<c>TMTpool1</c>, <c>TMT_pool5</c>). <b>PXD061609.subset</b>: 108
/// rows, 6 files, TMTpro 18, three SAX arms of ONE plex in <c>comment[sample preparation batch]</c>,
/// each arm numbering its fractions from 1. What is asserted is the projection, and enough of mzLib's
/// answer to prove the right call was made; the rules themselves are tested in mzLib.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SdrfIsobaricDesignTests
{
    private const string Pxd008841 = "PXD008841.subset.sdrf.tsv";
    private const string Pxd061609 = "PXD061609.subset.sdrf.tsv";
    private const string PoolPattern = @"TMT_?pool(?<plex>\d+)";
    private const string Batch = "comment[sample preparation batch]";

    private static readonly string[] Columns =
    {
        "full_path", "plex", "plex_id", "sample_name", "channel", "condition",
        "biological_replicate", "fraction", "technical_replicate",
    };

    [Test]
    public void IsobaricDesign_ReadsPxd008841ByItsFileNamePattern()
    {
        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
            "--kit", "TMT10", "--plex-file-name-pattern", PoolPattern);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True, data.GetProperty("report").GetString());
            Assert.That(data.GetProperty("kit").GetString(), Is.EqualTo("TMT10"));
            Assert.That(data.GetProperty("kit_type").GetString(), Is.EqualTo("TMT10"));
            Assert.That(Strings(data.GetProperty("kit_channels")), Has.Length.EqualTo(10));
            Assert.That(data.GetProperty("plex_file_name_pattern").GetString(), Is.EqualTo(PoolPattern));
            Assert.That(data.GetProperty("plex_column").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("single_plex").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("plex_source").GetString(), Does.Contain(PoolPattern).And.Contain("declared"));
            Assert.That(data.GetProperty("plex_ids").GetRawText(), Is.EqualTo("{\"1\":1,\"5\":2}"));
            Assert.That(data.GetProperty("plex_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("file_key_column").GetString(), Is.EqualTo("comment[data file]"));
            Assert.That(Strings(data.GetProperty("condition_columns")), Is.EqualTo(new[] { "factor value[disease]" }),
                "the document's only factor column, chosen by mzLib");
            Assert.That(data.GetProperty("condition_columns_declared").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("searched_files_given").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(4));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(40));
            Assert.That(data.GetProperty("refusals").GetArrayLength(), Is.Zero);
            Assert.That(data.GetProperty("notes").GetArrayLength(), Is.Zero);
            Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("report").GetString(), Does.Contain("4 file(s), 2 plex(es), kit TMT10"));
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(Columns));
        });

        JsonElement columns = data.GetProperty("columns");
        string[] channels = Strings(columns.GetProperty("channel"));
        Assert.That(channels.Take(10), Is.EqualTo(Strings(data.GetProperty("kit_channels"))),
            "a file's channels in reporter m/z order, as the kit lists them");
        Assert.That(Strings(columns.GetProperty("sample_name"))[9], Is.EqualTo("pool"), "131N is the bridge channel");
        Assert.That(Ints(columns.GetProperty("fraction")).Distinct().OrderBy(x => x), Is.EqualTo(new[] { 1, 2 }),
            "1-based, as the SDRF writes it");
        Assert.That(Ints(columns.GetProperty("technical_replicate")).Distinct(), Is.EqualTo(new[] { 1 }));
        Assert.That(Ints(columns.GetProperty("biological_replicate")).Distinct(), Is.EqualTo(new[] { 1 }));
        Assert.That(Ints(columns.GetProperty("plex_id")).Distinct(), Is.EqualTo(new[] { 1, 2 }));
        Assert.That(Strings(columns.GetProperty("full_path"))[0],
            Is.EqualTo("HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01.raw"), "the SDRF's file name");
    }

    [Test]
    public void IsobaricDesign_TheTableIsMzLibsDesignCellForCell()
    {
        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd061609),
            "--kit", "TMT18", "--plex-column", Batch);

        Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT18, out IsobaricMassTag? tag), Is.True);
        SdrfIsobaricDesign mzLib = SdrfIsobaricDesign.Read(Fixture(Pxd061609),
            new SdrfIsobaricDesignOptions { Tag = tag, PlexColumn = Batch });
        Assert.That(mzLib.IsValid, Is.True, mzLib.Report());

        var expected = mzLib.Files
            .SelectMany(f => f.Channels.Select(c => string.Join("|", f.FilePath, f.Plex, mzLib.PlexIds[f.Plex],
                c.SampleName, c.Label, c.Condition, c.BiologicalReplicate, f.Fraction, f.TechnicalReplicate)))
            .ToArray();
        JsonElement columns = data.GetProperty("columns");
        int rows = data.GetProperty("record_count").GetInt32();
        var actual = Enumerable.Range(0, rows)
            .Select(i => string.Join("|", Columns.Select(name => Cell(columns.GetProperty(name)[i]))))
            .ToArray();

        Assert.That(actual, Is.EqualTo(expected));
        Assert.That(rows, Is.EqualTo(108), "6 files x 18 channels");
        Assert.That(data.GetProperty("plex_ids").GetRawText(), Is.EqualTo("{\"HighSalt\":1,\"LowSalt\":2,\"NoSAX\":3}"),
            "1..N in case-insensitive name order, as MetaMorpheus numbers them");
        Assert.That(Strings(data.GetProperty("notes")), Is.EqualTo(mzLib.Notes));
        Assert.That(data.GetProperty("report").GetString(), Is.EqualTo(mzLib.Report()));
    }

    [Test]
    public void IsobaricDesign_OneSaxPlexIsRefusedEveryReasonAtOnce()
    {
        // The trap MAP-09 names: the batch column holds three arms of one plex, and read as ONE plex
        // every arm's fraction 1 collides with the others'.
        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd061609),
            "--kit", "TMT18", "--single-plex", "tube1");

        string[] refusals = Strings(data.GetProperty("refusals"));
        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(refusals, Has.Length.EqualTo(36));
        Assert.That(refusals, Has.All.Contain("more than once"));
        Assert.That(data.GetProperty("single_plex").GetString(), Is.EqualTo("tube1"));
        Assert.That(data.GetProperty("file_count").GetInt32(), Is.Zero);
        Assert.That(data.GetProperty("record_count").GetInt32(), Is.Zero);
        Assert.That(Strings(data.GetProperty("columns").GetProperty("full_path")), Is.Empty,
            "a refused design carries no rows, so nothing can be used by mistake");
        Assert.That(data.GetProperty("plex_ids").GetRawText(), Is.EqualTo("{}"));
    }

    [Test]
    public void IsobaricDesign_NoPlexSourceIsMzLibsRefusalAndWritesNothing()
    {
        string folder = Path.Combine(Path.GetTempPath(), $"isobaric-refused-{Guid.NewGuid():N}");

        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
            "--out", folder);

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(Strings(data.GetProperty("refusals")).Single(), Does.Contain("No plex source was declared"));
        Assert.That(data.GetProperty("plex_source").GetString(), Is.EqualTo("(not declared)"));
        Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Directory.Exists(folder), Is.False, "a refusal must never leave a design for MetaMorpheus to find");
    }

    [Test]
    public void IsobaricDesign_NoKitIsMzLibsRefusalNeverAGuess()
    {
        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
            "--plex-file-name-pattern", PoolPattern);

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(Strings(data.GetProperty("refusals")).Single(), Does.Contain("No isobaric tag was given"));
        Assert.That(data.GetProperty("kit").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(data.GetProperty("kit_type").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(data.GetProperty("kit_channels").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void IsobaricDesign_ResolvesTheKitExactlyAsIsobaricKitsDoes()
    {
        foreach (string name in new[] { "tmt10", "TMT10 on K" })
        {
            JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
                "--kit", name, "--plex-file-name-pattern", PoolPattern);
            Assert.That(data.GetProperty("kit").GetString(), Is.EqualTo(name), "echoed as given");
            Assert.That(data.GetProperty("kit_type").GetString(), Is.EqualTo("TMT10"));
        }

        string design = InvokeExpectingError("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
            "--kit", "TMT10plex").GetProperty("message").GetString()!;
        string kits = InvokeExpectingError("isobaric", "kits", "--kit", "TMT10plex").GetProperty("message").GetString()!;
        Assert.That(design, Is.EqualTo(kits), "one resolution path, one message");
    }

    [Test]
    public void IsobaricDesign_OnlyAnnotatedChannelsAreRows()
    {
        // TMT11 has 131C; PXD008841 annotates ten channels. The kit says so, the table does not
        // invent a row for the eleventh.
        JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
            "--kit", "TMT11", "--plex-file-name-pattern", PoolPattern);

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True);
        Assert.That(Strings(data.GetProperty("kit_channels")), Has.Length.EqualTo(11).And.Contain("131C"));
        Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(40));
        Assert.That(Strings(data.GetProperty("columns").GetProperty("channel")), Has.None.EqualTo("131C"));
    }

    [Test]
    public void IsobaricDesign_TwoPlexSourcesAndABadPatternAreRefusals()
    {
        JsonElement two = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
            "--plex-file-name-pattern", PoolPattern, "--single-plex", "P1");
        Assert.That(Strings(two.GetProperty("refusals")).Single(), Does.Contain("2 plex sources were declared"));

        JsonElement bad = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
            "--plex-file-name-pattern", "pool(");
        Assert.That(Strings(bad.GetProperty("refusals")), Has.Some.Contain("not a valid regular expression"));

        JsonElement column = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
            "--plex-column", "comment[plex]");
        Assert.That(Strings(column.GetProperty("refusals")), Has.Some.Contain("'comment[plex]' is not in the SDRF"));
    }

    [Test]
    public void IsobaricDesign_SearchedFilesDropOtherRowsAndSupplyThePaths()
    {
        string folder = Path.Combine(Path.GetTempPath(), "isobaric-searched");
        string[] searched =
        {
            Path.Combine(folder, "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr01.raw"),
            Path.Combine(folder, "HJOS2U_20140410_TMTpool1_300ugIPG37-49_7of15ul_fr02.raw"),
        };

        JsonElement data = InvokeWithStdin(string.Join("\n", searched),
            "sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
            "--plex-file-name-pattern", PoolPattern, "--searched-files-stdin");

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True, data.GetProperty("report").GetString());
        Assert.That(data.GetProperty("searched_files_given").GetBoolean(), Is.True);
        Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(2));
        Assert.That(Strings(data.GetProperty("columns").GetProperty("full_path")).Distinct(), Is.EqualTo(searched),
            "the searched paths, not the SDRF's bare names");
        Assert.That(data.GetProperty("plex_ids").GetRawText(), Is.EqualTo("{\"1\":1}"));
        string[] notes = Strings(data.GetProperty("notes"));
        Assert.That(notes, Has.Length.EqualTo(20), "every row of the two unsearched files");
        Assert.That(notes, Has.All.Contain("dropped: the search does not read that file"));
    }

    [Test]
    public void IsobaricDesign_OutWritesTmtDesignTxtInsideTheFolderItCreates()
    {
        string root = Path.Combine(Path.GetTempPath(), $"isobaric-out-{Guid.NewGuid():N}");
        string folder = Path.Combine(root, "design");
        try
        {
            JsonElement data = Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841), "--kit", "TMT10",
                "--plex-file-name-pattern", PoolPattern, "--out", folder);

            string expectedPath = Path.Combine(folder, "TmtDesign.txt");
            JsonElement written = data.GetProperty("written");
            Assert.That(written.GetProperty("path").GetString(), Is.EqualTo(expectedPath));
            Assert.That(written.GetProperty("file_count").GetInt32(), Is.EqualTo(4));
            Assert.That(written.GetProperty("record_count").GetInt32(), Is.EqualTo(40));

            string[] lines = File.ReadAllLines(expectedPath);
            Assert.That(lines[0], Is.EqualTo(SdrfIsobaricDesign.TmtDesignHeader), "mzLib's header, MetaMorpheus's format");
            Assert.That(lines, Has.Length.EqualTo(41));
            JsonElement columns = data.GetProperty("columns");
            // TmtDesign.txt: File, Plex, Sample Name, TMT Channel, Condition, Biological Replicate,
            // Fraction, Technical Replicate; the wire adds plex_id after plex and converts nothing.
            string[] first = lines[1].Split('\t');
            string[] wire = Columns.Where(c => c != "plex_id").Select(c => Cell(columns.GetProperty(c)[0])).ToArray();
            Assert.That(first, Is.EqualTo(wire), "the file and the wire hold the same 1-based numbers");
        }
        finally
        {
            if (Directory.Exists(root))
                Directory.Delete(root, recursive: true);
        }
    }

    [Test]
    public void IsobaricDesign_OutNamingAFileIsAUsageErrorBeforeAnyRead()
    {
        string file = Path.Combine(Path.GetTempPath(), $"TmtDesign-{Guid.NewGuid():N}.txt");
        File.WriteAllText(file, "not a folder");
        try
        {
            JsonElement error = InvokeExpectingError("sdrf", "isobaric-design", "--path", "absent.sdrf.tsv",
                "--kit", "TMT10", "--out", file);
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("names a folder"));
        }
        finally
        {
            File.Delete(file);
        }
    }

    [Test]
    public void IsobaricDesign_UsageErrors()
    {
        string path = Fixture(Pxd008841);
        Assert.Multiple(() =>
        {
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design")), Is.EqualTo("usage"));
            Assert.That(InvokeExpectingError("sdrf", "isobaric-design", "--path", "absent.sdrf.tsv")
                .GetProperty("message").GetString(), Does.Contain("not found"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--kit")), Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--kit", "TMT10",
                "--plex-column")), Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--kit", "TMT10",
                "--plex-file-name-pattern", " ")), Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--kit", "TMT10",
                "--single-plex")), Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--kit", "TMT10",
                "--condition-columns", "factor value[disease]\t ")), Is.EqualTo("usage"));
            Assert.That(Type(InvokeExpectingError("sdrf", "isobaric-design", "--path", path, "--out")), Is.EqualTo("usage"));
        });
    }

    [Test]
    public void IsobaricDesign_CarriesItsCaveats()
    {
        string[] caveats = Strings(Invoke("sdrf", "isobaric-design", "--path", Fixture(Pxd008841),
            "--kit", "TMT10", "--plex-file-name-pattern", PoolPattern).GetProperty("caveats"));

        Assert.That(caveats, Has.Length.EqualTo(8));
        Assert.That(caveats[0], Does.StartWith("1-BASED").And.Contain("sdrf design is 0-based"));
        Assert.That(caveats, Has.Some.Contain("TmtDesign.txt"));
    }

    // ---- helpers ------------------------------------------------------------------------------

    private static string Fixture(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "pkg", "python", "tests", "fixtures", name);
    }

    private static string Cell(JsonElement value) =>
        value.ValueKind == JsonValueKind.Number ? value.GetInt32().ToString() : value.GetString()!;

    private static string Type(JsonElement error) => error.GetProperty("type").GetString()!;

    private static string[] Strings(JsonElement array) =>
        array.EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static int[] Ints(JsonElement array) =>
        array.EnumerateArray().Select(x => x.GetInt32()).ToArray();

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(null, args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) => Unwrap(Envelope(null, args), expectOk: false);

    private static JsonElement InvokeWithStdin(string stdin, params string[] args) => Unwrap(Envelope(stdin, args), expectOk: true);

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
