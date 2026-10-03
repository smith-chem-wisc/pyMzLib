using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>sdrf design</c>, the projection of mzLib's <c>SdrfLabelFreeDesign</c>.
/// </summary>
/// <remarks>
/// The two real documents are mzLib's own design fixtures (aging's hand-written SDRFs), copied into
/// pyMzLib's fixture directory so the Python tests share them. <b>PXD067622</b> is a valid
/// 24-file, two-factor design; <b>PXD049018</b> writes its treatment as <c>not available</c>, which
/// mzLib refuses. <b>PXD067622_studywide</b> is PXD067622 with each biological replicate replaced by
/// the study-wide index in its file name (WT_DMSO1..CA_FA24), which is how a drafted SDRF numbers
/// them and what mzLib ranks back to 1..3 within each condition. What is asserted is the projection
/// and enough of mzLib's answer to prove the right call was made; the rules are tested in mzLib.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SdrfDesignTests
{
    private const string Genotype = "factor value[genotype]";
    private const string Treatment = "factor value[treatment]";
    private static readonly string BothFactors = $"{Genotype}\t{Treatment}";

    [Test]
    public void Design_ReadsAValidTwoFactorDesign()
    {
        JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
            "--condition-columns", BothFactors);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("file_key_column").GetString(), Is.EqualTo("comment[data file]"));
            Assert.That(Strings(data.GetProperty("condition_columns")), Is.EqualTo(new[] { Genotype, Treatment }));
            Assert.That(data.GetProperty("condition_columns_declared").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("searched_files_given").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(24));
            Assert.That(data.GetProperty("refusals").GetArrayLength(), Is.Zero);
            Assert.That(data.GetProperty("notes").GetArrayLength(), Is.Zero, "the SDRF already numbers 1..3 per condition");
            Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("report").GetString(), Does.Contain("24 file(s), 8 condition(s)"));
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(new[]
                { "full_path", "file_name", "condition", "biological_replicate", "technical_replicate", "fraction" }));
        });

        JsonElement columns = data.GetProperty("columns");
        string[] conditions = Strings(columns.GetProperty("condition"));
        Assert.That(conditions.Distinct().Count(), Is.EqualTo(8));
        Assert.That(conditions, Does.Contain("SPRTN-TurboID WT_DMSO (vehicle)"), "the declared columns joined with '_'");
        Assert.That(Ints(columns.GetProperty("biological_replicate")).Distinct().OrderBy(x => x), Is.EqualTo(new[] { 0, 1, 2 }),
            "0-based, as mzLib's SpectraFileInfo holds them");
        Assert.That(Ints(columns.GetProperty("fraction")).Distinct(), Is.EqualTo(new[] { 0 }));
        string[] fileNames = Strings(columns.GetProperty("file_name"));
        string[] paths = Strings(columns.GetProperty("full_path"));
        Assert.That(fileNames.Zip(paths).All(p => p.Second == p.First + ".raw"), Is.True,
            "file_name is the path without its extension, the key the quant verbs use");
    }

    [Test]
    public void Design_AnUndeclaredSecondFactorIsARefusalNotAnError()
    {
        JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"));

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(data.GetProperty("condition_columns_declared").GetBoolean(), Is.False);
        Assert.That(data.GetProperty("file_count").GetInt32(), Is.Zero);
        Assert.That(Strings(data.GetProperty("columns").GetProperty("full_path")), Is.Empty,
            "a refused design carries no rows, so nothing can be used by mistake");
        string refusal = Strings(data.GetProperty("refusals")).Single();
        Assert.That(refusal, Does.Contain("Several factor value columns and none declared")
            .And.Contain(Genotype).And.Contain(Treatment));
        Assert.That(data.GetProperty("report").GetString(), Does.Contain("REFUSED"));
    }

    [Test]
    public void Design_ListsEveryRefusalAtOnce()
    {
        JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD049018.sdrf.tsv"),
            "--condition-columns", BothFactors);

        string[] refusals = Strings(data.GetProperty("refusals"));
        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(refusals, Has.Length.EqualTo(20), "one per row, each naming its file");
        Assert.That(refusals[0], Does.Contain("MSB67868ABand_01.raw").And.Contain($"'{Treatment}' is 'not available'"));
    }

    [Test]
    public void Design_ReportsTheReplicateRenumberingInNotes()
    {
        JsonElement studyWide = Invoke("sdrf", "design", "--path", Fixture("PXD067622_studywide.sdrf.tsv"),
            "--condition-columns", BothFactors);
        JsonElement reference = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
            "--condition-columns", BothFactors);

        string[] notes = Strings(studyWide.GetProperty("notes"));
        Assert.That(studyWide.GetProperty("is_valid").GetBoolean(), Is.True);
        Assert.That(notes, Has.Length.EqualTo(7), "every condition but WT_DMSO, which is already 1..3");
        Assert.That(notes, Does.Contain(
            "Condition 'SPRTN-TurboID CA_formaldehyde 1 mM, 1 h': biological replicates renumbered 22 -> 1, 23 -> 2, 24 -> 3."));
        Assert.That(studyWide.GetProperty("columns").GetRawText(), Is.EqualTo(reference.GetProperty("columns").GetRawText()),
            "ranked within each condition, the design is the hand-numbered one");
    }

    [Test]
    public void Design_SearchedFilesDropUnreadRowsAndSupplyThePaths()
    {
        string directory = Path.Combine(Path.GetTempPath(), "sdrf-design-searched");
        string[] searched =
        {
            Path.Combine(directory, "20240830_HF_LC3_MAA_RK_12032_WT_DMSO1.raw"),
            Path.Combine(directory, "20240830_HF_LC3_MAA_RK_12032_WT_DMSO2.raw"),
            Path.Combine(directory, "20240830_HF_LC3_MAA_RK_12032_WT_DMSO3.raw"),
        };

        JsonElement data = InvokeWithStdin(string.Join("\n", searched),
            "sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"), "--condition-columns", BothFactors,
            "--searched-files-stdin");

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True, data.GetProperty("report").GetString());
        Assert.That(data.GetProperty("searched_files_given").GetBoolean(), Is.True);
        Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(3));
        Assert.That(Strings(data.GetProperty("columns").GetProperty("full_path")), Is.EqualTo(searched),
            "the searched paths, not the SDRF's bare names");
        string[] notes = Strings(data.GetProperty("notes"));
        Assert.That(notes, Has.Length.EqualTo(21), "every row the search does not read is reported");
        Assert.That(notes[0], Does.Contain("dropped: the search does not read that file"));
    }

    [Test]
    public void Design_OutWritesMetaMorpheusOneBasedDesign()
    {
        string output = Path.Combine(Path.GetTempPath(), $"design-{Guid.NewGuid():N}", "ExperimentalDesign.tsv");
        Directory.CreateDirectory(Path.GetDirectoryName(output)!);
        try
        {
            JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
                "--condition-columns", BothFactors, "--out", output);

            JsonElement written = data.GetProperty("written");
            Assert.That(written.GetProperty("path").GetString(), Is.EqualTo(output));
            Assert.That(written.GetProperty("file_count").GetInt32(), Is.EqualTo(24));

            string[] lines = File.ReadAllLines(output);
            Assert.That(lines[0], Is.EqualTo("FileName\tCondition\tBiorep\tFraction\tTechrep"));
            Assert.That(lines, Has.Length.EqualTo(25));
            string[] first = lines[1].Split('\t');
            int wireBiorep = data.GetProperty("columns").GetProperty("biological_replicate")[0].GetInt32();
            Assert.That(first[0], Is.EqualTo(Strings(data.GetProperty("columns").GetProperty("full_path"))[0]));
            Assert.That(int.Parse(first[2]), Is.EqualTo(wireBiorep + 1), "the file is 1-based, the wire 0-based");
        }
        finally
        {
            Directory.Delete(Path.GetDirectoryName(output)!, recursive: true);
        }
    }

    [Test]
    public void Design_ARefusedDesignWritesNoFile()
    {
        string output = Path.Combine(Path.GetTempPath(), $"refused-{Guid.NewGuid():N}.tsv");

        JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD049018.sdrf.tsv"),
            "--condition-columns", BothFactors, "--out", output);

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(File.Exists(output), Is.False, "a refusal must never leave a design for MetaMorpheus to find");
    }

    [TestCase("ExperimentalDesign.txt")]
    [TestCase("ExperimentalDesign")]
    [TestCase("design.csv")]
    public void Design_OutMustBeTsvAndIsCheckedBeforeAnyRead(string output)
    {
        // The SDRF does not exist: the --out check must come first (PYB-1).
        JsonElement error = InvokeExpectingError("sdrf", "design", "--path", "absent.sdrf.tsv", "--out", output);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(".tsv"));
    }

    [Test]
    public void Design_OutUpperCaseTsvIsAccepted()
    {
        string output = Path.Combine(Path.GetTempPath(), $"design-{Guid.NewGuid():N}.TSV");
        try
        {
            JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
                "--condition-columns", BothFactors, "--out", output);
            Assert.That(data.GetProperty("written").GetProperty("path").GetString(), Is.EqualTo(output));
        }
        finally
        {
            File.Delete(output);
        }
    }

    [Test]
    public void Design_UsageErrors()
    {
        Assert.Multiple(() =>
        {
            Assert.That(InvokeExpectingError("sdrf", "design").GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(InvokeExpectingError("sdrf", "design", "--path", "absent.sdrf.tsv")
                .GetProperty("message").GetString(), Does.Contain("not found"));
            Assert.That(InvokeExpectingError("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
                "--condition-columns", $"{Genotype}\t ").GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(InvokeExpectingError("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
                "--condition-columns").GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(InvokeExpectingError("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
                "--out").GetProperty("type").GetString(), Is.EqualTo("usage"));
        });
    }

    [Test]
    public void Design_AnUnknownConditionColumnIsMzLibsRefusal()
    {
        JsonElement data = Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
            "--condition-columns", "factor value[dose]");

        Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
        Assert.That(Strings(data.GetProperty("refusals")), Has.Some.Contain("factor value[dose]"));
    }

    [Test]
    public void Design_CarriesItsCaveats()
    {
        string[] caveats = Strings(Invoke("sdrf", "design", "--path", Fixture("PXD067622.sdrf.tsv"),
            "--condition-columns", BothFactors).GetProperty("caveats"));

        Assert.That(caveats, Has.Length.EqualTo(6));
        Assert.That(caveats[0], Does.StartWith("LABEL-FREE ONLY"));
        Assert.That(caveats, Has.Some.Contain("0-based"));
    }

    // ---- helpers ------------------------------------------------------------------------------

    private static string Fixture(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "pkg", "python", "tests", "fixtures", name);
    }

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
