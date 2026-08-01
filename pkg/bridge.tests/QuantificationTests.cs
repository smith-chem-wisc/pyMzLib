using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using MassSpectrometry;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the FlashLFQ verb's argument boundary — the spectra-design parsing that turns stdin
/// lines into <see cref="SpectraFileInfo"/>s.
/// </summary>
/// <remarks>
/// The quantification itself is mzLib's and is covered by mzLib's own FlashLFQ tests. What lives
/// only in the bridge, and is where a bug would land, is this translation: the design defaults, the
/// mzML-only guard, and the up-front errors that turn a deep engine failure into one clear message.
/// These need no engine and no real spectra — only files that exist, so 0-byte placeholders suffice
/// (existence is checked here; the mzML is parsed later, by the engine).
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class QuantificationTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-quant-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void Cleanup()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    /// <summary>Creates an empty file at <paramref name="relativePath"/> and returns its full path.</summary>
    private string TouchMzml(string relativePath)
    {
        string path = Path.Combine(_tempDirectory, relativePath);
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        File.WriteAllText(path, string.Empty);
        return path;
    }

    [Test]
    public void BarePaths_GetDefaultDesign_EachItsOwnBiologicalReplicate()
    {
        string a = TouchMzml("a.mzML");
        string b = TouchMzml("b.mzML");

        List<SpectraFileInfo> files = Quantification.BuildSpectraFiles(new[] { a, b });

        Assert.That(files, Has.Count.EqualTo(2));
        Assert.That(files[0].FilenameWithoutExtension, Is.EqualTo("a"));
        Assert.That(files[0].Condition, Is.EqualTo(string.Empty));
        Assert.That(files[0].BiologicalReplicate, Is.EqualTo(0));
        Assert.That(files[1].BiologicalReplicate, Is.EqualTo(1), "each bare path becomes its own biorep");
        Assert.That(files[0].TechnicalReplicate, Is.EqualTo(0));
        Assert.That(files[0].Fraction, Is.EqualTo(0));
    }

    [Test]
    public void DesignFields_AreParsedInOrder()
    {
        string a = TouchMzml("a.mzML");

        SpectraFileInfo file = Quantification.BuildSpectraFiles(new[] { $"{a}\tcontrol\t2\t1\t3" })[0];

        Assert.That(file.Condition, Is.EqualTo("control"));
        Assert.That(file.BiologicalReplicate, Is.EqualTo(2));
        Assert.That(file.TechnicalReplicate, Is.EqualTo(1));
        Assert.That(file.Fraction, Is.EqualTo(3));
    }

    [Test]
    public void EmptyMiddleField_FallsBackToDefault_ButLaterFieldStillReads()
    {
        string a = TouchMzml("a.mzML");

        // condition set, biorep and techrep blank, fraction = 5.
        SpectraFileInfo file = Quantification.BuildSpectraFiles(new[] { $"{a}\tcond\t\t\t5" })[0];

        Assert.That(file.Condition, Is.EqualTo("cond"));
        Assert.That(file.BiologicalReplicate, Is.EqualTo(0), "blank biorep falls back to the line index");
        Assert.That(file.TechnicalReplicate, Is.EqualTo(0));
        Assert.That(file.Fraction, Is.EqualTo(5));
    }

    [Test]
    public void NonMzml_IsRejected()
    {
        string raw = TouchMzml("a.raw");

        Assert.Throws<Program.UsageException>(() => Quantification.BuildSpectraFiles(new[] { raw }));
    }

    [Test]
    public void MissingFile_IsRejected()
    {
        string missing = Path.Combine(_tempDirectory, "not-here.mzML");

        Assert.Throws<Program.UsageException>(() => Quantification.BuildSpectraFiles(new[] { missing }));
    }

    [Test]
    public void DuplicateBaseName_IsRejected()
    {
        // Same base name, different directories: FlashLFQ matches ids to runs by base name, so this
        // would be ambiguous and must fail up front rather than silently quantify against one.
        string first = TouchMzml("x.mzML");
        string second = TouchMzml(Path.Combine("sub", "x.mzML"));

        Assert.Throws<Program.UsageException>(() => Quantification.BuildSpectraFiles(new[] { first, second }));
    }

    [Test]
    public void NoSpectra_IsRejected()
    {
        Assert.Throws<Program.UsageException>(() => Quantification.BuildSpectraFiles(Array.Empty<string>()));
    }

    [Test]
    public void NonIntegerDesignField_IsRejected()
    {
        string a = TouchMzml("a.mzML");

        Assert.Throws<Program.UsageException>(
            () => Quantification.BuildSpectraFiles(new[] { $"{a}\tcond\tnotanumber" }));
    }

    // -----------------------------------------------------------------------------------------
    // End-to-end: run the real verb against the mzLib FlashLFQ test data. This is what exercises
    // the engine call, the console redirect, the result projection (peptides/proteins/peaks), and
    // WriteResults — none of which the boundary tests above reach.
    // -----------------------------------------------------------------------------------------

    /// <summary>
    /// The FlashLFQ test data in the pinned mzLib worktree, located relative to this source file so
    /// it resolves both locally and in CI (where mzLib is cloned to the same relative path).
    /// </summary>
    private static string TestDataDirectory([CallerFilePath] string thisFile = "")
    {
        // thisFile = <root>/pkg/bridge.tests/QuantificationTests.cs → up to <root>, then into mzLib.
        string root = Directory.GetParent(thisFile)!.Parent!.Parent!.FullName;
        return Path.Combine(root, "code", "mzLib", "mzLib", "Test", "FlashLFQ", "TestData");
    }

    [Test]
    public void QuantifyRealData_RunsEndToEnd_ProducesPeptidesProteinsAndMbrPeaks()
    {
        string data = TestDataDirectory();
        string psms = Path.Combine(data, "AllPSMs.psmtsv");
        string run3 = Path.Combine(data, "20100614_Velos1_TaGe_SA_K562_3.mzML");
        string run4 = Path.Combine(data, "20100614_Velos1_TaGe_SA_K562_4.mzML");
        if (!File.Exists(psms) || !File.Exists(run3) || !File.Exists(run4))
            Assert.Ignore("FlashLFQ test data not present in the mzLib worktree.");

        string outDir = Path.Combine(_tempDirectory, "out");
        object result = RunFlashLfq(
            $"{run3}\n{run4}\n",
            "quant", "flashlfq", "--psms", psms, "--mbr", "--shared-peptides",
            "--ppm", "10", "--mbr-ppm", "10", "--mbr-q", "0.05", "--threads", "1", "--out", outDir);

        using JsonDocument doc = JsonDocument.Parse(JsonSerializer.Serialize(result, Program.JsonOptions));
        JsonElement root = doc.RootElement;

        // The set of quantified peptides/proteins comes from the identifications, so it is stable;
        // peak and MBR counts depend on the ML-based PEP model, so they are asserted loosely to stay
        // robust across platforms — exact numbers are pinned by the Python and real-run checks.
        Assert.That(root.GetProperty("identification_count").GetInt32(), Is.EqualTo(594));
        Assert.That(root.GetProperty("peptides").GetArrayLength(), Is.GreaterThan(300));
        Assert.That(root.GetProperty("proteins").GetArrayLength(), Is.GreaterThan(900));

        JsonElement peaks = root.GetProperty("peaks");
        Assert.That(peaks.GetArrayLength(), Is.GreaterThan(500));
        int mbrPeaks = 0;
        foreach (JsonElement peak in peaks.EnumerateArray())
        {
            if (peak.GetProperty("detection_type").GetString() == "MBR")
                mbrPeaks++;
        }
        Assert.That(mbrPeaks, Is.GreaterThan(50), "the peaks surface must carry the MBR transfers");

        Assert.That(File.Exists(Path.Combine(outDir, "QuantifiedPeptides.tsv")), Is.True);
        Assert.That(File.Exists(Path.Combine(outDir, "QuantifiedPeaks.tsv")), Is.True);
    }

    [Test]
    public void PsmReferencingAnUnprovidedRun_IsRejected()
    {
        string data = TestDataDirectory();
        string psms = Path.Combine(data, "AllPSMs.psmtsv");
        string run3 = Path.Combine(data, "20100614_Velos1_TaGe_SA_K562_3.mzML");
        if (!File.Exists(psms) || !File.Exists(run3))
            Assert.Ignore("FlashLFQ test data not present in the mzLib worktree.");

        // The PSM file names run_3 AND run_4; provide only run_3, so run_4 is referenced-but-missing.
        Assert.Throws<Program.UsageException>(
            () => RunFlashLfq($"{run3}\n", "quant", "flashlfq", "--psms", psms));
    }

    [Test]
    public void MissingPsmFile_IsRejected()
    {
        string mzml = TouchMzml("a.mzML");
        string missingPsm = Path.Combine(_tempDirectory, "not-here.psmtsv");

        Assert.Throws<Program.UsageException>(
            () => RunFlashLfq($"{mzml}\n", "quant", "flashlfq", "--psms", missingPsm));
    }

    /// <summary>Runs the verb with the given stdin, restoring the console afterwards.</summary>
    private static object RunFlashLfq(string stdin, params string[] args)
    {
        TextReader originalIn = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            return Quantification.FlashLfq(new Program.Arguments(args));
        }
        finally
        {
            Console.SetIn(originalIn);
        }
    }

    // -----------------------------------------------------------------------------------------
    // Median polish: roll a QuantifiedPeptides.tsv up to proteins without re-running peak-finding.
    // Unlike the FlashLFQ engine, this whole path — parse the table, rebuild the object graph, run
    // the (mzLib) median polish, project the result — lives in the bridge, so it is exercised here
    // end-to-end against small tables written into the temp directory. No mzML and no engine run.
    // -----------------------------------------------------------------------------------------

    /// <summary>The standard peptide-table columns before the per-run intensity/detection blocks.</summary>
    private const string PeptideTableLead = "Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism";

    /// <summary>
    /// Writes a QuantifiedPeptides.tsv for runs <c>run_1..run_N</c> and returns its path. Each row is
    /// (sequence, protein group, then one intensity per run); every intensity is an MS/MS detection.
    /// </summary>
    private string WritePeptideTable(string[] runs, params (string Seq, string Protein, double[] Intensities)[] rows)
    {
        var header = new System.Text.StringBuilder(PeptideTableLead);
        foreach (string run in runs) header.Append("\tIntensity_").Append(run);
        foreach (string run in runs) header.Append("\tDetection Type_").Append(run);

        var lines = new List<string> { header.ToString() };
        foreach (var (seq, protein, intensities) in rows)
        {
            var line = new System.Text.StringBuilder();
            line.Append(seq).Append('\t').Append(seq).Append('\t').Append(protein).Append("\tGENE\tHomo sapiens");
            foreach (double i in intensities)
                line.Append('\t').Append(i.ToString(System.Globalization.CultureInfo.InvariantCulture));
            foreach (double i in intensities)
                line.Append('\t').Append(i > 0 ? "MSMS" : "NotDetected");
            lines.Add(line.ToString());
        }

        string path = Path.Combine(_tempDirectory, "QuantifiedPeptides.tsv");
        File.WriteAllLines(path, lines);
        return path;
    }

    /// <summary>Runs the median-polish verb with the given stdin, returning its result as parsed JSON.</summary>
    private static JsonElement RunMedianPolish(string stdin, params string[] args)
    {
        TextReader originalIn = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            object result = Quantification.MedianPolish(new Program.Arguments(args));
            return JsonDocument.Parse(JsonSerializer.Serialize(result, Program.JsonOptions)).RootElement;
        }
        finally
        {
            Console.SetIn(originalIn);
        }
    }

    private static Dictionary<string, JsonElement> ProteinsByName(JsonElement root) =>
        root.GetProperty("proteins").EnumerateArray().ToDictionary(p => p.GetProperty("protein_group").GetString()!);

    private static double? Intensity(JsonElement protein, string sample)
    {
        JsonElement value = protein.GetProperty("intensities").GetProperty(sample);
        return value.ValueKind == JsonValueKind.Null ? null : value.GetDouble();
    }

    [Test]
    public void MedianPolish_NoDesign_LabelsEachRunAndTracksFoldChange()
    {
        // Two runs at one level, two at double: a protein's intensity should track that 2x step, and
        // with no design each run is its own sample keyed by its base name.
        string[] runs = { "run_1", "run_2", "run_3", "run_4" };
        string path = WritePeptideTable(runs,
            ("PEPTIDEK", "P1", new[] { 1000.0, 1100, 2000, 2200 }),
            ("AAAAAR", "P1", new[] { 500.0, 550, 1000, 1050 }),
            ("LLLLLK", "P1", new[] { 2000.0, 2100, 4000, 4100 }));

        JsonElement root = RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path);
        JsonElement p1 = ProteinsByName(root)["P1"];

        Assert.That(Intensity(p1, "run_1"), Is.GreaterThan(0));
        Assert.That(Intensity(p1, "run_3")! / Intensity(p1, "run_1")!, Is.EqualTo(2.0).Within(0.1),
            "run_3 is 2x run_1 across every peptide, so the protein intensity should double too");
    }

    [Test]
    public void MedianPolish_Design_GroupsByConditionAndBiorepLabels()
    {
        string[] runs = { "run_1", "run_2" };
        string path = WritePeptideTable(runs,
            ("PEPTIDEK", "P1", new[] { 1000.0, 2000 }),
            ("AAAAAR", "P1", new[] { 500.0, 1000 }));

        string design = "run_1\tcontrol\t0\nrun_2\ttreated\t0\n";
        JsonElement root = RunMedianPolish(design, "quant", "median-polish", "--peptides", path);

        // Real conditions give condition_(biorep+1) sample labels.
        JsonElement p1 = ProteinsByName(root)["P1"];
        Assert.That(p1.GetProperty("intensities").TryGetProperty("control_1", out _), Is.True);
        Assert.That(p1.GetProperty("intensities").TryGetProperty("treated_1", out _), Is.True);
    }

    [Test]
    public void MedianPolish_SharedPeptides_LiftsAProteinWithOnlySharedPeptides()
    {
        // Four runs with a signal that varies across them, so the shared peptide resolves to a real
        // number rather than the NaN a degenerate (identical-across-runs) matrix would give. P2 has
        // its own unique peptides; P3's only peptide is the one it shares with P2.
        string[] runs = { "run_1", "run_2", "run_3", "run_4" };
        string path = WritePeptideTable(runs,
            ("CCCCCR", "P2", new[] { 800.0, 900, 850, 820 }),
            ("DDDDDK", "P2", new[] { 1600.0, 1700, 1650, 1620 }),
            ("EEEEEK", "P2;P3", new[] { 400.0, 410, 420, 430 }));

        // Without shared peptides P3 has no unique peptide and cannot be quantified: it stays 0.
        JsonElement without = RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path);
        Assert.That(Intensity(ProteinsByName(without)["P3"], "run_1"), Is.EqualTo(0.0));

        // With shared peptides the shared peptide quantifies P3 to a real intensity.
        JsonElement with = RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path, "--shared-peptides");
        Assert.That(Intensity(ProteinsByName(with)["P3"], "run_1"), Is.GreaterThan(0));
    }

    [Test]
    public void MedianPolish_WritesQuantifiedProteinsTsv_WhenOutGiven()
    {
        string[] runs = { "run_1", "run_2" };
        string path = WritePeptideTable(runs, ("PEPTIDEK", "P1", new[] { 1000.0, 1000 }));
        string outDir = Path.Combine(_tempDirectory, "out");

        RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path, "--out", outDir);

        Assert.That(File.Exists(Path.Combine(outDir, "QuantifiedProteins.tsv")), Is.True);
    }

    [Test]
    public void MedianPolish_DesignNamingAnUnknownRun_IsRejected()
    {
        string path = WritePeptideTable(new[] { "run_1" }, ("PEPTIDEK", "P1", new[] { 1000.0 }));

        Assert.Throws<Program.UsageException>(
            () => RunMedianPolish("not_a_run\tcontrol\n", "quant", "median-polish", "--peptides", path));
    }

    [Test]
    public void MedianPolish_DesignMissingARunInTheTable_IsRejected()
    {
        string path = WritePeptideTable(new[] { "run_1", "run_2" }, ("PEPTIDEK", "P1", new[] { 1000.0, 1000 }));

        // Design mentions run_1 but not run_2: an ambiguous replicate assignment, so it must fail.
        Assert.Throws<Program.UsageException>(
            () => RunMedianPolish("run_1\tcontrol\n", "quant", "median-polish", "--peptides", path));
    }

    [Test]
    public void MedianPolish_TableMissingARequiredColumn_IsRejected()
    {
        string path = Path.Combine(_tempDirectory, "QuantifiedPeptides.tsv");
        // No "Protein Groups" column.
        File.WriteAllLines(path, new[]
        {
            "Sequence\tBase Sequence\tIntensity_run_1\tDetection Type_run_1",
            "PEPTIDEK\tPEPTIDEK\t1000\tMSMS",
        });

        Assert.Throws<Program.UsageException>(
            () => RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path));
    }

    [Test]
    public void MedianPolish_TableWithNoIntensityColumns_IsRejected()
    {
        string path = Path.Combine(_tempDirectory, "QuantifiedPeptides.tsv");
        File.WriteAllLines(path, new[] { PeptideTableLead, "PEPTIDEK\tPEPTIDEK\tP1\tGENE\tHomo sapiens" });

        Assert.Throws<Program.UsageException>(
            () => RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", path));
    }

    [Test]
    public void MedianPolish_MissingPeptidesFile_IsRejected()
    {
        string missing = Path.Combine(_tempDirectory, "not-here.tsv");

        Assert.Throws<Program.UsageException>(
            () => RunMedianPolish(string.Empty, "quant", "median-polish", "--peptides", missing));
    }

    [Test]
    public void MedianPolish_StripsAUtf8BomFromTheDesignStdin()
    {
        string path = WritePeptideTable(new[] { "run_1" }, ("PEPTIDEK", "P1", new[] { 1000.0 }));

        // A BOM prefixing the first design line must not turn "run_1" into "﻿run_1", which would
        // then match no column and be rejected. With the BOM stripped this quantifies cleanly.
        JsonElement root = RunMedianPolish("﻿run_1\tcontrol\n", "quant", "median-polish", "--peptides", path);
        Assert.That(ProteinsByName(root)["P1"].GetProperty("intensities").TryGetProperty("control_1", out _), Is.True);
    }
}
