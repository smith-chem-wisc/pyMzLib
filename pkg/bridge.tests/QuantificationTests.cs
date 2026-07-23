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
}
