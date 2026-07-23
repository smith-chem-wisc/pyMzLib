using System.Diagnostics.CodeAnalysis;
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
}
