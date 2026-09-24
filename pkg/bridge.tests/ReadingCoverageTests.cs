using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Readers;

namespace MzLibBridge.Tests;

/// <summary>
/// The exhaustive-coverage guarantee: every file type mzLib recognises is readable through the
/// bridge, and the three typed views project the interfaces they claim to.
/// </summary>
/// <remarks>
/// <para>
/// The headline test here is
/// <see cref="EverySupportedFileType_IsReadableThroughReadRecords"/>, which is driven by
/// <see cref="Enum.GetValues{TEnum}"/> over <see cref="SupportedFileType"/> rather than by a list
/// maintained here. That matters: a file type added to mzLib arrives in this suite as a failing
/// test demanding a fixture, not as a silent hole in the coverage claim. The same reasoning that
/// makes <c>readers formats</c> enumerate rather than transcribe applies to the promise that it
/// can read all of them.
/// </para>
/// <para>
/// These tests need real bytes, unlike most of <see cref="ReadingTests"/>: dispatch can be faked
/// with an empty file of the right name, but parsing cannot. Each is ignored rather than failed
/// when its fixture is absent, and <see cref="EveryFileTypeHasAFixture"/> fails if the map itself
/// develops a hole — so a missing fixture can never quietly reduce what "exhaustive" means.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReadingCoverageTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-coverage-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void RemoveTempDirectory()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    // ---- the fixture map -----------------------------------------------------------------------

    /// <summary>The pinned mzLib worktree's test tree, resolved relative to this source file.</summary>
    private static string TestRoot([CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "code", "mzLib", "mzLib", "Test");
    }

    /// <summary>
    /// One real file per <see cref="SupportedFileType"/>, relative to mzLib's own test tree.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Deliberately the SMALLEST usable fixture for each type: the whole map is about 1.2 MB
    /// excluding timsTOF, so the suite reads real data without becoming slow.
    /// </para>
    /// <para>
    /// <b><c>Tsv_Dinosaur</c> is the one type with no directly usable fixture.</b> mzLib's own file
    /// is named <c>DinoSnippet.features.tsv</c> — plural — and mzLib dispatches Dinosaur on
    /// <c>.feature.tsv</c>, so <see cref="SupportedFileTypeExtensions.ParseFileType"/> cannot
    /// resolve it and falls through to the generic <c>.tsv</c> branch. It is also the only enum
    /// member absent from mzLib's own <c>TestSupportedFileExtensions</c> cases, which is presumably
    /// how the naming slipped through. The fixture is therefore copied to a correctly-named
    /// temporary file by <see cref="FixtureFor"/> rather than skipped, so the coverage claim holds
    /// for all thirty-six; the upstream fixture is tracked in bridge/UPSTREAM.md.
    /// </para>
    /// </remarks>
    private static readonly Dictionary<SupportedFileType, string> Fixtures = new()
    {
        [SupportedFileType.Ms1Feature] = "FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature",
        [SupportedFileType.Ms2Feature] = "FileReadingTests/ExternalFileTypes/Ms2Feature_FlashDeconvOpenMs3.0.0_ms2.feature",
        [SupportedFileType.TopFDMzrt] = "FileReadingTests/ExternalFileTypes/mzrt_TopFDv1.6.2.mzrt.csv",
        [SupportedFileType.Ms1Tsv_FlashDeconv] = "FileReadingTests/ExternalFileTypes/Ms1Tsv_FlashDeconvOpenMs3.0.0_ms1.tsv",
        [SupportedFileType.Tsv_FlashDeconv] = "FileReadingTests/ExternalFileTypes/Tsv_FlashDeconvOpenMs3.0.0.tsv",
        [SupportedFileType.Tsv_Dinosaur] = "FileReadingTests/ExternalFileTypes/DinoSnippet.features.tsv",
        [SupportedFileType.ThermoRaw] = "DataFiles/ScanDescriptionTestData.raw",
        [SupportedFileType.MzML] = "DataFiles/sliced_ethcd.mzML",
        [SupportedFileType.Mgf] = "DataFiles/withZeros.mgf",
        [SupportedFileType.Ms1Align] = "FileReadingTests/ExternalFileTypes/Ms1Align_FlashDeconvOpenMs3.0.0_ms1.msalign",
        [SupportedFileType.Ms2Align] = "FileReadingTests/ExternalFileTypes/Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign",
        [SupportedFileType.psmtsv] = "FileReadingTests/SearchResults/ExcelEditedPeptide.psmtsv",
        [SupportedFileType.osmtsv] = "Transcriptomics/TestData/OsmWithCustomMIons.osmtsv",
        [SupportedFileType.ToppicPrsm] = "FileReadingTests/ExternalFileTypes/ToppicPrsm_TopPICv1.6.2_prsm.tsv",
        [SupportedFileType.ToppicPrsmSingle] = "FileReadingTests/ExternalFileTypes/ToppicPrsmSingle_TopPICv1.6.2_prsm_single.tsv",
        [SupportedFileType.ToppicProteoform] = "FileReadingTests/ExternalFileTypes/ToppicProteofrom_TopPICv1.6.2_proteoform.tsv",
        [SupportedFileType.ToppicProteoformSingle] = "FileReadingTests/ExternalFileTypes/ToppicProteoformSingle_TopPICv1.5.3_proteoform_single.tsv",
        [SupportedFileType.MsFraggerPsm] = "FileReadingTests/ExternalFileTypes/FraggerPsm_FragPipev21.1_psm.tsv",
        [SupportedFileType.MsFraggerPeptide] = "FileReadingTests/ExternalFileTypes/FraggerPeptide_FragPipev21.1individual_peptide.tsv",
        [SupportedFileType.MsFraggerProtein] = "FileReadingTests/ExternalFileTypes/FraggerProtein_FragPipev21.1individual_protein.tsv",
        [SupportedFileType.FlashLFQQuantifiedPeak] = "FileReadingTests/ExternalFileTypes/FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv",
        [SupportedFileType.MsPathFinderTTargets] = "FileReadingTests/ExternalFileTypes/MsPathFinderT_TargetResults_IcTarget.tsv",
        [SupportedFileType.MsPathFinderTDecoys] = "FileReadingTests/ExternalFileTypes/MsPathFinderT_DecoyResults_IcDecoy.tsv",
        [SupportedFileType.MsPathFinderTAllResults] = "FileReadingTests/ExternalFileTypes/MsPathFinderT_AllResults_IcTda.tsv",
        [SupportedFileType.CruxResult] = "FileReadingTests/ExternalFileTypes/crux.txt",
        [SupportedFileType.ExperimentAnnotation] = "FileReadingTests/ExternalFileTypes/EditedMSFraggerResults/experiment_annotation.tsv",
        [SupportedFileType.BrukerD] = "DataFiles/centroid_1x_MS1_4x_autoMS2.d",
        [SupportedFileType.BrukerTimsTof] = "DataFiles/timsTOF_snippet.d",
        [SupportedFileType.CasanovoMzTab] = "FileReadingTests/ExternalFileTypes/Casanovo_5.0.0.mztab",
        [SupportedFileType.DiaNnReport] = "FileReadingTests/ExternalFileTypes/DiaNn_LongFormat_report.tsv",
        // The smallest of mzLib's three SDRF corpus files; the other two (PXD026824, PXD059974)
        // exercise the validator and the cross-document lint, which are not this suite's subject.
        [SupportedFileType.Sdrf] = "FileReadingTests/ExternalFileTypes/PXD000070.sdrf.tsv",
        // mzLib's only Pytheas file is 5.5 MB; FixtureFor trims it to its head (see there).
        [SupportedFileType.PytheasResult] = "FileReadingTests/ExternalFileTypes/match_output_Lumos_Orbi.txt",
        // mzLib #1313. The same MS-GF+ PRIDE cut in both forms, so the two types are compared on the
        // same records; the .gz is 3 KB. MS-GF+ is the writer whose scores and q-value the reader
        // maps (#1306).
        [SupportedFileType.MzIdentML] = "DataFiles/PXD078927_msgf_1_1_0.mzid",
        [SupportedFileType.MzIdentMLGz] = "DataFiles/PXD078927_msgf_1_1_0.mzid.gz",
        // mzLib #1347. MetaMorpheus 1.1.11 output, the only version mzLib ships a fixture for.
        [SupportedFileType.MetaMorpheusQuantifiedProteinGroups] = "FileReadingTests/ExternalFileTypes/MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv",
        [SupportedFileType.FlashLFQQuantifiedPeptide] = "FileReadingTests/ExternalFileTypes/MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv",
    };

    /// <summary>
    /// The two file types read through vendor native libraries, which are Windows-x64 only.
    /// </summary>
    /// <remarks>
    /// Bruker's <c>baf2sql_c.dll</c> and <c>timsdata.dll</c> are P/Invoked, so these two genuinely
    /// cannot be read on Linux or macOS. Ignoring them off-Windows is honest; silently dropping
    /// them from the coverage claim would not be, which is why they are named here rather than
    /// omitted from <see cref="Fixtures"/>.
    /// </remarks>
    private static readonly HashSet<SupportedFileType> WindowsOnly =
    [
        SupportedFileType.BrukerD,
        SupportedFileType.BrukerTimsTof,
    ];

    private static IEnumerable<SupportedFileType> AllFileTypes() => Enum.GetValues<SupportedFileType>();

    [Test]
    public void EveryFileTypeHasAFixture()
    {
        // Enumerated, not transcribed: a type added to mzLib fails here demanding a fixture rather
        // than quietly shrinking what the exhaustive-coverage claim covers.
        List<SupportedFileType> missing = AllFileTypes().Where(t => !Fixtures.ContainsKey(t)).ToList();

        Assert.That(missing, Is.Empty,
            "Every SupportedFileType needs a fixture for the exhaustive-coverage test. Missing: " +
            string.Join(", ", missing));
    }

    [Test]
    public void FixtureRoot_Exists_SoTheCoverageSuiteCannotSilentlyVanish()
    {
        Assert.That(Directory.Exists(TestRoot()), Is.True,
            $"mzLib test root not found at {TestRoot()} — every coverage test would ignore and the " +
            "run would look clean. Check the pinned mzLib worktree.");
    }

    /// <summary>The fixture path for a type, ignoring the test when it is absent or unusable here.</summary>
    private string FixtureFor(SupportedFileType fileType)
    {
        if (WindowsOnly.Contains(fileType) && !OperatingSystem.IsWindows())
            Assert.Ignore($"{fileType} is read through Windows-x64 vendor native libraries.");

        string path = Path.Combine(TestRoot(), Fixtures[fileType]);
        if (!File.Exists(path) && !Directory.Exists(path))
            Assert.Ignore($"mzLib fixture not present in the worktree: {path}");

        // See the Fixtures remarks: mzLib's Dinosaur fixture carries a name mzLib itself cannot
        // dispatch, so it is copied to a dispatchable one rather than dropped from the sweep.
        if (fileType == SupportedFileType.Tsv_Dinosaur)
        {
            string renamed = Path.Combine(_tempDirectory, "DinoSnippet.feature.tsv");
            File.Copy(path, renamed, overwrite: true);
            return renamed;
        }

        // mzLib's Pytheas fixture is 5.5 MB, over four times the rest of the map combined. Its first
        // 40 lines are the whole '#' parameter header plus the first several PRECURSOR_ION groups,
        // which is every line shape the reader parses. Each match line is self-contained, so a cut
        // between groups leaves a valid file.
        if (fileType == SupportedFileType.PytheasResult)
        {
            string head = Path.Combine(_tempDirectory, "match_output_head.txt");
            File.WriteAllLines(head, File.ReadLines(path).Take(40));
            return head;
        }

        return path;
    }

    // ---- the exhaustive-coverage guarantee -----------------------------------------------------

    [Test]
    [TestCaseSource(nameof(AllFileTypes))]
    public void EverySupportedFileType_IsReadableThroughReadRecords(SupportedFileType fileType)
    {
        string path = FixtureFor(fileType);

        JsonElement data = Invoke("readers", "read-records", "--path", path, "--limit", "2");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("file_type").GetString(), Is.EqualTo(fileType.ToString()),
                "read-records must report the type mzLib dispatched, not a guess from the extension.");
            Assert.That(data.GetProperty("record_type").GetString(), Is.Not.Empty,
                "Every read names the mzLib record class its columns came from.");
            Assert.That(data.GetProperty("column_names").GetArrayLength(), Is.GreaterThan(0),
                "A format with no columns would mean the projection found nothing to report.");
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.GreaterThan(0),
                "Every fixture holds records; a zero count means the file parsed to nothing.");
        });
    }

    [Test]
    [TestCaseSource(nameof(AllFileTypes))]
    public void EverySupportedFileType_ReportsColumnsThatMatchItsRecordCount(SupportedFileType fileType)
    {
        string path = FixtureFor(fileType);

        JsonElement data = Invoke("readers", "read-records", "--path", path, "--limit", "2");

        int returned = data.GetProperty("returned_count").GetInt32();
        JsonElement columns = data.GetProperty("columns");

        foreach (JsonProperty column in columns.EnumerateObject())
        {
            // A column shorter than returned_count is the silent-truncation failure this whole
            // module exists to refuse: it would drop rows from one field and no other.
            Assert.That(column.Value.GetArrayLength(), Is.EqualTo(returned),
                $"Column '{column.Name}' has {column.Value.GetArrayLength()} values but " +
                $"{returned} records were returned.");
        }
    }

    [Test]
    [TestCaseSource(nameof(AllFileTypes))]
    public void EverySupportedFileType_WritesATableWithOneHeaderPerColumn(SupportedFileType fileType)
    {
        string path = FixtureFor(fileType);
        string output = Path.Combine(_tempDirectory, $"{fileType}.tsv");

        JsonElement data = Invoke("readers", "read-records", "--path", path, "--limit", "2", "--out", output);

        Assert.That(File.Exists(output), Is.True, "--out must write the table it reports.");

        string[] header = File.ReadLines(output).First().Split('\t');
        Assert.Multiple(() =>
        {
            Assert.That(header, Is.EqualTo(data.GetProperty("column_names").EnumerateArray()
                    .Select(c => c.GetString()).ToArray()),
                "The written header and the reported column_names must not be able to disagree.");
            Assert.That(data.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null),
                "Writing to disk must omit the inline payload, or both are materialised.");
        });
    }

    // ---- the projection's honesty rules --------------------------------------------------------

    [Test]
    public void CompositeFieldsAreNamedAsExcluded_NotSilentlyDropped()
    {
        // ToppicPrsm carries List<AlternativeToppicId>, which has no faithful column shape.
        JsonElement data = Invoke("readers", "read-records",
            "--path", FixtureFor(SupportedFileType.ToppicPrsm), "--limit", "1");

        List<string?> excluded = data.GetProperty("excluded_fields").EnumerateArray()
            .Select(e => e.GetProperty("field").GetString()).ToList();

        Assert.That(excluded, Does.Contain("alternative_identifications"),
            "A field that cannot cross the wire must be named with its reason. A column that simply " +
            "vanished is indistinguishable from a field the format does not have.");
        Assert.That(
            data.GetProperty("excluded_fields")[0].GetProperty("reason").GetString(),
            Is.Not.Empty, "Every exclusion carries the reason for it.");
    }

    [TestCase(SupportedFileType.MetaMorpheusQuantifiedProteinGroups, "sample_groups")]
    [TestCase(SupportedFileType.FlashLFQQuantifiedPeptide, "samples")]
    [TestCase(SupportedFileType.MzIdentML, "scores")]
    public void DictionaryFieldsAreExcludedAsDictionaries(SupportedFileType fileType, string field)
    {
        // mzLib 1.0.592's new readers keep their most important values in read-only dictionaries:
        // per-sample intensity, spectral count and occupancy for a protein group (#1347), per-run
        // intensity for a peptide (#1347), and the engine's scores for an mzIdentML match (#1313).
        // read-records cannot project them, so it must say so, and say it is a dictionary:
        // IReadOnlyDictionary<K, V> is not a non-generic IDictionary, and before the generic check
        // these were excluded as "a list of composite values", which misdescribes them.
        JsonElement data = Invoke("readers", "read-records", "--path", FixtureFor(fileType), "--limit", "1");

        JsonElement? exclusion = data.GetProperty("excluded_fields").EnumerateArray()
            .Cast<JsonElement?>()
            .FirstOrDefault(e => e!.Value.GetProperty("field").GetString() == field);

        Assert.That(exclusion, Is.Not.Null, $"{fileType}.{field} must be named in excluded_fields");
        Assert.That(exclusion!.Value.GetProperty("reason").GetString(), Does.StartWith("a dictionary"));
    }

    [Test]
    public void AcronymsAndPluralsBecomeReadableColumnNames()
    {
        JsonElement data = Invoke("readers", "read-records",
            "--path", FixtureFor(SupportedFileType.ToppicPrsm), "--limit", "1");

        List<string?> columns = data.GetProperty("column_names").EnumerateArray()
            .Select(c => c.GetString()).ToList();

        Assert.Multiple(() =>
        {
            // EValue, not E_Value: a single leading capital is its own word.
            Assert.That(columns, Does.Contain("e_value"));
            // MIScore, not m_i_score: a run of capitals is one word.
            Assert.That(columns, Does.Contain("mi_score"));
            // FixedPTMs, not fixed_pt_ms: a pluralising 's' belongs to the acronym before it.
            Assert.That(columns, Does.Contain("fixed_ptms"));
        });
    }

    [Test]
    public void NonFiniteNumbersCrossAsNull_RatherThanBreakingTheEnvelope()
    {
        // Regression: an mzML whose scan window is unbounded reports infinity for its bound, and
        // System.Text.Json refuses to serialise it — so the whole read failed with a message about
        // JSON rather than about the file. tester.mzML is such a file.
        string path = Path.Combine(TestRoot(), "DataFiles", "tester.mzML");
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present: {path}");

        JsonElement data = Invoke("readers", "read-spectra", "--path", path, "--limit", "1");

        Assert.That(data.GetProperty("record_count").GetInt32(), Is.GreaterThan(0),
            "A non-finite value must null the cell, not fail the read.");
    }

    [Test]
    public void NegativeModeMgf_ReportsANegativeChargeAndNegativePolarity()
    {
        // Regression for mzLib #1164. An MGF CHARGE line carries its sign as a TRAILING character
        // ("CHARGE=2-"), which the reader previously dropped: a negative-mode precursor arrived as
        // charge +2 with positive polarity. Nothing in the envelope disclosed it, and the sign is
        // not recoverable downstream -- a caller computing a neutral mass from selected_ion_mz and
        // selected_ion_charge_state_guess got an answer wrong by two proton masses that looked
        // entirely ordinary.
        //
        // Pinned here rather than only upstream because both wrong values are columns this verb
        // publishes, so the projection is what a pyMzLib caller actually sees.
        //
        // The fixture is written here rather than taken from the mzLib worktree ON PURPOSE.
        // mzLib ships negativeModeCharge.mgf, but it arrived WITH the fix, so a worktree-sourced
        // test cannot fail against a pre-#1164 mzLib -- it skips for a missing fixture and reports
        // green. Six lines inline make this fail where it should. Verified: against pin 5ba13155
        // this asserts -2 and gets 2.
        string path = Path.Combine(Path.GetTempPath(), $"pymzlib-negative-{Guid.NewGuid():N}.mgf");
        File.WriteAllText(path, """
            BEGIN IONS
            TITLE=negative mode scan
            PEPMASS=571.806916 999999
            CHARGE=2-
            RTINSECONDS=16
            SCANS=1
            110.0719 3823.8
            130.0863 5897.2
            END IONS
            BEGIN IONS
            TITLE=positive mode scan
            PEPMASS=684.312500 888888
            CHARGE=3+
            RTINSECONDS=32
            SCANS=2
            115.0866 2100.4
            310.4567 4400.1
            END IONS
            """);

        try
        {
            JsonElement data = Invoke("readers", "read-spectra", "--path", path);
            JsonElement columns = data.GetProperty("columns");

            string[] polarity = columns.GetProperty("polarity").EnumerateArray()
                .Select(x => x.GetString()!).ToArray();
            int?[] charge = columns.GetProperty("selected_ion_charge_state_guess").EnumerateArray()
                .Select(x => x.ValueKind == JsonValueKind.Null ? (int?)null : x.GetInt32()).ToArray();

            // One negative-mode scan followed by one positive-mode scan, so the sign is pinned in
            // both directions: a reader that simply negated everything would pass on scan 1 alone.
            Assert.Multiple(() =>
            {
                Assert.That(charge, Is.EqualTo(new int?[] { -2, 3 }),
                    "CHARGE=2- must read as -2 - mzLib #1164 has regressed");
                Assert.That(polarity, Is.EqualTo(new[] { "Negative", "Positive" }));
            });
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void AParallelReaderFailureReportsItsRealCause_NotTheAggregateWrapper()
    {
        // mzLib parallelises its mzML reader, so "profile mode is unsupported" arrives wrapped in
        // an AggregateException. Unwrapped, the caller sees what actually went wrong; wrapped, they
        // see "One or more errors occurred."
        string path = Path.Combine(TestRoot(), "DataFiles", "tiny.pwiz.1.1.mzML");
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present: {path}");

        JsonElement error = InvokeExpectingError("readers", "read-spectra", "--path", path);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.Not.EqualTo("AggregateException"),
                "The wrapper type tells a caller nothing about what failed.");
            Assert.That(error.GetProperty("message").GetString(), Does.Not.Contain("One or more errors"),
                "The message must name the real cause.");
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("profile"));
        });
    }

    // ---- the typed views -----------------------------------------------------------------------

    [Test]
    public void ReadFeatures_ExpandsMs1FeatureRowsAcrossChargeStates()
    {
        JsonElement records = Invoke("readers", "read-records",
            "--path", FixtureFor(SupportedFileType.Ms1Feature), "--limit", "1");
        JsonElement features = Invoke("readers", "read-features",
            "--path", FixtureFor(SupportedFileType.Ms1Feature), "--limit", "1");

        Assert.That(features.GetProperty("record_count").GetInt32(),
            Is.GreaterThan(records.GetProperty("record_count").GetInt32()),
            "mzLib expands each _ms1.feature row into one feature per charge state, so the feature " +
            "view must report more rows than the file has. If these ever match, the caveat saying " +
            "so has become a lie.");
    }

    [Test]
    public void ReadFeatures_DoesNotExpandDinosaurRows()
    {
        // The counterpart, and the reason the expansion caveat is per-format:
        // DinosaurTsvFile.GetMs1Features is `=> Results`, one feature per row.
        string path = FixtureFor(SupportedFileType.Tsv_Dinosaur);

        JsonElement records = Invoke("readers", "read-records", "--path", path, "--limit", "1");
        JsonElement features = Invoke("readers", "read-features", "--path", path, "--limit", "1");

        Assert.That(features.GetProperty("record_count").GetInt32(),
            Is.EqualTo(records.GetProperty("record_count").GetInt32()),
            "Dinosaur features are one-for-one with its rows, so the expansion caveat must not be " +
            "emitted for it.");

        List<string?> caveats = features.GetProperty("caveats").EnumerateArray()
            .Select(c => c.GetString()).ToList();
        Assert.That(caveats.Any(c => c!.Contains("CHARGE STATE")), Is.False,
            "The expansion caveat is false for Dinosaur and must not appear.");
    }

    [Test]
    public void Ms1FeatureRetentionTimeUnitIsUnknown_BecauseTopFdChangedItMidVersion()
    {
        JsonElement data = Invoke("readers", "read-features",
            "--path", FixtureFor(SupportedFileType.Ms1Feature), "--limit", "1");

        Assert.That(data.GetProperty("retention_time_unit").GetString(), Is.EqualTo("unknown"),
            "TopFD wrote seconds through v1.6.2 and minutes from v1.7.0 without changing the file " +
            "type. Claiming either would be a guess; mzLib's own deconvolution code guesses here " +
            "and this deliberately does not.");
    }

    [Test]
    public void DiaNnRetentionTimeUnitIsMinutes_BecauseDiaNnWritesMinutes()
    {
        JsonElement data = Invoke("readers", "read-results",
            "--path", FixtureFor(SupportedFileType.DiaNnReport), "--limit", "1");

        Assert.That(data.GetProperty("retention_time_unit").GetString(), Is.EqualTo("minutes"),
            "DIA-NN writes retention times in minutes and mzLib converts nothing, which its own " +
            "reader states at DiaNnPrecursor.RetentionTime. Letting this fall through to 'unknown' " +
            "would understate what is actually known, and would make the guide's claim that every " +
            "quantifiable format reports minutes false.");
    }

    [Test]
    public void Ms1FeatureFixturesShowTheUnitGenuinelyChanged()
    {
        // The evidence for the caveat above, pinned so it cannot rot into folklore: the same file
        // type, two TopFD versions, retention times three orders of magnitude apart.
        string older = Path.Combine(TestRoot(),
            "FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.6.2_ms1.feature");
        string newer = Path.Combine(TestRoot(),
            "FileReadingTests/ExternalFileTypes/Ms1Feature_TopFDv1.7.0_ms1.feature");
        if (!File.Exists(older) || !File.Exists(newer))
            Assert.Ignore("Both TopFD fixtures are needed to pin the unit change.");

        double FirstStart(string path) => Invoke("readers", "read-features", "--path", path, "--limit", "1")
            .GetProperty("columns").GetProperty("retention_time_start")[0].GetDouble();

        Assert.Multiple(() =>
        {
            Assert.That(FirstStart(older), Is.GreaterThan(600),
                "v1.6.2 writes seconds — a value beyond any plausible gradient length in minutes.");
            Assert.That(FirstStart(newer), Is.LessThan(600),
                "v1.7.0 writes minutes.");
        });
    }

    [Test]
    public void CasanovoIsDecoyIsNull_BecauseDeNovoSequencingHasNoDecoys()
    {
        JsonElement data = Invoke("readers", "read-matches",
            "--path", FixtureFor(SupportedFileType.CasanovoMzTab), "--limit", "1");

        Assert.That(data.GetProperty("columns").GetProperty("is_decoy")[0].ValueKind,
            Is.EqualTo(JsonValueKind.Null),
            "mzLib leaves Casanovo's IsDecoy at its default false and never assigns it, so false " +
            "would mean 'unknown' and let a caller filter on a fabricated column — the same trap " +
            "the quantifiable view already refuses for MSFragger.");
    }

    [Test]
    public void FlashLFQPeaks_AbsentMbrScoreIsNull_NotZero()
    {
        // mzLib #1345 made QuantifiedPeak.MBRScore a double?. The 1.0.549 fixture's MBR Score cells
        // are blank on its MSMS peaks; they read 0 before, a number that says "scored, and scored
        // nothing", and must read null now.
        JsonElement data = Invoke("readers", "read-records",
            "--path", FixtureFor(SupportedFileType.FlashLFQQuantifiedPeak));

        foreach (JsonElement score in data.GetProperty("columns").GetProperty("mbr_score").EnumerateArray())
            Assert.That(score.ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void FlashLFQPeaks_CurrentFormatWithoutAnMbrScoreColumn_Reads()
    {
        // mzLib #1345: current FlashLFQ (and MetaMorpheus 1.1.11) peaks tables have no MBR Score
        // column and add seven others. They threw HeaderValidationException through 1.0.591.
        string path = Path.Combine(TestRoot(), "FileReadingTests", "ExternalFileTypes",
            "FlashLFQ_MzLib1.0.591_QuantifiedPeaks.tsv");
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present: {path}");

        JsonElement data = Invoke("readers", "read-records", "--path", path);
        List<string?> columns = data.GetProperty("column_names").EnumerateArray().Select(c => c.GetString()).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.GreaterThan(0));
            Assert.That(columns, Is.SupersetOf(new[]
            {
                "organism", "peak_fwhm", "peak_fwhm_status", "pip_q_value", "pip_pep",
                "decoy_peptide", "random_rt",
            }));
        });
    }

    [TestCase(SupportedFileType.MzIdentML)]
    [TestCase(SupportedFileType.MzIdentMLGz)]
    public void MzIdentMLIsDecoyIsNull_AndItsCaveatsSaySo(SupportedFileType fileType)
    {
        // mzIdentML's isDecoy attribute is optional and defaults to false, so a writer that omits it
        // reads as a target. Kept off the DecoysAreReported allowlist for the same reason as
        // Casanovo; read-records still carries mzLib's own boolean.
        JsonElement data = Invoke("readers", "read-matches", "--path", FixtureFor(fileType), "--limit", "2");

        List<string> caveats = data.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()!).ToList();
        Assert.Multiple(() =>
        {
            foreach (JsonElement flag in data.GetProperty("columns").GetProperty("is_decoy").EnumerateArray())
                Assert.That(flag.ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(caveats, Has.Some.StartsWith("is_decoy is null for this format. mzIdentML"));
            Assert.That(caveats, Has.Some.Contains("skipped_count"),
                "mzLib drops items it cannot represent; the caveat must point at where they are reported.");
            Assert.That(data.GetProperty("absent_fields").EnumerateArray().Select(f => f.GetString()),
                Does.Contain("is_decoy"),
                "A null the format cannot fill is absent, and must be named as such.");
        });
    }

    [Test]
    public void MsPathFinderTReportsRealDecoyFlags()
    {
        // The counterpart to the Casanovo test, and it has to assert the flag is RIGHT rather than
        // merely non-null: "not null" would pass for an implementation that fabricated false for
        // every row, which is the exact under-report this test is named for. The decoy fixture's
        // every protein name starts with XXX, so every flag must be true.
        JsonElement decoys = Invoke("readers", "read-matches",
            "--path", FixtureFor(SupportedFileType.MsPathFinderTDecoys));
        JsonElement targets = Invoke("readers", "read-matches",
            "--path", FixtureFor(SupportedFileType.MsPathFinderTTargets));

        Assert.Multiple(() =>
        {
            foreach (JsonElement flag in decoys.GetProperty("columns").GetProperty("is_decoy").EnumerateArray())
                Assert.That(flag.ValueKind, Is.EqualTo(JsonValueKind.True),
                    "every row of the decoy fixture is a decoy");

            foreach (JsonElement flag in targets.GetProperty("columns").GetProperty("is_decoy").EnumerateArray())
                Assert.That(flag.ValueKind, Is.EqualTo(JsonValueKind.False),
                    "every row of the target fixture is a target");
        });
    }

    [Test]
    public void ReadSpectra_FiltersByMsOrderBeforeTheWindow()
    {
        string path = FixtureFor(SupportedFileType.MzML);

        JsonElement all = Invoke("readers", "read-spectra", "--path", path);
        JsonElement ms2 = Invoke("readers", "read-spectra", "--path", path, "--ms-order", "2");

        Assert.Multiple(() =>
        {
            Assert.That(ms2.GetProperty("scan_count").GetInt32(),
                Is.EqualTo(all.GetProperty("scan_count").GetInt32()),
                "scan_count reports the file's real total, so a filter that matched nothing can " +
                "never look like an empty file.");
            // Asserted before the loop below, which has an empty body — and so proves nothing —
            // if the filter happened to match no scans at all.
            Assert.That(ms2.GetProperty("record_count").GetInt32(), Is.GreaterThan(0),
                "the fixture must contain MS2 scans for this test to mean anything");
            Assert.That(ms2.GetProperty("record_count").GetInt32(),
                Is.LessThan(all.GetProperty("record_count").GetInt32()),
                "the fixture must also contain non-MS2 scans, or the filter is untested");
            foreach (JsonElement order in ms2.GetProperty("columns").GetProperty("ms_order").EnumerateArray())
                Assert.That(order.GetInt32(), Is.EqualTo(2));
        });
    }

    [Test]
    public void ReadSpectra_OmitsPeaksUnlessAsked()
    {
        string path = FixtureFor(SupportedFileType.MzML);

        JsonElement without = Invoke("readers", "read-spectra", "--path", path, "--limit", "1");
        JsonElement with = Invoke("readers", "read-spectra", "--path", path, "--limit", "1", "--peaks");

        Assert.Multiple(() =>
        {
            Assert.That(without.GetProperty("peaks_included").GetBoolean(), Is.False);
            Assert.That(without.GetProperty("columns").TryGetProperty("mz", out _), Is.False,
                "Peaks must be absent by default: a mid-size mzML would otherwise serialise " +
                "hundreds of megabytes for the ordinary 'what is in this file' call.");
            Assert.That(with.GetProperty("peaks_included").GetBoolean(), Is.True);
            Assert.That(with.GetProperty("columns").GetProperty("mz")[0].GetArrayLength(),
                Is.EqualTo(with.GetProperty("columns").GetProperty("peak_count")[0].GetInt32()),
                "The peak array and the reported peak count must agree.");
        });
    }

    [Test]
    public void AViewVerbRejectsAFileWithoutThatView_AndNamesWhatToUseInstead()
    {
        JsonElement error = InvokeExpectingError("readers", "read-features",
            "--path", FixtureFor(SupportedFileType.psmtsv));

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"),
                "Asking for a view a file does not have is the caller's mistake, not a fault.");
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("quantifiable"),
                "The message names the views the file DOES have.");
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("read-records"),
                "…and the verb that can read it regardless.");
        });
    }

    [TestCase("read-records")]
    [TestCase("read-features")]
    [TestCase("read-matches")]
    [TestCase("read-spectra")]
    public void EveryReadVerbRefusesToOverwriteItsInput(string verb)
    {
        // The rule was argued out for read-results; sharing one Window parser is what stops the
        // four new verbs from each having to remember it.
        string path = FixtureFor(SupportedFileType.psmtsv);

        JsonElement error = InvokeExpectingError("readers", verb, "--path", path, "--out", path);

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("must differ from"),
            $"{verb} must refuse to write its projection over the file it is reading.");
    }

    [TestCase("read-records")]
    [TestCase("read-features")]
    [TestCase("read-matches")]
    [TestCase("read-spectra")]
    public void EveryReadVerbRejectsAnOptionGivenWithoutAValue(string verb)
    {
        string path = FixtureFor(SupportedFileType.psmtsv);

        JsonElement error = InvokeExpectingError("readers", verb, "--path", path, "--limit");

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("no value"),
            $"{verb} must reject '--limit' with no value rather than silently returning everything.");
    }

    [Test]
    public void AFlashDeconvFeatureFileReportsNullIntensityRatherThanFabricatedZero()
    {
        // A within-type schema variant, which one fixture per SupportedFileType cannot cover:
        // Apex_intensity is [Optional] and the FLASHDeconv/OpenMS _ms1.feature layout omits it, so
        // mzLib's `IntensityApex ?? 0` hands back a whole column of zeros that look exactly like
        // measurements of nothing. The TopFD fixture in the map has the column and hides this.
        string path = Path.Combine(TestRoot(),
            "FileReadingTests/ExternalFileTypes/Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature");
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present: {path}");

        JsonElement data = Invoke("readers", "read-features", "--path", path, "--limit", "5");

        Assert.Multiple(() =>
        {
            foreach (JsonElement intensity in data.GetProperty("columns").GetProperty("intensity").EnumerateArray())
                Assert.That(intensity.ValueKind, Is.EqualTo(JsonValueKind.Null),
                    "a fabricated zero must not be handed back as a measurement");

            Assert.That(
                data.GetProperty("caveats").EnumerateArray()
                    .Any(caveat => caveat.GetString()!.Contains("intensity is NULL")),
                Is.True,
                "and the reason must be stated, not left as an unexplained empty column");
        });
    }

    [Test]
    public void ATopFdFeatureFileStillReportsRealIntensities()
    {
        // The counterpart. Nulling every intensity would be an equally serious over-correction, and
        // this is the fixture that proves the null above is conditional rather than blanket.
        JsonElement data = Invoke("readers", "read-features",
            "--path", FixtureFor(SupportedFileType.Ms1Feature), "--limit", "3");

        foreach (JsonElement intensity in data.GetProperty("columns").GetProperty("intensity").EnumerateArray())
            Assert.That(intensity.ValueKind, Is.EqualTo(JsonValueKind.Number));
    }

    [Test]
    public void TheTwoVerbsAgreeOnTheAbsentSentinelForTheSameColumn()
    {
        // read-records reaches the SAME IQuantifiableRecord properties that read-results nulls, so
        // without the scoped exception one verb would answer null and the other -1 for the same
        // column of the same file, and the -1 would enter a mean.
        string path = FixtureFor(SupportedFileType.psmtsv);

        JsonElement uniform = Invoke("readers", "read-results", "--path", path);
        JsonElement native = Invoke("readers", "read-records", "--path", path);

        JsonElement[] fromResults = uniform.GetProperty("columns").GetProperty("retention_time")
            .EnumerateArray().ToArray();
        JsonElement[] fromRecords = native.GetProperty("columns").GetProperty("retention_time")
            .EnumerateArray().ToArray();

        Assert.That(fromRecords, Has.Length.EqualTo(fromResults.Length));
        for (int row = 0; row < fromResults.Length; row++)
        {
            Assert.That(fromRecords[row].ValueKind, Is.EqualTo(fromResults[row].ValueKind),
                $"row {row}: the two verbs disagree about whether retention_time is present");
        }
    }

    [Test]
    [TestCaseSource(nameof(SpectraFileTypes))]
    public void EverySpectraFormatReportsItsOwnCaveats(SupportedFileType fileType)
    {
        // Driven by the type list rather than run once on mzML: the caveats differ per format and
        // each one is a claim about that format specifically. A single mzML case would leave the
        // msalign neutral-mass warning, the MGF derived-scan-window warning and the vendor
        // native-library warning all unexercised — and those are the three most likely to be wrong.
        string path = FixtureFor(fileType);

        JsonElement data = Invoke("readers", "read-spectra", "--path", path, "--limit", "1");

        List<string> caveats = data.GetProperty("caveats").EnumerateArray()
            .Select(caveat => caveat.GetString()!).ToList();

        Assert.Multiple(() =>
        {
            // Every format gets the peaks-omitted note, because peaks are opt-in for all of them.
            Assert.That(caveats.Any(caveat => caveat.Contains("Peaks are not included")), Is.True);

            switch (fileType)
            {
                case SupportedFileType.Mgf:
                    // mzLib #1165 gave MGF a real MS level: MSLEVEL when the writer supplied one,
                    // otherwise MS2 for a block with a precursor and MS1 for one without. The old claim
                    // that MGF has no MS1 scans is therefore retired rather than adjusted -- it was the
                    // assertion that went wrong, not the reader, and loosening it would have hidden that.
                    // What still has to reach a caller is that ms_order is inferred rather than fixed,
                    // since a caller filtering on it would otherwise assume every MGF scan is MS2.
                    Assert.That(caveats.Any(c => c.Contains("ms_order is no longer always 2")), Is.True);
                    // Derived from the observed peaks rather than recorded by the format.
                    Assert.That(caveats.Any(c => c.Contains("DERIVED")), Is.True);
                    break;

                case SupportedFileType.Ms1Align:
                case SupportedFileType.Ms2Align:
                    Assert.That(caveats.Any(c => c.Contains("DECONVOLVED")), Is.True);
                    // ...and that its scan window is on a different axis from its mz column.
                    Assert.That(caveats.Any(c => c.Contains("neutral")), Is.True);
                    break;

                case SupportedFileType.ThermoRaw:
                    Assert.That(caveats.Any(c => c.Contains("RawFileReader")), Is.True);
                    break;

                case SupportedFileType.BrukerD:
                    Assert.That(caveats.Any(c => c.Contains("Windows-x64")), Is.True);
                    break;

                case SupportedFileType.BrukerTimsTof:
                    Assert.That(caveats.Any(c => c.Contains("Windows-x64")), Is.True);
                    // The mobility axis is collapsed into scans, which is worth saying out loud.
                    Assert.That(caveats.Any(c => c.Contains("mobility")), Is.True);
                    break;

                case SupportedFileType.MzML:
                    // mzML is the one format with no format-specific hazard to report.
                    break;

                default:
                    Assert.Fail($"{fileType} is not a spectra type");
                    break;
            }
        });
    }

    private static IEnumerable<SupportedFileType> SpectraFileTypes() =>
    [
        SupportedFileType.MzML,
        SupportedFileType.Mgf,
        SupportedFileType.Ms1Align,
        SupportedFileType.Ms2Align,
        SupportedFileType.ThermoRaw,
        SupportedFileType.BrukerD,
        SupportedFileType.BrukerTimsTof,
    ];

    [Test]
    [TestCaseSource(nameof(AllFileTypes))]
    public void EveryFormatsCaveatsAreAsciiAndNonEmpty(SupportedFileType fileType)
    {
        // A caveat is read in a terminal and quoted into three bindings' documentation. A stray
        // em-dash renders as mojibake in a Windows console, and an empty string is a caveat that
        // was meant to say something.
        string path = FixtureFor(fileType);

        foreach (string verb in new[] { "read-features", "read-matches", "read-spectra" })
        {
            JsonElement envelope = Envelope(["readers", verb, "--path", path, "--limit", "1"]);
            if (!envelope.GetProperty("ok").GetBoolean())
                continue;

            foreach (JsonElement caveat in envelope.GetProperty("data").GetProperty("caveats").EnumerateArray())
            {
                string text = caveat.GetString()!;
                Assert.That(text, Is.Not.Empty);
                Assert.That(text.All(character => character < 128), Is.True,
                    $"non-ASCII in a {fileType} {verb} caveat: {text}");
            }
        }
    }

    // ---- the caveats' own citations --------------------------------------------------------------

    /// <summary>
    /// Every <c>File.cs:NNN</c> a caveat cites must point at a line that mentions the thing the
    /// caveat is about.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The caveats are the most dangerous text this bridge emits: they read authoritatively, they
    /// are quoted in three bindings' documentation, and nothing else checks them. Two had already
    /// gone stale against the pinned mzLib — mzLib PR #1116 inserted fourteen lines into
    /// <c>MsFraggerPsm.cs</c> and both citations below it silently moved — which is exactly the rot
    /// this guards.
    /// </para>
    /// <para>
    /// It checks the citation is <i>anchored</i>, not that the whole claim is true; a sentence
    /// cannot be verified mechanically. But a citation pointing at a blank line is proof the claim
    /// was not re-checked, and that is worth failing for.
    /// </para>
    /// </remarks>
    [Test]
    public void EveryCitedMzLibLineStillMentionsWhatTheCaveatSaysItDoes()
    {
        // file name -> the token the cited line must contain, per caveat.
        (string Citation, string Token)[] citations =
        [
            ("MsFraggerPsm.cs:231", "IsDecoy"),
            ("MsFraggerPsm.cs:233", "MonoisotopicMass"),
            ("SpectrumMatchFromTsv.cs:119", "MonoisotopicMass"),
            ("SpectrumMatchFromTsv.cs:194", "FullSequence"),
            ("SpectrumMatchTsvReader.cs:71", "catch"),
            ("PsmFromTsvFile.cs:17", "warnings"),
            ("Ms1Feature.cs:84", "ChargeState"),
            ("Ms1Feature.cs:86", "IntensityApex"),
            ("Ms1Feature.cs:91", "SingleChargeMs1Feature"),
            ("DinosaurFeature.cs:18", "Intensity"),
            ("DinosaurTsvFile.cs:15", "Results"),
            ("MsPathFinderTResult.cs:89", "Accession"),
            ("MsPathFinderTResult.cs:92", "XXX"),
            ("CasanovoMzTabRecord.cs:84", "IsDecoy"),
            ("CasanovoMzTabFile.cs:116", "OneBasedScanNumber"),
            ("CasanovoMzTabFile.cs:124", "Modification"),
            ("Mgf.cs:297", "PRECURSORSCAN"),
            ("Mgf.cs:346", "MzRange"),
            ("Mgf.cs:350", "msLevel"),
            ("MsAlign.cs:526", "MzRange"),
            ("MzIdentMLResultFile.cs:123", "skipped"),
            ("MzIdentMLResultFile.cs:159", "OneBasedScanNumber"),
            ("MzIdentMLResultFile.cs:173", "IsDecoy"),
            ("MzIdentMLResultFile.cs:177", "Rank"),
            ("MzIdentMLResultFile.cs:179", "Scores"),
            ("MzIdentMLRecord.cs:45", "Accession"),
        ];

        var wrong = new List<string>();
        foreach ((string citation, string token) in citations)
        {
            string[] parts = citation.Split(':');
            string? file = Directory
                .EnumerateFiles(Path.Combine(TestRoot(), "..", "Readers"), parts[0], SearchOption.AllDirectories)
                .FirstOrDefault();
            if (file is null)
            {
                wrong.Add($"{citation}: no such file under Readers/");
                continue;
            }

            string[] lines = File.ReadAllLines(file);
            int index = int.Parse(parts[1]) - 1;
            if (index < 0 || index >= lines.Length)
            {
                wrong.Add($"{citation}: the file has only {lines.Length} lines");
                continue;
            }

            if (!lines[index].Contains(token, StringComparison.Ordinal))
                wrong.Add($"{citation}: expected '{token}', line reads: {lines[index].Trim()}");
        }

        Assert.That(wrong, Is.Empty,
            "A caveat cites an mzLib line that no longer says what the caveat claims:\n  " +
            string.Join("\n  ", wrong));
    }

    // ---- harness -------------------------------------------------------------------------------

    private static JsonElement Invoke(params string[] args)
    {
        JsonElement envelope = Envelope(args);
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.True,
            $"Expected success, got: {envelope}");
        return envelope.GetProperty("data");
    }

    private static JsonElement InvokeExpectingError(params string[] args)
    {
        JsonElement envelope = Envelope(args);
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.False,
            $"Expected a failure, got: {envelope}");
        return envelope.GetProperty("error");
    }

    private static JsonElement Envelope(string[] args)
    {
        StringWriter captured = new();
        TextWriter previous = Console.Out;
        Console.SetOut(captured);
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.SerializeToElement(
                new { ok = true, data }, Program.JsonOptions);
        }
        catch (Program.UsageException usage)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = "usage", message = usage.Message } },
                Program.JsonOptions);
        }
        catch (Exception exception)
        {
            return JsonSerializer.SerializeToElement(
                new
                {
                    ok = false,
                    error = new
                    {
                        type = Program.ClassifyError(exception),
                        message = Program.Unwrap(exception).Message,
                    },
                },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetOut(previous);
        }
    }
}
