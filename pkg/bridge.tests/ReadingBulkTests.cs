using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using System.Text.RegularExpressions;

namespace MzLibBridge.Tests;

/// <summary>
/// The multi-file form of every reader verb, the per-file block, <c>absent_fields</c>, the
/// <c>source</c> block, mzIdentML's confidence fields, and the three mzLib 1.0.592 quantification
/// views (bridge <c>design/BULK.md</c>).
/// </summary>
/// <remarks>
/// <para>
/// The property these tests exist for above all is that <b>the output never depends on
/// <c>--threads</c></b>: every bulk verb is run at 1 and at 4 threads, and the two envelopes, and
/// the two <c>--out</c> tables, must be byte-identical. That is what lets the default be 1 without
/// being a correctness choice, and it is the property a naive <c>Parallel.ForEach</c> breaks first.
/// </para>
/// <para>
/// Real mzLib fixtures from the pinned worktree, as in <see cref="ReadingCoverageTests"/>; each test
/// that needs one is ignored rather than failed when the worktree lacks it, and
/// <see cref="ReadingCoverageTests.FixtureRoot_Exists_SoTheCoverageSuiteCannotSilentlyVanish"/>
/// fails if the whole tree is missing.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReadingBulkTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-bulk-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void RemoveTempDirectory()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    // ---- fixtures ------------------------------------------------------------------------------

    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    private static string TestRoot() => Path.Combine(RepoRoot(), "code", "mzLib", "mzLib", "Test");

    private static string External(string name) => Fixture(Path.Combine("FileReadingTests", "ExternalFileTypes", name));

    private static string Data(string name) => Fixture(Path.Combine("DataFiles", name));

    private static string Fixture(string relative)
    {
        string path = Path.Combine(TestRoot(), relative);
        if (!File.Exists(path) && !Directory.Exists(path))
            Assert.Ignore($"mzLib fixture not present in the worktree: {path}");
        return path;
    }

    private static string ProteinGroups => External("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv");

    private static string Peptides => External("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv");

    private static string Mzid => Data("PXD078927_msgf_1_1_0.mzid");

    /// <summary>A copy of a fixture under another name, so a batch can hold it twice without being a duplicate.</summary>
    private string CopyOf(string path, string name)
    {
        string copy = Path.Combine(_tempDirectory, name);
        File.Copy(path, copy, overwrite: true);
        return copy;
    }

    /// <summary>For each bulk verb, a batch of real inputs of more than one file (and format, where the view has several).</summary>
    private IEnumerable<(string Verb, string[] Paths, string[] Extra)> Batches() =>
    [
        ("read-results", [Fixture("FileReadingTests/SearchResults/ExcelEditedPeptide.psmtsv"),
            External("FraggerPsm_FragPipev21.1_psm.tsv"), External("DiaNn_LongFormat_report.tsv")], []),
        ("read-records", [Mzid, Data("PXD078927_msgf_1_1_0.mzid.gz"), Data("PXD000783_scaffold_mods_1_1_0.mzid")], []),
        ("read-features", [External("Ms1Feature_TopFDv1.6.2_ms1.feature"),
            External("Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature")], []),
        ("read-matches", [External("MsPathFinderT_AllResults_IcTda.tsv"), External("MsPathFinderT_TargetResults_IcTarget.tsv"),
            External("Casanovo_5.0.0.mztab"), Mzid], ["--scores"]),
        ("read-spectra", [Data("sliced_ethcd.mzML"), Data("withZeros.mgf"), Data("ScanDescriptionTestData.raw"),
            External("Ms2Align_FlashDeconvOpenMs3.0.0_ms2.msalign")], ["--peaks"]),
        ("read-protein-groups", [ProteinGroups, CopyOf(ProteinGroups, "copy_AllQuantifiedProteinGroups.tsv")], []),
        ("read-quantified-peptides", [Peptides, CopyOf(Peptides, "copy_AllQuantifiedPeptides.tsv")], []),
        ("read-occupancy", [ProteinGroups, CopyOf(ProteinGroups, "other_AllQuantifiedProteinGroups.tsv")], []),
    ];

    private static IEnumerable<string> BulkVerbs() =>
    [
        "read-results", "read-records", "read-features", "read-matches", "read-spectra",
        "read-protein-groups", "read-quantified-peptides", "read-occupancy",
    ];

    private (string[] Paths, string[] Extra) BatchFor(string verb) =>
        Batches().Where(b => b.Verb == verb).Select(b => (b.Paths, b.Extra)).Single();

    // ---- the thread-count property -------------------------------------------------------------

    [Test]
    [TestCaseSource(nameof(BulkVerbs))]
    public void OutputIsByteIdenticalAtOneAndFourThreads(string verb)
    {
        (string[] paths, string[] extra) = BatchFor(verb);

        string one = Raw(Stdin(paths), ["readers", verb, "--paths-stdin", "--threads", "1", .. extra]);
        string four = Raw(Stdin(paths), ["readers", verb, "--paths-stdin", "--threads", "4", .. extra]);

        Assert.That(four, Is.EqualTo(one),
            $"{verb}: the envelope changed with the thread count. Files must be emitted in input order " +
            "whatever order they finish in.");
        Assert.That(JsonDocument.Parse(one).RootElement.GetProperty("ok").GetBoolean(), Is.True, one);
    }

    [Test]
    [TestCaseSource(nameof(BulkVerbs))]
    public void WrittenTableIsByteIdenticalAtOneAndFourThreads(string verb)
    {
        (string[] paths, string[] extra) = BatchFor(verb);
        string outOne = Path.Combine(_tempDirectory, "one.tsv");
        string outFour = Path.Combine(_tempDirectory, "four.tsv");

        JsonElement one = Invoke(Stdin(paths), ["readers", verb, "--paths-stdin", "--threads", "1", "--out", outOne, .. extra]);
        Invoke(Stdin(paths), ["readers", verb, "--paths-stdin", "--threads", "4", "--out", outFour, .. extra]);

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllBytes(outFour), Is.EqualTo(File.ReadAllBytes(outOne)));
            Assert.That(one.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null));
            string[] header = File.ReadLines(outOne).First().Split('\t');
            Assert.That(header.Take(2), Is.EqualTo(new[] { "source_index", "source_path" }));
            Assert.That(header, Is.EqualTo(one.GetProperty("column_names").EnumerateArray().Select(c => c.GetString())));
            Assert.That(File.ReadLines(outOne).Count() - 1,
                Is.EqualTo(one.GetProperty("output").GetProperty("row_count").GetInt32()));
            Assert.That(one.GetProperty("row_count").GetInt32(), Is.Zero, "rows went to disk");
        });
    }

    [Test]
    [TestCaseSource(nameof(BulkVerbs))]
    public void EveryInputHasAFilesEntryWithTheSameKeysAsASingleFileRead(string verb)
    {
        (string[] paths, string[] extra) = BatchFor(verb);

        JsonElement bulk = Invoke(Stdin(paths), ["readers", verb, "--paths-stdin", .. extra]);
        JsonElement single = Invoke(null, ["readers", verb, "--path", paths[0], .. extra]);
        JsonElement[] files = bulk.GetProperty("files").EnumerateArray().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(files, Has.Length.EqualTo(paths.Length));
            for (int i = 0; i < paths.Length; i++)
            {
                Assert.That(files[i].GetProperty("path").GetString(), Is.EqualTo(Path.GetFullPath(paths[i])),
                    "files[] is in input order");
                Assert.That(files[i].GetProperty("error").ValueKind, Is.EqualTo(JsonValueKind.Null));
            }

            // BULK.md section 3: every per-file fact a single read carries at the top level lives in
            // files[i] under the same name, so a binding reuses one type for both.
            var fileKeys = files[0].EnumerateObject().Select(p => p.Name).ToHashSet();
            var singleKeys = single.EnumerateObject().Select(p => p.Name).ToHashSet();
            Assert.That(fileKeys.Except(singleKeys), Is.Empty, "a per-file key missing from the single-file answer");

            Assert.That(bulk.GetProperty("record_count").GetInt64(),
                Is.EqualTo(files.Sum(f => f.GetProperty("record_count").GetInt64())));
            Assert.That(bulk.GetProperty("row_count").GetInt32(),
                Is.EqualTo(bulk.GetProperty("columns").GetProperty("source_index").GetArrayLength()));
            Assert.That(bulk.GetProperty("file_count").GetInt32(), Is.EqualTo(paths.Length));
            Assert.That(bulk.GetProperty("read_count").GetInt32(), Is.EqualTo(paths.Length));
            Assert.That(bulk.GetProperty("failed_count").GetInt32(), Is.Zero);
        });
    }

    [Test]
    public void SourceColumnsNumberEachRowByItsInputAndKeepEachFilesOwnOrder()
    {
        (string[] paths, _) = BatchFor("read-spectra");

        JsonElement bulk = Invoke(Stdin(paths), ["readers", "read-spectra", "--paths-stdin"]);
        JsonElement columns = bulk.GetProperty("columns");
        int[] sources = columns.GetProperty("source_index").EnumerateArray().Select(v => v.GetInt32()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(sources, Is.Ordered, "rows are grouped by input, in input order");
            for (int i = 0; i < paths.Length; i++)
            {
                JsonElement single = Invoke(null, ["readers", "read-spectra", "--path", paths[i]]);
                int[] expected = single.GetProperty("columns").GetProperty("one_based_scan_number")
                    .EnumerateArray().Select(v => v.GetInt32()).ToArray();
                int[] actual = columns.GetProperty("one_based_scan_number").EnumerateArray()
                    .Where((_, row) => sources[row] == i).Select(v => v.GetInt32()).ToArray();
                Assert.That(actual, Is.EqualTo(expected), $"input {i}: its rows, in its own order");
            }
        });
    }

    // ---- failures ------------------------------------------------------------------------------

    [Test]
    public void OnErrorFail_StopsAtTheFirstBadInput_AndNamesIt()
    {
        string missing = Path.Combine(_tempDirectory, "missing.mzML");
        JsonElement error = InvokeExpectingError(Stdin([Data("sliced_ethcd.mzML"), missing]),
            ["readers", "read-spectra", "--paths-stdin"]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"),
                "a missing input is still the caller's mistake, in a batch as alone");
            Assert.That(error.GetProperty("message").GetString(), Does.StartWith("Input 1 (").And.Contains("missing.mzML"));
        });
    }

    [Test]
    public void OnErrorFail_KeepsMzLibsOwnTypeForAParseFailure()
    {
        string broken = Path.Combine(_tempDirectory, "broken.mzid");
        File.WriteAllText(broken, "<not mzIdentML");

        JsonElement error = InvokeExpectingError(Stdin([Mzid, broken]), ["readers", "read-records", "--paths-stdin"]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.Not.EqualTo("usage").And.Not.EqualTo("SourceFailure"));
            Assert.That(error.GetProperty("message").GetString(), Does.StartWith("Input 1 ("));
        });
    }

    [Test]
    public void OnErrorSkip_RecordsTheFailureInItsEntry_AndReadsTheRest()
    {
        string missing = Path.Combine(_tempDirectory, "missing.mzML");
        string output = Path.Combine(_tempDirectory, "batch.tsv");

        JsonElement bulk = Invoke(Stdin([missing, Data("sliced_ethcd.mzML")]),
            ["readers", "read-spectra", "--paths-stdin", "--on-error", "skip", "--out", output]);
        JsonElement[] files = bulk.GetProperty("files").EnumerateArray().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(bulk.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(bulk.GetProperty("read_count").GetInt32(), Is.EqualTo(1));
            Assert.That(bulk.GetProperty("on_error").GetString(), Is.EqualTo("skip"));
            JsonElement error = files[0].GetProperty("error");
            Assert.That(error.GetProperty("kind").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("missing.mzML"));
            Assert.That(files[0].GetProperty("record_count").ValueKind, Is.EqualTo(JsonValueKind.Null),
                "a file that was not read has no count, not a count of zero");
            Assert.That(files[0].EnumerateObject().Select(p => p.Name),
                Is.EqualTo(files[1].EnumerateObject().Select(p => p.Name)),
                "a failed entry carries the same keys as a read one, so one type describes both");
            Assert.That(File.ReadLines(output).Skip(1).All(line => line.StartsWith("1\t")), Is.True,
                "only the readable input has rows");
        });
    }

    [Test]
    public void AFailedBatchDoesNotLeaveAHalfWrittenTable()
    {
        string output = Path.Combine(_tempDirectory, "partial.tsv");
        InvokeExpectingError(Stdin([Data("sliced_ethcd.mzML"), Path.Combine(_tempDirectory, "missing.mzML")]),
            ["readers", "read-spectra", "--paths-stdin", "--out", output]);

        Assert.That(File.Exists(output), Is.False, "a table that stopped part-way must not look complete");
    }

    [Test]
    public void AnInputThatIsNotThisViewIsAUsageFailureOfThatInput()
    {
        JsonElement bulk = Invoke(Stdin([External("Casanovo_5.0.0.mztab"), Mzid]),
            ["readers", "read-features", "--paths-stdin", "--on-error", "skip"]);

        Assert.That(bulk.GetProperty("files").EnumerateArray().Select(f => f.GetProperty("error").GetProperty("kind").GetString()),
            Is.All.EqualTo("usage"));
    }

    // ---- usage ---------------------------------------------------------------------------------

    private static IEnumerable<TestCaseData> UsageCases()
    {
        yield return new TestCaseData(new[] { "--path", "x.mzML" }, "mutually exclusive").SetName("path and paths-stdin");
        yield return new TestCaseData(new[] { "--limit", "3" }, "window ONE file").SetName("limit with paths-stdin");
        yield return new TestCaseData(new[] { "--offset", "3" }, "window ONE file").SetName("offset with paths-stdin");
        yield return new TestCaseData(new[] { "--threads", "0" }, "--threads must be").SetName("threads zero");
        yield return new TestCaseData(new[] { "--threads", "-2" }, "--threads must be").SetName("threads below -1");
        yield return new TestCaseData(new[] { "--threads", "many" }, "must be an integer").SetName("threads not a number");
        yield return new TestCaseData(new[] { "--on-error", "ignore" }, "'fail' or 'skip'").SetName("unknown on-error");
        yield return new TestCaseData(new[] { "--threads" }, "has no value").SetName("threads without a value");
    }

    [Test]
    [TestCaseSource(nameof(UsageCases))]
    public void BadBatchOptionsAreUsageErrors(string[] extra, string expected)
    {
        JsonElement error = InvokeExpectingError(Stdin([Data("sliced_ethcd.mzML")]),
            ["readers", "read-spectra", "--paths-stdin", .. extra]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    [Test]
    public void AnEmptyListIsAUsageError()
    {
        JsonElement error = InvokeExpectingError("\n  \n", ["readers", "read-spectra", "--paths-stdin"]);
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("no paths arrived"));
    }

    [Test]
    public void ARepeatedPathIsRefused_RatherThanCountedTwice()
    {
        string mzml = Data("sliced_ethcd.mzML");
        JsonElement error = InvokeExpectingError(Stdin([mzml, "  " + mzml + "  "]), ["readers", "read-spectra", "--paths-stdin"]);
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("repeats input 0"));
    }

    [Test]
    public void OutMayNotOverwriteAnyInput()
    {
        string copy = CopyOf(Data("sliced_ethcd.mzML"), "run.mzML");
        JsonElement error = InvokeExpectingError(Stdin([Data("withZeros.mgf"), copy]),
            ["readers", "read-spectra", "--paths-stdin", "--out", copy]);
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("overwrite input 1"));
    }

    [Test]
    public void OnErrorSkipWithOnePathIsRefused_BecauseItWouldHideTheOnlyFailure()
    {
        JsonElement error = InvokeExpectingError(null,
            ["readers", "read-spectra", "--path", Data("sliced_ethcd.mzML"), "--on-error", "skip"]);
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("needs --paths-stdin"));
    }

    [Test]
    public void ThreadsWithOnePathIsAcceptedAndChangesNothing()
    {
        string mzml = Data("sliced_ethcd.mzML");
        Assert.That(Raw(null, ["readers", "read-spectra", "--path", mzml, "--threads", "4"]),
            Is.EqualTo(Raw(null, ["readers", "read-spectra", "--path", mzml])));
    }

    [Test]
    public void ReadRecordsRefusesABatchOfMixedRecordTypes_BeforeReadingAnything()
    {
        JsonElement error = InvokeExpectingError(Stdin([Mzid, ProteinGroups, Data("PXD078927_msgf_1_1_0.mzid.gz")]),
            ["readers", "read-records", "--paths-stdin"]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(),
                Does.Contain("MzIdentMLRecord (inputs 0, 2)").And.Contain("ProteinGroupFromTsv (inputs 1)"));
        });
    }

    [Test]
    public void ReadRecordsSkipsAnUnidentifiableInputWhenAsked()
    {
        JsonElement bulk = Invoke(Stdin([Path.Combine(_tempDirectory, "gone.mzid"), Mzid]),
            ["readers", "read-records", "--paths-stdin", "--on-error", "skip"]);

        Assert.Multiple(() =>
        {
            Assert.That(bulk.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(bulk.GetProperty("column_names").EnumerateArray().Select(c => c.GetString()),
                Does.Contain("q_value"), "the batch takes its columns from the inputs that can be read");
        });
    }

    // ---- identify ------------------------------------------------------------------------------

    [Test]
    public void IdentifyTakesAListAndAnswersEachInTurn()
    {
        string missing = Path.Combine(_tempDirectory, "nope.psmtsv");
        string[] paths = [Mzid, ProteinGroups, missing, Data("sliced_ethcd.mzML")];

        string one = Raw(Stdin(paths), ["readers", "identify", "--paths-stdin", "--on-error", "skip"]);
        string four = Raw(Stdin(paths), ["readers", "identify", "--paths-stdin", "--on-error", "skip", "--threads", "4"]);
        JsonElement data = JsonDocument.Parse(one).RootElement.GetProperty("data");
        JsonElement[] files = data.GetProperty("files").EnumerateArray().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(four, Is.EqualTo(one));
            Assert.That(files.Select(f => f.GetProperty("file_type").GetString()),
                Is.EqualTo(new[] { "MzIdentML", "MetaMorpheusQuantifiedProteinGroups", null, "MzML" }));
            Assert.That(files[2].GetProperty("error").GetProperty("kind").GetString(), Is.EqualTo("usage"));
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(files[0].EnumerateObject().Select(p => p.Name),
                Is.EquivalentTo(Invoke(null, ["readers", "identify", "--path", Mzid]).EnumerateObject().Select(p => p.Name)));
        });
    }

    [Test]
    public void IdentifyHasNoTableToWrite()
    {
        JsonElement error = InvokeExpectingError(Stdin([Mzid]),
            ["readers", "identify", "--paths-stdin", "--out", Path.Combine(_tempDirectory, "x.tsv")]);
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("no table"));
    }

    [Test]
    public void IdentifyStopsOnAMissingInputByDefault()
    {
        JsonElement error = InvokeExpectingError(Stdin([Mzid, Path.Combine(_tempDirectory, "nope.mzML")]),
            ["readers", "identify", "--paths-stdin"]);
        Assert.That(error.GetProperty("message").GetString(), Does.StartWith("Input 1 ("));
    }

    // ---- version -------------------------------------------------------------------------------

    [Test]
    public void VersionListsExactlyTheVerbsTheSwitchRoutes()
    {
        // Read the switch the way bridge/tools/check_verbs.py does, so the generated list, the
        // bridge's own dispatch and the spec checker all agree on one set.
        string program = File.ReadAllText(Path.Combine(RepoRoot(), "pkg", "bridge", "Program.cs"));
        string body = program[program.IndexOf("return arguments.Verb switch", StringComparison.Ordinal)..];
        body = body[..body.IndexOf("_ =>", StringComparison.Ordinal)];
        string[] switchKeys = Regex.Matches(body, "^\\s*\"([^\"]+)\"\\s*=>", RegexOptions.Multiline)
            .Select(m => m.Groups[1].Value).ToArray();

        string[] listed = Invoke(null, ["version"]).GetProperty("verbs").EnumerateArray().Select(v => v.GetString()!).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(listed, Is.EqualTo(switchKeys));
            Assert.That(listed, Does.Contain("readers read-protein-groups"));
            Assert.That(listed, Does.Contain("readers read-quantified-peptides"));
            Assert.That(listed, Does.Contain("readers read-occupancy"));
        });
    }

    [Test]
    public void AnUnknownCommandListsEveryVerb()
    {
        JsonElement error = InvokeExpectingError(null, ["readers", "read-everything"]);
        string message = error.GetProperty("message").GetString()!;

        Assert.That(Program.Verbs.All(verb => message.Contains(verb, StringComparison.Ordinal)), Is.True, message);
    }

    // ---- absent_fields -------------------------------------------------------------------------

    [Test]
    public void AFlashLfqPeaksTableWithoutAnMbrScoreColumnNamesItAbsent()
    {
        // mzLib #1345: current FlashLFQ tables have no MBR Score column, and mbr_score reads null
        // for every row. Null alone cannot say "the file has no such column"; absent_fields does.
        JsonElement current = Invoke(null, ["readers", "read-records", "--path", External("FlashLFQ_MzLib1.0.591_QuantifiedPeaks.tsv")]);
        JsonElement older = Invoke(null, ["readers", "read-records", "--path", External("FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv")]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(current, "absent_fields"), Is.EqualTo(new[] { "mbr_score" }));
            Assert.That(current.GetProperty("columns").GetProperty("mbr_score").EnumerateArray()
                .Select(v => v.ValueKind), Is.All.EqualTo(JsonValueKind.Null));
            Assert.That(Strings(older, "absent_fields"), Is.EqualTo(new[]
                { "organism", "peak_fwhm", "peak_fwhm_status", "pip_q_value", "pip_pep", "decoy_peptide", "random_rt" }),
                "the seven columns #1345 added are absent from a table written before them");
        });
    }

    [Test]
    public void AnMsPathFinderTargetsFileHasNoQValue_SoItsFabricatedZeroIsNulled()
    {
        // mzLib types MsPathFinderTResult.QValue as a non-nullable double, so a file without the
        // column reads 0 - a perfect q-value for every match.
        JsonElement targets = Invoke(null, ["readers", "read-matches", "--path", External("MsPathFinderT_TargetResults_IcTarget.tsv")]);
        JsonElement tda = Invoke(null, ["readers", "read-matches", "--path", External("MsPathFinderT_AllResults_IcTda.tsv")]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(targets, "absent_fields"), Is.EqualTo(new[] { "q_value", "rank", "pass_threshold" }));
            Assert.That(targets.GetProperty("columns").GetProperty("q_value").EnumerateArray().Select(v => v.ValueKind),
                Is.All.EqualTo(JsonValueKind.Null));
            Assert.That(Strings(tda, "absent_fields"), Does.Not.Contain("q_value"));
            Assert.That(tda.GetProperty("columns").GetProperty("q_value").EnumerateArray().Select(v => v.ValueKind),
                Is.All.EqualTo(JsonValueKind.Number));
        });
    }

    [Test]
    public void MsFraggerIsDecoyIsAbsent_AndPsmtsvIsNot()
    {
        JsonElement fragger = Invoke(null, ["readers", "read-results", "--path", External("FraggerPsm_FragPipev21.1_psm.tsv")]);
        JsonElement psmtsv = Invoke(null, ["readers", "read-results", "--path", Fixture("FileReadingTests/SearchResults/ExcelEditedPeptide.psmtsv")]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(fragger, "absent_fields"), Is.EqualTo(new[] { "is_decoy" }));
            Assert.That(Strings(psmtsv, "absent_fields"), Is.Empty);
        });
    }

    [Test]
    public void AFlashDeconvFeatureFileNamesIntensityAbsent()
    {
        JsonElement flash = Invoke(null, ["readers", "read-features", "--path", External("Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature")]);
        JsonElement topfd = Invoke(null, ["readers", "read-features", "--path", External("Ms1Feature_TopFDv1.6.2_ms1.feature")]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(flash, "absent_fields"), Is.EqualTo(new[] { "intensity", "number_of_isotopes" }));
            Assert.That(Strings(topfd, "absent_fields"), Is.EqualTo(new[] { "number_of_isotopes" }));
        });
    }

    [Test]
    public void ReadRecordsPointsEachExcludedDictionaryAtTheVerbThatCarriesIt()
    {
        JsonElement groups = Invoke(null, ["readers", "read-records", "--path", ProteinGroups, "--limit", "1"]);
        JsonElement mzid = Invoke(null, ["readers", "read-records", "--path", Mzid, "--limit", "1"]);

        string? VerbFor(JsonElement data, string field) => data.GetProperty("excluded_fields").EnumerateArray()
            .Single(e => e.GetProperty("field").GetString() == field).GetProperty("verb").GetString();

        Assert.Multiple(() =>
        {
            Assert.That(VerbFor(groups, "sample_groups"), Is.EqualTo("readers read-protein-groups"));
            Assert.That(VerbFor(mzid, "scores"), Is.EqualTo("readers read-matches"));
            Assert.That(VerbFor(mzid, "match"), Is.Null);
        });
    }

    // ---- read-spectra source -------------------------------------------------------------------

    [Test]
    public void SpectraReportTheInstrumentAndWhenAcquisitionStarted()
    {
        JsonElement mzml = Invoke(null, ["readers", "read-spectra", "--path", Data("sliced_ethcd.mzML"), "--limit", "0"])
            .GetProperty("source");
        JsonElement raw = Invoke(null, ["readers", "read-spectra", "--path", Data("ScanDescriptionTestData.raw"), "--limit", "0"])
            .GetProperty("source");
        JsonElement mgf = Invoke(null, ["readers", "read-spectra", "--path", Data("withZeros.mgf"), "--limit", "0"])
            .GetProperty("source");

        Assert.Multiple(() =>
        {
            Assert.That(mzml.GetProperty("instrument_model").GetString(), Is.EqualTo("Orbitrap Fusion"));
            Assert.That(mzml.GetProperty("instrument_model_accession").GetString(), Is.EqualTo("MS:1002416"));
            Assert.That(mzml.GetProperty("instrument_serial_number").GetString(), Is.EqualTo("FSN10189"));
            Assert.That(mzml.GetProperty("acquisition_start_time").GetString(), Is.EqualTo("2021-03-16T17:09:07Z"));
            Assert.That(mzml.GetProperty("acquisition_start_time_is_utc").GetBoolean(), Is.True);

            // Thermo records the model by name only, and its clock is the instrument PC's.
            Assert.That(raw.GetProperty("instrument_model").GetString(), Is.EqualTo("Orbitrap Fusion Lumos"));
            Assert.That(raw.GetProperty("instrument_model_accession").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(raw.GetProperty("acquisition_start_time").GetString(), Does.Not.EndWith("Z"));
            Assert.That(raw.GetProperty("acquisition_start_time_is_utc").GetBoolean(), Is.False);

            Assert.That(mgf.GetProperty("instrument_model").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(mgf.GetProperty("acquisition_start_time").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    // ---- mzIdentML confidence ------------------------------------------------------------------

    [Test]
    public void MzIdentMLCarriesItsQValueRankAndThreshold()
    {
        JsonElement data = Invoke(null, ["readers", "read-matches", "--path", Mzid]);
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data, "absent_fields"), Is.EqualTo(new[] { "is_decoy" }));
            Assert.That(columns.GetProperty("q_value").EnumerateArray().Select(v => v.ValueKind), Is.All.EqualTo(JsonValueKind.Number));
            Assert.That(columns.GetProperty("rank").EnumerateArray().Select(v => v.GetInt32()), Is.All.GreaterThanOrEqualTo(1));
            Assert.That(data.GetProperty("skipped_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("scores_included").GetBoolean(), Is.False);
        });
    }

    [Test]
    public void ScoresMakeTheMatchTableLong_OneRowPerMatchAndScore()
    {
        JsonElement narrow = Invoke(null, ["readers", "read-matches", "--path", Mzid, "--limit", "2"]);
        JsonElement longForm = Invoke(null, ["readers", "read-matches", "--path", Mzid, "--limit", "2", "--scores"]);
        JsonElement columns = longForm.GetProperty("columns");
        string[] names = columns.GetProperty("score_name").EnumerateArray().Select(v => v.GetString()!).ToArray();
        int[] matchIndex = columns.GetProperty("match_index").EnumerateArray().Select(v => v.GetInt32()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(narrow.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
            Assert.That(narrow.GetProperty("row_count").GetInt32(), Is.EqualTo(2));
            Assert.That(longForm.GetProperty("returned_count").GetInt32(), Is.EqualTo(2),
                "returned_count counts MATCHES, the unit limit and offset count in");
            Assert.That(longForm.GetProperty("row_count").GetInt32(), Is.EqualTo(names.Length).And.GreaterThan(2));
            Assert.That(names, Does.Contain("MS-GF:SpecEValue"));
            Assert.That(matchIndex.Distinct(), Is.EqualTo(new[] { 0, 1 }));
            Assert.That(longForm.GetProperty("truncated").GetBoolean(), Is.True, "limit still counts matches");
        });
    }

    [Test]
    public void ScoresOnAFormatWithoutThemKeepOneRowPerMatch_AndSayWhy()
    {
        JsonElement data = Invoke(null, ["readers", "read-matches", "--path", External("Casanovo_5.0.0.mztab"), "--scores"]);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(data.GetProperty("record_count").GetInt32()));
            Assert.That(Strings(data, "absent_fields"), Does.Contain("score_name").And.Contain("score_value"));
        });
    }

    [Test]
    public void ACrosslinkSearchReportsWhatWasSkippedAndWhy()
    {
        JsonElement data = Invoke(null, ["readers", "read-matches", "--path", Data("OpenxQuest_example_1_2_0.mzid")]);
        JsonElement[] skipped = data.GetProperty("skipped").EnumerateArray().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("skipped_count").GetInt32(), Is.EqualTo(skipped.Length).And.GreaterThan(0));
            Assert.That(skipped.Select(s => s.GetProperty("reason").GetString()), Is.All.EqualTo("crosslink identification"));
            Assert.That(skipped[0].GetProperty("spectrum_identification_item_id").GetString(), Is.Not.Empty);
        });
    }

    [Test]
    public void AFormatWithoutASkipListReportsNullNotZero()
    {
        JsonElement data = Invoke(null, ["readers", "read-matches", "--path", External("Casanovo_5.0.0.mztab")]);
        Assert.That(data.GetProperty("skipped_count").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    // ---- the #1347 views -----------------------------------------------------------------------

    [Test]
    public void ProteinGroupsAreOneRowPerGroupAndSampleGroup_WithMzLibsValues()
    {
        JsonElement data = Invoke(null, ["readers", "read-protein-groups", "--path", ProteinGroups]);
        JsonElement columns = data.GetProperty("columns");
        string[] labels = Strings(data, "sample_labels");
        int groups = data.GetProperty("record_count").GetInt32();

        // Straight from mzLib, to compare cell for cell.
        var file = new Readers.ProteinGroupFromTsvFile(ProteinGroups);
        file.LoadResults();
        Readers.ProteinGroupFromTsv first = file.Results[0];

        Assert.Multiple(() =>
        {
            Assert.That(labels, Has.Length.EqualTo(18));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(groups));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(groups * labels.Length));
            Assert.That(columns.GetProperty("sample_label").EnumerateArray().Take(labels.Length).Select(v => v.GetString()),
                Is.EqualTo(labels), "samples in header order");
            Assert.That(columns.GetProperty("protein_group_name")[0].GetString(), Is.EqualTo(first.ProteinGroupName));
            for (int i = 0; i < labels.Length; i++)
            {
                Readers.SampleGroupMeasurement sample = first.SampleGroups[labels[i]];
                JsonElement intensity = columns.GetProperty("intensity")[i];
                Assert.That(intensity.ValueKind == JsonValueKind.Null ? (double?)null : intensity.GetDouble(), Is.EqualTo(sample.Intensity));
                JsonElement count = columns.GetProperty("spectral_count")[i];
                Assert.That(count.ValueKind == JsonValueKind.Null ? (int?)null : count.GetInt32(), Is.EqualTo(sample.SpectralCount));
            }

            Assert.That(data.GetProperty("excluded_fields").EnumerateArray().Select(e => e.GetProperty("verb").GetString()),
                Is.All.EqualTo("readers read-occupancy"));
            Assert.That(Strings(data, "absent_fields"), Is.Empty);
        });
    }

    [Test]
    public void QuantifiedPeptidesKeepFlashLfqsZero_AndSayWhatItMeans()
    {
        JsonElement data = Invoke(null, ["readers", "read-quantified-peptides", "--path", Peptides]);
        JsonElement columns = data.GetProperty("columns");
        var zeroes = columns.GetProperty("intensity").EnumerateArray()
            .Select((v, i) => (v, detection: columns.GetProperty("detection_type")[i].GetString()))
            .Where(p => p.v.ValueKind == JsonValueKind.Number && p.v.GetDouble() == 0)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(zeroes, Is.Not.Empty, "the fixture has unquantified peptides");
            Assert.That(zeroes.Select(z => z.detection), Is.All.Not.EqualTo("MSMS"),
                "a written 0 is a peptide FlashLFQ did not quantify");
            Assert.That(Strings(data, "absent_fields"), Is.EqualTo(new[] { "peak_order", "retention_time" }),
                "MetaMorpheus is not IsoTracker: no peak order, no per-sample retention time");
            Assert.That(data.GetProperty("caveats")[0].GetString(), Does.StartWith("intensity 0 is NOT a measured zero"));
        });
    }

    [Test]
    public void OccupancyIsOneRowPerSite_FromBothCells()
    {
        JsonElement data = Invoke(null, ["readers", "read-occupancy", "--path", ProteinGroups]);
        JsonElement columns = data.GetProperty("columns");
        string[] bases = columns.GetProperty("basis").EnumerateArray().Select(v => v.GetString()!).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(bases.Distinct().Order(), Is.EqualTo(new[] { "count", "intensity" }));
            Assert.That(data.GetProperty("truncated_cell_count").GetInt32(), Is.Zero);
            Assert.That(Strings(data, "failed_fields"), Is.Empty);
            int trimethyl = Array.FindIndex(columns.GetProperty("modification").EnumerateArray().ToArray(),
                v => v.GetString() == "N6,N6,N6-trimethyllysine on K");
            Assert.That(trimethyl, Is.GreaterThanOrEqualTo(0), "a name with commas survives the parse whole");
            Assert.That(columns.GetProperty("position")[trimethyl].GetInt32(), Is.EqualTo(52));
        });
    }

    [Test]
    public void AMalformedOccupancyCellIsNamedInFailedFields_AndATruncatedOneIsCounted()
    {
        string original = File.ReadAllText(ProteinGroups);
        const string countCell = "pos52[N6,N6,N6-trimethyllysine on K,info:fraction=1.00(1/1)]";
        const string intensityCell = "pos52[N6,N6,N6-trimethyllysine on K,info:fraction=1.0000(3.559E+06/3.559E+06)]";
        Assume.That(original, Does.Contain(countCell).And.Contain(intensityCell));

        string edited = Path.Combine(_tempDirectory, "edited_AllQuantifiedProteinGroups.tsv");
        File.WriteAllText(edited, original
            .Replace(countCell, "this is not an occupancy cell")
            .Replace(intensityCell, "Output too long for Excel"));

        JsonElement data = Invoke(null, ["readers", "read-occupancy", "--path", edited]);

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data, "failed_fields"), Is.EqualTo(new[] { "count_occupancy: FormatException" }));
            Assert.That(data.GetProperty("truncated_cell_count").GetInt32(), Is.GreaterThan(0));
            Assert.That(data.GetProperty("row_count").GetInt32(),
                Is.LessThan(Invoke(null, ["readers", "read-occupancy", "--path", ProteinGroups]).GetProperty("row_count").GetInt32()),
                "neither cell contributes a row");
        });
    }

    [Test]
    [TestCase("read-protein-groups")]
    [TestCase("read-quantified-peptides")]
    [TestCase("read-occupancy")]
    public void TheQuantViewsRefuseOtherFiles_AndNameWhatTheyRead(string verb)
    {
        JsonElement error = InvokeExpectingError(null, ["readers", verb, "--path", Mzid]);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("'MzIdentML' files cannot be read by " + verb));
        });
    }

    [Test]
    public void AGroupTableWithoutSampleColumnsStillListsEachGroup()
    {
        // Strip every per-sample column: the groups remain, one row each, with the sample fields absent.
        string[] lines = File.ReadAllLines(ProteinGroups);
        string[] header = lines[0].Split('\t');
        int[] keep = Enumerable.Range(0, header.Length)
            .Where(i => !Regex.IsMatch(header[i], "^(SpectralCount|Intensity|CountOccupancy|IntensityOccupancy)_")).ToArray();
        string stripped = Path.Combine(_tempDirectory, "bare_AllQuantifiedProteinGroups.tsv");
        File.WriteAllLines(stripped, lines.Select(line => string.Join('\t', keep.Select(i => line.Split('\t')[i]))));

        JsonElement groups = Invoke(null, ["readers", "read-protein-groups", "--path", stripped]);
        JsonElement occupancy = Invoke(null, ["readers", "read-occupancy", "--path", stripped]);

        Assert.Multiple(() =>
        {
            Assert.That(groups.GetProperty("row_count").GetInt32(), Is.EqualTo(groups.GetProperty("record_count").GetInt32()),
                "with no sample columns each group is one row");
            Assert.That(Strings(groups, "absent_fields"), Is.EqualTo(new[] { "sample_label", "spectral_count", "intensity" }));
            Assert.That(occupancy.GetProperty("row_count").GetInt32(), Is.Zero);
            Assert.That(Strings(occupancy, "absent_fields"), Does.Contain("fraction"));
        });
    }

    // ---- harness -------------------------------------------------------------------------------

    private static string Stdin(IEnumerable<string> paths) => string.Join('\n', paths) + "\n";

    private static string[] Strings(JsonElement data, string name) =>
        data.GetProperty(name).EnumerateArray().Select(v => v.GetString()!).ToArray();

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

    /// <summary>
    /// Runs a verb as Main would and returns the envelope text exactly as it would be written, so a
    /// byte-for-byte comparison is a comparison of what a caller receives.
    /// </summary>
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
