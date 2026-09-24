using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>sdrf validate</c>, <c>lint</c>, <c>assess</c>, <c>samples</c> and <c>parse-age</c>.
/// </summary>
/// <remarks>
/// <para>
/// Most run against three small documents in pyMzLib's own fixture directory rather than mzLib's,
/// because mzLib's three SDRF fixtures carry no ages, no conflicting sample and no skeleton, and
/// those are the cases these verbs exist for. <b>sdrf_cohort</b>: six samples in two fractions each,
/// one age per precision plus a bare number and a mis-cased reserved word, and sample S6 whose two
/// rows disagree about its disease. <b>sdrf_cohort_partner</b>: a second study written with one
/// drift of every kind against the first, and a searched-data-file name that differs from the
/// first's only by case — which must NOT be reported (mzLib #1335). <b>sdrf_skeleton</b>: reserved
/// words and one replicate number, missing most required columns.
/// </para>
/// <para>
/// What is asserted is the projection — shape, order, nulls, bulk behaviour — plus enough of
/// mzLib's answer to prove the right mzLib call was made. The rules themselves are mzLib's and are
/// tested there.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SdrfAnalysisTests
{
    // ---- sdrf validate ------------------------------------------------------------------------

    [Test]
    public void Validate_ReportsMzLibsFindingsAsATable()
    {
        JsonElement data = Invoke("sdrf", "validate", "--path", Fixture("sdrf_cohort.sdrf.tsv"));
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("error_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("warning_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(12));
            Assert.That(Strings(columns.GetProperty("rule")), Is.All.EqualTo("ReservedWordCase"));
            Assert.That(Strings(columns.GetProperty("severity")), Is.All.EqualTo("Warning"));
            // S5's two rows, 0-based 8 and 9, which are lines 10 and 11 of the file.
            Assert.That(Ints(columns.GetProperty("row_index")), Is.EqualTo(new[] { 8, 9 }));
            Assert.That(Ints(columns.GetProperty("line_number")), Is.EqualTo(new[] { 10, 11 }));
        });
    }

    [Test]
    public void Validate_ADocumentLevelFindingHasNullRowAndLine()
    {
        // Missing required columns are about the document, not a row. Null here - not 0, not -1 -
        // is what stops them being read as "line 2".
        JsonElement data = Invoke("sdrf", "validate", "--path", Fixture("sdrf_skeleton.sdrf.tsv"));
        JsonElement columns = data.GetProperty("columns");
        int i = Strings(columns.GetProperty("rule")).ToList().IndexOf("RequiredColumn");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("is_valid").GetBoolean(), Is.False);
            Assert.That(i, Is.GreaterThanOrEqualTo(0));
            Assert.That(columns.GetProperty("row_index")[i].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(columns.GetProperty("line_number")[i].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(columns.GetProperty("column_name")[i].GetString(), Is.Not.Empty);
        });
    }

    [Test]
    public void Validate_Bulk_OneTableInInputOrderWithAFileBlockEach()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_skeleton.sdrf.tsv")}\n{Fixture("sdrf_cohort.sdrf.tsv")}",
            "sdrf", "validate", "--paths-stdin");
        int[] sources = Ints(data.GetProperty("columns").GetProperty("source_index"));
        JsonElement files = data.GetProperty("files");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("file_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("valid_count").GetInt32(), Is.EqualTo(1));
            Assert.That(sources, Is.Ordered, "rows must follow input order");
            Assert.That(sources.Distinct(), Is.EqualTo(new[] { 0, 1 }));
            Assert.That(files[0].GetProperty("is_valid").GetBoolean(), Is.False);
            Assert.That(files[1].GetProperty("is_valid").GetBoolean(), Is.True);
            Assert.That(files[1].GetProperty("error").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(Strings(data.GetProperty("column_names")).Take(2),
                Is.EqualTo(new[] { "source_index", "source_path" }));
        });
    }

    [Test]
    public void Validate_Bulk_OutputDoesNotDependOnThreads()
    {
        // BULK.md §2: byte-identical output at any thread count is a tested property.
        string stdin = string.Join('\n', new[]
        {
            Fixture("sdrf_cohort.sdrf.tsv"), Fixture("sdrf_skeleton.sdrf.tsv"),
            Fixture("sdrf_cohort_partner.sdrf.tsv"), Mzlib("PXD059974.sdrf.tsv"), Mzlib("PXD000070.sdrf.tsv"),
        });

        foreach (string verb in new[] { "validate", "assess", "samples" })
        {
            string one = InvokeWithStdin(stdin, "sdrf", verb, "--paths-stdin", "--threads", "1").GetRawText();
            string four = InvokeWithStdin(stdin, "sdrf", verb, "--paths-stdin", "--threads", "4").GetRawText();
            string all = InvokeWithStdin(stdin, "sdrf", verb, "--paths-stdin", "--threads", "-1").GetRawText();
            Assert.That(four, Is.EqualTo(one), verb);
            Assert.That(all, Is.EqualTo(one), verb);
        }
    }

    [Test]
    public void Bulk_SkipRecordsTheFailureAndCarriesOn()
    {
        string missing = Path.Combine(Path.GetTempPath(), "no-such-file.sdrf.tsv");
        JsonElement data = InvokeWithStdin($"{missing}\n{Fixture("sdrf_cohort.sdrf.tsv")}",
            "sdrf", "validate", "--paths-stdin", "--on-error", "skip");
        JsonElement failed = data.GetProperty("files")[0];

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("read_count").GetInt32(), Is.EqualTo(1));
            Assert.That(failed.GetProperty("error").GetProperty("kind").GetString(), Is.EqualTo("usage"));
            Assert.That(failed.GetProperty("error").GetProperty("message").GetString(), Does.Contain("no-such-file"));
            Assert.That(failed.GetProperty("is_valid").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(Ints(data.GetProperty("columns").GetProperty("source_index")), Is.All.EqualTo(1));
        });
    }

    [Test]
    public void Bulk_SkipRecordsAnUnreadableFileAsACorrectnessFailure()
    {
        // mzLib refuses an empty file with an MzLibException; under skip that is recorded with its
        // .NET type, as a correctness failure rather than a usage one.
        string empty = Path.Combine(Path.GetTempPath(), $"pymzlib-empty-{Guid.NewGuid():N}.sdrf.tsv");
        File.WriteAllText(empty, "");
        try
        {
            JsonElement data = InvokeWithStdin($"{empty}\n{Fixture("sdrf_cohort.sdrf.tsv")}",
                "sdrf", "assess", "--paths-stdin", "--on-error", "skip");
            JsonElement error = data.GetProperty("files")[0].GetProperty("error");

            Assert.Multiple(() =>
            {
                Assert.That(error.GetProperty("kind").GetString(), Is.EqualTo("correctness"));
                Assert.That(error.GetProperty("message").GetString(), Does.StartWith("MzLibException:"));
                Assert.That(data.GetProperty("files")[1].GetProperty("verdict").GetString(), Is.EqualTo("Informative"));
            });
        }
        finally
        {
            File.Delete(empty);
        }
    }

    [Test]
    public void Bulk_FailRaisesTheFirstFailureInInputOrder()
    {
        string empty = Path.Combine(Path.GetTempPath(), $"pymzlib-empty-{Guid.NewGuid():N}.sdrf.tsv");
        File.WriteAllText(empty, "");
        try
        {
            JsonElement error = InvokeWithStdinExpectingError($"{Fixture("sdrf_cohort.sdrf.tsv")}\n{empty}",
                "sdrf", "validate", "--paths-stdin", "--threads", "2");
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("MzLibException"));
        }
        finally
        {
            File.Delete(empty);
        }
    }

    [Test]
    public void Bulk_FailReportsTheMissingFileBeforeReadingAnything()
    {
        JsonElement error = InvokeWithStdinExpectingError(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\n{Path.Combine(Path.GetTempPath(), "absent.sdrf.tsv")}",
            "sdrf", "samples", "--paths-stdin");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("absent.sdrf.tsv"));
        });
    }

    [TestCase(new[] { "--path", "x", "--paths-stdin" }, "not both")]
    [TestCase(new string[0], "--paths-stdin")]
    [TestCase(new[] { "--paths-stdin", "--threads", "0" }, "--threads")]
    [TestCase(new[] { "--paths-stdin", "--threads", "-2" }, "--threads")]
    [TestCase(new[] { "--paths-stdin", "--on-error", "ignore" }, "'fail' or 'skip'")]
    [TestCase(new[] { "--paths-stdin", "--on-error" }, "needs a value")]
    public void Bulk_BadOptionsAreUsageErrors(string[] options, string expected)
    {
        JsonElement error = InvokeWithStdinExpectingError(Fixture("sdrf_cohort.sdrf.tsv"),
            new[] { "sdrf", "validate" }.Concat(options).ToArray());

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    [Test]
    public void Bulk_SkipWithOnePathIsAUsageError()
    {
        JsonElement error = InvokeExpectingError("sdrf", "assess", "--path", Fixture("sdrf_cohort.sdrf.tsv"),
            "--on-error", "skip");
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("applies to --paths-stdin"));
    }

    [Test]
    public void Bulk_ADuplicatePathIsRefused()
    {
        string path = Fixture("sdrf_cohort.sdrf.tsv");
        JsonElement error = InvokeWithStdinExpectingError($"{path}\n{path}", "sdrf", "assess", "--paths-stdin");
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("appears twice"));
    }

    [Test]
    public void Bulk_EmptyStdinIsAUsageError()
    {
        JsonElement error = InvokeWithStdinExpectingError("", "sdrf", "samples", "--paths-stdin");
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("stdin held no paths"));
    }

    [Test]
    public void Single_AMissingPathIsAUsageError()
    {
        JsonElement error = InvokeExpectingError("sdrf", "validate", "--path",
            Path.Combine(Path.GetTempPath(), "gone.sdrf.tsv"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("gone.sdrf.tsv"));
    }

    // ---- sdrf assess --------------------------------------------------------------------------

    [TestCase("sdrf_cohort.sdrf.tsv", "Informative")]
    [TestCase("sdrf_skeleton.sdrf.tsv", "Skeleton")]
    public void Assess_ReturnsMzLibsVerdict(string file, string verdict)
    {
        JsonElement data = Invoke("sdrf", "assess", "--path", Fixture(file));
        Assert.That(data.GetProperty("verdict").GetString(), Is.EqualTo(verdict));
    }

    [Test]
    public void Assess_ShowsTheEvidenceEachCheckWasDecidedFrom()
    {
        JsonElement data = Invoke("sdrf", "assess", "--path", Mzlib("PXD000070.sdrf.tsv"));
        JsonElement columns = data.GetProperty("columns");
        string[] roles = Strings(columns.GetProperty("role"));

        Assert.Multiple(() =>
        {
            // PXD000070 has no factor value column at all, so the first check cannot pass.
            Assert.That(data.GetProperty("verdict").GetString(), Is.EqualTo("Partial"));
            Assert.That(data.GetProperty("factor_value_varies").GetBoolean(), Is.False);
            Assert.That(roles, Has.None.EqualTo("factor_value"));
            Assert.That(roles.Last(), Is.EqualTo("biological_replicate"));
            // Organism is deliberately not a sample characteristic here: a search can fill it.
            Assert.That(Strings(columns.GetProperty("column_name")), Has.None.EqualTo("characteristics[organism]"));
        });
    }

    [Test]
    public void Assess_Bulk_CountsVerdicts()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\n{Fixture("sdrf_skeleton.sdrf.tsv")}\n{Mzlib("PXD000070.sdrf.tsv")}",
            "sdrf", "assess", "--paths-stdin");
        JsonElement counts = data.GetProperty("verdict_counts");

        Assert.Multiple(() =>
        {
            Assert.That(counts.GetProperty("informative").GetInt32(), Is.EqualTo(1));
            Assert.That(counts.GetProperty("partial").GetInt32(), Is.EqualTo(1));
            Assert.That(counts.GetProperty("skeleton").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("files")[1].GetProperty("verdict").GetString(), Is.EqualTo("Skeleton"));
        });
    }

    // ---- sdrf lint ----------------------------------------------------------------------------

    [Test]
    public void Lint_FindsOneDriftOfEachKindTheFixturesPlanted()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\tcohort\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}\tpartner",
            "sdrf", "lint");
        string[] kinds = Strings(data.GetProperty("columns").GetProperty("kind")).Distinct().OrderBy(k => k).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("finding_count").GetInt32(), Is.EqualTo(4));
            Assert.That(kinds, Is.EqualTo(new[]
                { "AccessionNameConflict", "ColumnNameVariant", "MixedTermAndFreeText", "ValueCaseVariant" }));
            Assert.That(Caveats(data), Has.None.Contains("WHERE THE FILES SIT"));
        });
    }

    [Test]
    public void Lint_DoesNotReportSearchedDataFileNamesAsDrift()
    {
        // mzLib #1335: comment[searched data file] is unique per row by construction. The fixtures
        // write S1_F1.mzML and s1_f1.mzML, which would otherwise be a ValueCaseVariant.
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\tcohort\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}\tpartner",
            "sdrf", "lint");

        Assert.That(Strings(data.GetProperty("columns").GetProperty("column_name")),
            Has.None.EqualTo("comment[searched data file]"));
    }

    [Test]
    public void Lint_OneRowPerVariant_MajorityFirst_DocumentsAsAList()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\tcohort\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}\tpartner",
            "sdrf", "lint");
        JsonElement columns = data.GetProperty("columns");
        int[] findings = Ints(columns.GetProperty("finding_index"));
        int[] ranks = Ints(columns.GetProperty("variant_rank"));
        int column = Strings(columns.GetProperty("kind")).ToList().IndexOf("ColumnNameVariant");

        Assert.Multiple(() =>
        {
            Assert.That(findings, Is.EqualTo(new[] { 0, 0, 1, 1, 2, 2, 3, 3 }));
            Assert.That(ranks, Is.EqualTo(new[] { 0, 1, 0, 1, 0, 1, 0, 1 }));
            Assert.That(columns.GetProperty("documents")[0].ValueKind, Is.EqualTo(JsonValueKind.Array));
            Assert.That(columns.GetProperty("column_name")[column].ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void Lint_WithoutLabels_WarnsThatNamesDependOnWhereTheFilesSit()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}", "sdrf", "lint");
        Assert.That(Caveats(data), Has.Some.Contains("WHERE THE FILES SIT"));
    }

    [Test]
    public void Lint_WithNoStdin_IsAUsageError()
    {
        JsonElement error = InvokeWithStdinExpectingError("", "sdrf", "lint");
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("on stdin"));
    }

    // ---- sdrf samples -------------------------------------------------------------------------

    [Test]
    public void Samples_OneSamplePerSourceNameInDocumentOrder()
    {
        JsonElement data = Invoke("sdrf", "samples", "--path", Fixture("sdrf_cohort.sdrf.tsv"));
        string[] names = Strings(data.GetProperty("columns").GetProperty("source_name")).Distinct().ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("sample_count").GetInt32(), Is.EqualTo(6));
            Assert.That(names, Is.EqualTo(new[] { "S1", "S2", "S3", "S4", "S5", "S6" }));
            Assert.That(Ints(data.GetProperty("columns").GetProperty("sample_row_count")), Is.All.EqualTo(2));
            Assert.That(data.GetProperty("problems").GetArrayLength(), Is.Zero);
        });
    }

    [Test]
    public void Samples_AConflictingColumnIsNamedAndWithheld()
    {
        JsonElement data = Invoke("sdrf", "samples", "--path", Fixture("sdrf_cohort.sdrf.tsv"));
        JsonElement c = data.GetProperty("columns");
        int i = Strings(c.GetProperty("status")).ToList().IndexOf("conflicting");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("conflict_count").GetInt32(), Is.EqualTo(1));
            Assert.That(c.GetProperty("source_name")[i].GetString(), Is.EqualTo("S6"));
            Assert.That(c.GetProperty("column_name")[i].GetString(), Is.EqualTo("characteristics[disease]"));
            Assert.That(c.GetProperty("value")[i].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(c.GetProperty("position")[i].ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void Samples_ParsesAgesInYearsOnTheAgeRowsOnly()
    {
        JsonElement data = Invoke("sdrf", "samples", "--path", Fixture("sdrf_cohort.sdrf.tsv"));
        JsonElement c = data.GetProperty("columns");
        string[] columns = Strings(c.GetProperty("column_name"));
        string[] samples = Strings(c.GetProperty("source_name"));
        int Age(string sample) => Enumerable.Range(0, columns.Length)
            .Single(i => columns[i] == "characteristics[age]" && samples[i] == sample);
        int organism = Array.IndexOf(columns, "characteristics[organism]");

        Assert.Multiple(() =>
        {
            Assert.That(c.GetProperty("age_years")[Age("S2")].GetDouble(), Is.EqualTo(62.5));
            Assert.That(c.GetProperty("age_precision")[Age("S2")].GetString(), Is.EqualTo("Range"));
            // >=90Y: no upper bound, which JSON cannot write as infinity.
            Assert.That(c.GetProperty("age_max_years")[Age("S3")].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(c.GetProperty("age_precision")[Age("S3")].GetString(), Is.EqualTo("LowerBound"));
            Assert.That(c.GetProperty("age_refusal")[Age("S4")].GetString(), Is.EqualTo("no_unit"));
            Assert.That(c.GetProperty("age_refusal")[Age("S5")].GetString(), Is.EqualTo("reserved_word"));
            Assert.That(c.GetProperty("age_precision")[Age("S6")].GetString(), Is.EqualTo("Exact"));
            Assert.That(c.GetProperty("age_refusal")[organism].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(c.GetProperty("age_years")[organism].ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void Samples_ARepeatedColumnIsOneRowPerPosition()
    {
        // TmtChannelLevel has no repeats; build one: a sample column written twice.
        string path = Path.Combine(Path.GetTempPath(), $"pymzlib-repeat-{Guid.NewGuid():N}.sdrf.tsv");
        File.WriteAllText(path,
            "source name\tcharacteristics[organism part]\tcharacteristics[organism part]\tassay name\n" +
            "s1\tliver\tkidney\trun 1\n");
        try
        {
            JsonElement c = Invoke("sdrf", "samples", "--path", path).GetProperty("columns");
            string[] columns = Strings(c.GetProperty("column_name"));
            int[] rows = Enumerable.Range(0, columns.Length).Where(i => columns[i] == "characteristics[organism part]").ToArray();

            Assert.That(rows.Select(i => c.GetProperty("position")[i].GetInt32()), Is.EqualTo(new[] { 0, 1 }));
            Assert.That(rows.Select(i => c.GetProperty("value")[i].GetString()), Is.EqualTo(new[] { "liver", "kidney" }));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void Samples_ADocumentWithoutSourceNameReportsAProblemNotAnError()
    {
        string path = Path.Combine(Path.GetTempPath(), $"pymzlib-nosource-{Guid.NewGuid():N}.sdrf.tsv");
        File.WriteAllText(path, "assay name\tcomment[data file]\nrun 1\ta.raw\n");
        try
        {
            JsonElement data = Invoke("sdrf", "samples", "--path", path);
            Assert.That(data.GetProperty("sample_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("problems")[0].GetString(), Does.Contain("no source name column"));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void Samples_Bulk_SumsSamplesAndKeysRowsBySource()
    {
        JsonElement data = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}",
            "sdrf", "samples", "--paths-stdin");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("sample_count").GetInt32(), Is.EqualTo(8));
            Assert.That(data.GetProperty("files")[1].GetProperty("sample_count").GetInt32(), Is.EqualTo(2));
            Assert.That(Ints(data.GetProperty("columns").GetProperty("source_index")), Is.Ordered);
        });
    }

    // ---- sdrf parse-age -----------------------------------------------------------------------

    [Test]
    public void ParseAge_OneRowPerLine_BlankLinesKept()
    {
        JsonElement data = InvokeWithStdin("58Y\n\n>=90Y\n63\n6-8 weeks\n", "sdrf", "parse-age");
        JsonElement c = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("cell_count").GetInt32(), Is.EqualTo(5), "the trailing newline is not a cell");
            Assert.That(data.GetProperty("parsed_count").GetInt32(), Is.EqualTo(3));
            Assert.That(Strings(c.GetProperty("cell")), Is.EqualTo(new[] { "58Y", "", ">=90Y", "63", "6-8 weeks" }));
            Assert.That(c.GetProperty("refusal")[1].GetString(), Is.EqualTo("empty"));
            Assert.That(c.GetProperty("min_years")[2].GetDouble(), Is.EqualTo(90));
            Assert.That(c.GetProperty("max_years")[2].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(c.GetProperty("refusal")[3].GetString(), Is.EqualTo("no_unit"));
            Assert.That(c.GetProperty("follows_specification")[4].GetBoolean(), Is.False);
            Assert.That(c.GetProperty("years")[4].GetDouble(), Is.EqualTo(7d * 7 / 365.25).Within(1e-12));
        });
    }

    [TestCase("about forty", "unreadable")]
    [TestCase("Not Applicable", "reserved_word")]
    [TestCase("  ", "empty")]
    public void ParseAge_NamesWhyACellWasRefused(string cell, string refusal)
    {
        JsonElement c = InvokeWithStdin(cell, "sdrf", "parse-age").GetProperty("columns");
        Assert.Multiple(() =>
        {
            Assert.That(c.GetProperty("refusal")[0].GetString(), Is.EqualTo(refusal));
            Assert.That(c.GetProperty("years")[0].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(c.GetProperty("precision")[0].ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void ParseAge_WithNoStdin_IsAUsageError()
    {
        JsonElement error = InvokeWithStdinExpectingError("", "sdrf", "parse-age");
        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
    }

    // ---- sdrf read / pool with comment[searched data file] ------------------------------------

    [Test]
    public void ReadAndPool_CarryTheSearchedDataFileColumn()
    {
        // #1327 made mzLib's builder write it; read and pool must carry it like any other column.
        JsonElement read = Invoke("sdrf", "read", "--path", Fixture("sdrf_cohort.sdrf.tsv"));
        JsonElement pooled = InvokeWithStdin(
            $"{Fixture("sdrf_cohort.sdrf.tsv")}\ta\n{Fixture("sdrf_cohort_partner.sdrf.tsv")}\tb", "sdrf", "pool");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(read.GetProperty("column_names")), Contains.Item("comment[searched data file]"));
            Assert.That(Strings(pooled.GetProperty("column_names")), Contains.Item("comment[searched data file]"));
        });
    }

    // ---- helpers ------------------------------------------------------------------------------

    /// <summary>A fixture in pyMzLib's own Python test directory, shared with the Python tests.</summary>
    private static string Fixture(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        return Path.Combine(root, "pkg", "python", "tests", "fixtures", name);
    }

    /// <summary>An SDRF fixture from the pinned mzLib worktree, or an ignored test without one.</summary>
    private static string Mzlib(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        string path = Path.Combine(root, "code", "mzLib", "mzLib", "Test", "FileReadingTests",
            "ExternalFileTypes", name);
        if (!File.Exists(path))
            Assert.Ignore($"mzLib SDRF fixture not present in the worktree: {path}");
        return path;
    }

    private static string[] Strings(JsonElement array) =>
        array.EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static int[] Ints(JsonElement array) =>
        array.EnumerateArray().Select(x => x.GetInt32()).ToArray();

    private static string[] Caveats(JsonElement data) => Strings(data.GetProperty("caveats"));

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(null, args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) =>
        Unwrap(Envelope(null, args), expectOk: false);

    private static JsonElement InvokeWithStdin(string stdin, params string[] args) =>
        Unwrap(Envelope(stdin, args), expectOk: true);

    private static JsonElement InvokeWithStdinExpectingError(string stdin, params string[] args) =>
        Unwrap(Envelope(stdin, args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk),
            $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    /// <summary>Dispatches a verb with <paramref name="stdin"/> as the console, shaped as <c>Main</c> would.</summary>
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
            Console.SetIn(previousIn);
        }
    }
}
