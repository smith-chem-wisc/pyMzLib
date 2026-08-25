using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the <c>sdrf</c> verbs, read against mzLib's own SDRF fixtures.
/// </summary>
/// <remarks>
/// <para>
/// The three fixtures are chosen for what each one breaks. <b>PXD000070</b> repeats
/// <c>comment[modification parameters]</c> eight times, which is what rules out a name-keyed
/// payload. <b>PXD059974</b> is ragged — a 46-column header with 17 of its 23 rows carrying 42
/// cells — which is what rules out padding. <b>PXD026824</b> is a second experiment with a
/// different column set, which is what makes the pooled union non-trivial.
/// </para>
/// <para>
/// What is deliberately NOT tested here is whether a document is valid. mzLib owns that in
/// <c>SdrfValidator</c> and <c>SdrfDriftLint</c>, both currently <c>internal</c> to its Readers
/// assembly, so the bridge cannot reach them and must not grow its own second opinion.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SdrfTests
{
    // ---- sdrf read ----------------------------------------------------------------------------

    [Test]
    public void Read_ReturnsTheHeaderAndOneCellListPerRow()
    {
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv"));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("column_names").GetArrayLength(), Is.EqualTo(31));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(6));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(6));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("rows").GetArrayLength(), Is.EqualTo(6));
        });
    }

    [Test]
    public void Read_KeepsARepeatedColumnNameRatherThanCollapsingIt()
    {
        // The reason the payload is row-major. A name-to-values map would keep one of these eight
        // and silently drop the rest, and nothing downstream could tell that had happened.
        string[] names = Columns(Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv")));

        Assert.That(names.Count(n => n == "comment[modification parameters]"), Is.EqualTo(8));
    }

    [Test]
    public void Read_LeavesAControlledVocabularyCellWhole()
    {
        // An SDRF CV term is itself semicolon-delimited. read-records joins an SdrfRow's cells with
        // ";" and so cannot be split back apart; this verb exists because that is unrecoverable.
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv"));
        string[] names = Columns(data);
        int i = Array.IndexOf(names, "comment[modification parameters]");

        string cell = data.GetProperty("rows")[0][i].GetString()!;

        Assert.That(cell, Is.EqualTo("NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed"));
    }

    [Test]
    public void Read_PreservesARaggedRowRatherThanPaddingIt()
    {
        // mzLib preserves raggedness so the file round-trips byte for byte. Padding here would
        // invent cells that the document does not contain, inside a payload whose whole promise is
        // that cells are what was written.
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD059974.sdrf.tsv"));
        int width = data.GetProperty("column_names").GetArrayLength();

        int[] lengths = data.GetProperty("rows").EnumerateArray()
            .Select(r => r.GetArrayLength()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(width, Is.EqualTo(46));
            Assert.That(lengths, Has.Some.LessThan(width), "the ragged rows were padded");
            Assert.That(lengths.Count(l => l < width), Is.EqualTo(17));
        });
    }

    [Test]
    public void Read_CountsTheRaggedRowsOverTheWholeDocumentNotTheWindow()
    {
        // The caveat describes the DOCUMENT, so it must not change with --limit. A caveat that
        // silently narrowed to the visible window would report a clean file to anyone who paged.
        JsonElement windowed = Invoke("sdrf", "read", "--path", Sdrf("PXD059974.sdrf.tsv"), "--limit", "1");

        Assert.That(Caveats(windowed), Has.Some.Contains("17 of 22 rows are SHORT"));
    }

    [Test]
    public void Read_DoesNotClaimRaggednessForAFileThatHasNone()
    {
        // The corollary, and the one that makes the caveat worth reading: it is absent when the
        // document is rectangular, rather than being boilerplate on every response.
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv"));

        Assert.That(Caveats(data), Has.None.Contains("SHORT"));
    }

    [Test]
    public void Read_WindowsRowsAndSaysWhenItLeftSomeBehind()
    {
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv"),
            "--offset", "2", "--limit", "2");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(6));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("offset").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
        });
    }

    [Test]
    public void Read_AnOffsetPastTheEndIsNotACompleteAnswer()
    {
        // `offset + count < row_count` would call this complete. It is not: every row was left out.
        JsonElement data = Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv"), "--offset", "99");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
        });
    }

    [Test]
    public void Read_WithoutAPath_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError("sdrf", "read");

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
    }

    // ---- sdrf pool ----------------------------------------------------------------------------

    [Test]
    public void Pool_MergesEveryRowAndStampsProvenance()
    {
        JsonElement data = InvokeWithStdin(
            $"{Sdrf("PXD000070.sdrf.tsv")}\tmalaria\n{Sdrf("PXD026824.sdrf.tsv")}\tcolon",
            "sdrf", "pool");

        string[] names = Columns(data);
        int source = Array.IndexOf(names, "comment[source document]");
        string[] provenance = data.GetProperty("rows").EnumerateArray()
            .Select(r => r[source].GetString()!).Distinct().OrderBy(x => x).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("document_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(24));
            Assert.That(provenance, Is.EqualTo(new[] { "colon", "malaria" }));
        });
    }

    [Test]
    public void Pool_MergedColumnsAreTheUnionOfBothDocuments()
    {
        string[] first = Columns(Invoke("sdrf", "read", "--path", Sdrf("PXD000070.sdrf.tsv")));
        string[] second = Columns(Invoke("sdrf", "read", "--path", Sdrf("PXD026824.sdrf.tsv")));
        string[] merged = Columns(InvokeWithStdin(
            $"{Sdrf("PXD000070.sdrf.tsv")}\n{Sdrf("PXD026824.sdrf.tsv")}", "sdrf", "pool"));

        Assert.Multiple(() =>
        {
            Assert.That(merged, Is.SupersetOf(first.Distinct()));
            Assert.That(merged, Is.SupersetOf(second.Distinct()));
            Assert.That(merged, Contains.Item("comment[source document]"));
        });
    }

    [Test]
    public void Pool_WithoutLabels_WarnsThatProvenanceDependsOnWhereTheFilesSit()
    {
        // mzLib's default label is containing-folder/file-stem, which is not reproducible on
        // another machine. Silence here would let someone publish a table keyed on a local path.
        JsonElement data = InvokeWithStdin(
            $"{Sdrf("PXD000070.sdrf.tsv")}\n{Sdrf("PXD026824.sdrf.tsv")}", "sdrf", "pool");

        Assert.That(Caveats(data), Has.Some.Contains("not reproducible on another machine"));
    }

    [Test]
    public void Pool_WithLabels_DropsThatWarning()
    {
        JsonElement data = InvokeWithStdin(
            $"{Sdrf("PXD000070.sdrf.tsv")}\tmalaria\n{Sdrf("PXD026824.sdrf.tsv")}\tcolon",
            "sdrf", "pool");

        Assert.That(Caveats(data), Has.None.Contains("not reproducible on another machine"));
    }

    [Test]
    public void Pool_PartiallyLabelled_IsRefusedRatherThanMixingTheTwoSchemes()
    {
        // Half the rows stamped with a chosen name and half with a path-derived one is a table
        // whose provenance column means two different things, which nothing downstream can undo.
        JsonElement error = InvokeWithStdinExpectingError(
            $"{Sdrf("PXD000070.sdrf.tsv")}\tmalaria\n{Sdrf("PXD026824.sdrf.tsv")}",
            "sdrf", "pool");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("1 of 2"));
        });
    }

    [Test]
    public void Pool_WritesTheWholeDocument_NotTheWindowedSlice()
    {
        // --limit shapes the payload, never the file. A written document that silently held only
        // the visible rows would disagree with the one the caller asked to pool.
        string output = Path.Combine(Path.GetTempPath(), $"pymzlib-pool-{Guid.NewGuid():N}.sdrf.tsv");
        try
        {
            JsonElement data = InvokeWithStdin(
                $"{Sdrf("PXD000070.sdrf.tsv")}\tmalaria\n{Sdrf("PXD026824.sdrf.tsv")}\tcolon",
                "sdrf", "pool", "--limit", "2", "--out", output);

            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("written").GetProperty("row_count").GetInt32(), Is.EqualTo(24));

            // 24 rows plus the header line, read back from disk rather than trusted from the report.
            Assert.That(File.ReadAllLines(output), Has.Length.EqualTo(25));
        }
        finally
        {
            File.Delete(output);
        }
    }

    [Test]
    public void Pool_WhatItWritesCanBeReadBackAsSdrf()
    {
        // The round trip is the point of writing at all: a pooled table nothing can reopen is a
        // dead end. Reading it back also proves the merged header survived the write.
        string output = Path.Combine(Path.GetTempPath(), $"pymzlib-pool-{Guid.NewGuid():N}.sdrf.tsv");
        try
        {
            JsonElement pooled = InvokeWithStdin(
                $"{Sdrf("PXD000070.sdrf.tsv")}\tmalaria\n{Sdrf("PXD026824.sdrf.tsv")}\tcolon",
                "sdrf", "pool", "--out", output);

            JsonElement reread = Invoke("sdrf", "read", "--path", output);

            Assert.Multiple(() =>
            {
                Assert.That(Columns(reread), Is.EqualTo(Columns(pooled)));
                Assert.That(reread.GetProperty("row_count").GetInt32(),
                    Is.EqualTo(pooled.GetProperty("row_count").GetInt32()));
            });
        }
        finally
        {
            File.Delete(output);
        }
    }

    [Test]
    public void Pool_WithNoStdin_IsAUsageError()
    {
        JsonElement error = InvokeWithStdinExpectingError("", "sdrf", "pool");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("on stdin"));
        });
    }

    [Test]
    public void Pool_NamingAMissingFile_SaysWhichOne()
    {
        string missing = Path.Combine(Path.GetTempPath(), "no-such-file.sdrf.tsv");

        JsonElement error = InvokeWithStdinExpectingError(
            $"{Sdrf("PXD000070.sdrf.tsv")}\n{missing}", "sdrf", "pool");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("no-such-file"));
        });
    }

    // ---- helpers ------------------------------------------------------------------------------

    /// <summary>An SDRF fixture from the pinned mzLib worktree, or an ignored test without one.</summary>
    private static string Sdrf(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        string path = Path.Combine(root, "code", "mzLib", "mzLib", "Test", "FileReadingTests",
            "ExternalFileTypes", name);
        if (!File.Exists(path))
            Assert.Ignore($"mzLib SDRF fixture not present in the worktree: {path}");
        return path;
    }

    private static string[] Columns(JsonElement data) =>
        data.GetProperty("column_names").EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static string[] Caveats(JsonElement data) =>
        data.GetProperty("caveats").EnumerateArray().Select(x => x.GetString()!).ToArray();

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

    /// <summary>
    /// Dispatches a verb with <paramref name="stdin"/> standing in for the console, and shapes the
    /// result the way <c>Program.Main</c> would.
    /// </summary>
    private static JsonElement Envelope(string? stdin, string[] args)
    {
        TextReader previousIn = Console.In;
        if (stdin is not null)
            Console.SetIn(new StringReader(stdin));
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
