using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>peptidoform convert</c>, the projection of mzLib's <c>SequenceConversionService</c>.
/// </summary>
/// <remarks>
/// The accessions are mzLib's. These tests pin the cases the verb exists for (UniProt-sourced
/// modifications to Unimod), the ProForma gap the caveats describe (so an mzLib fix is noticed
/// here), what each <c>SequenceConversionHandlingMode</c> does to a row, and that rows come back in
/// input order at any thread count.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PeptidoformConvertTests
{
    private const string Probe =
        "[UniProt:N-acetylserine on S]SEQK\n" +
        "PEPK[UniProt:N6,N6-dimethyllysine on K]R\n" +
        "PEPM[Common Variable:Oxidation on M]K\n" +
        "PEPK[Made Up:Not a modification on K]R\n";

    [Test]
    public void Unimod_ResolvesUniProtAndMetaMorpheusModifications()
    {
        JsonElement data = Invoke(Probe, "peptidoform", "convert");
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("source_format").GetString(), Is.EqualTo("mzLib"));
            Assert.That(data.GetProperty("target_format").GetString(), Is.EqualTo("Unimod"));
            Assert.That(data.GetProperty("mode").GetString(), Is.EqualTo("ReturnNull"));
            Assert.That(Strings(columns.GetProperty("output")), Is.EqualTo(new[]
            {
                "[UNIMOD:1]SEQK",     // N-acetylserine
                "PEPK[UNIMOD:36]R",   // N6,N6-dimethyllysine: Dimethyl, not Ethyl (UNIMOD:280)
                "PEPM[UNIMOD:35]K",   // Oxidation on M
                null,
            }));
            Assert.That(Strings(columns.GetProperty("status")), Is.EqualTo(new[]
                { "converted", "converted", "converted", "failed" }));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(4));
            Assert.That(data.GetProperty("converted_count").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("warned_count").GetInt32(), Is.EqualTo(0));
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
        });
    }

    [Test]
    public void Unimod_UnresolvableModificationFailsItsRowAndIsNamed()
    {
        JsonElement columns = Invoke(Probe, "peptidoform", "convert").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("input")[3].GetString(), Is.EqualTo("PEPK[Made Up:Not a modification on K]R"));
            Assert.That(columns.GetProperty("output")[3].ValueKind, Is.EqualTo(JsonValueKind.Null));
            // mzLib's Unimod serializer under ReturnNull records the item but no reason code.
            Assert.That(columns.GetProperty("failure_reason")[3].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(Strings(columns.GetProperty("incompatible_items")[3]),
                Is.EqualTo(new[] { "Made Up:Not a modification on K @3(K)" }));
            Assert.That(columns.GetProperty("incompatible_items")[0].GetArrayLength(), Is.Zero);
        });
    }

    [Test]
    public void Envelope_ListsTheRegisteredFormatsAndColumns()
    {
        JsonElement data = Invoke("PEPTIDE\n", "peptidoform", "convert");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data.GetProperty("source_formats")), Is.EqualTo(new[] { "MassShift", "Modomics", "ProForma", "mzLib" }));
            Assert.That(Strings(data.GetProperty("target_formats")), Is.EqualTo(new[]
                { "Chronologer", "Essential", "MassShift", "ProForma", "Unimod", "mzLib" }));
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(new[]
                { "input", "output", "status", "failure_reason", "incompatible_items", "warnings", "errors" }));
            Assert.That(data.GetProperty("caveats").GetArrayLength(), Is.EqualTo(6));
            Assert.That(data.GetProperty("columns").GetProperty("output")[0].GetString(), Is.EqualTo("PEPTIDE"));
        });
    }

    [Test]
    public void ProForma_LeavesUniProtModificationsUnresolved()
    {
        // The mzLib gap the first caveat describes. When mzLib fixes it this test fails, and the
        // caveat, the spec and the guide come out together.
        JsonElement data = Invoke(Probe, "peptidoform", "convert", "--to", "ProForma");
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("target_format").GetString(), Is.EqualTo("ProForma"));
            Assert.That(columns.GetProperty("output")[0].GetString(), Is.EqualTo("[UniProt:N-acetylserine on S]-SEQK"));
            Assert.That(columns.GetProperty("output")[1].GetString(), Is.EqualTo("PEPK[UniProt:N6,N6-dimethyllysine on K]R"));
            Assert.That(columns.GetProperty("output")[2].GetString(), Is.EqualTo("PEPM[UNIMOD:35]K"));
            Assert.That(columns.GetProperty("output")[3].GetString(), Is.EqualTo("PEPK[Made Up:Not a modification on K]R"));
            Assert.That(Strings(columns.GetProperty("status")), Is.All.EqualTo("converted"));
        });
    }

    [TestCase("RemoveIncompatibleElements", true)]
    [TestCase("UsePrimarySequence", false)]
    public void DroppingModes_KeepTheRowAndSayWhatWasDropped(string mode, bool mzLibWarns)
    {
        JsonElement data = Invoke(Probe, "peptidoform", "convert", "--mode", mode);
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("mode").GetString(), Is.EqualTo(mode));
            Assert.That(columns.GetProperty("output")[3].GetString(), Is.EqualTo("PEPKR"));
            Assert.That(columns.GetProperty("status")[3].GetString(), Is.EqualTo("converted_with_warnings"));
            Assert.That(Strings(columns.GetProperty("incompatible_items")[3]),
                Is.EqualTo(new[] { "Made Up:Not a modification on K @3(K)" }));
            Assert.That(columns.GetProperty("warnings")[3].GetArrayLength() > 0, Is.EqualTo(mzLibWarns));
            Assert.That(data.GetProperty("warned_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.Zero);
        });
    }

    [TestCase(1)]
    [TestCase(4)]
    public void ThrowException_FailsTheCallOnTheFirstFailingInputInOrder(int threads)
    {
        string stdin = "PEPTIDE\n" + "PEPK[Made Up:First on K]R\n" + string.Concat(Enumerable.Repeat("PEPTIDE\n", 50)) +
            "PEPK[Made Up:Second on K]R\n";
        JsonElement error = InvokeExpectingError(stdin, "peptidoform", "convert", "--mode", "throwexception",
            "--threads", threads.ToString());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("'PEPK[Made Up:First on K]R'")
            .And.Contain("IncompatibleModifications").And.Contain("ReturnNull"));
    }

    [Test]
    public void Unparseable_FailsItsRowWithMzLibsReason()
    {
        JsonElement columns = Invoke("PEP[unclosed\n", "peptidoform", "convert").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("status")[0].GetString(), Is.EqualTo("failed"));
            Assert.That(columns.GetProperty("failure_reason")[0].GetString(), Is.EqualTo("UnknownFormat"));
            Assert.That(Strings(columns.GetProperty("errors")[0]), Has.Some.Contains("Unclosed bracket"));
        });
    }

    [Test]
    public void AmbiguousFullSequence_IsJoinedByMzLibWithAWarning()
    {
        // Documented in the caveats: mzLib's parser skips '|' rather than refusing the input.
        JsonElement columns = Invoke("PEPTIDE|PEPTLDE\n", "peptidoform", "convert").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("output")[0].GetString(), Is.EqualTo("PEPTIDEPEPTLDE"));
            Assert.That(columns.GetProperty("status")[0].GetString(), Is.EqualTo("converted_with_warnings"));
            Assert.That(Strings(columns.GetProperty("warnings")[0]), Has.Some.Contains("'|'"));
        });
    }

    [Test]
    public void FormatNames_MatchCaseInsensitivelyAndEchoMzLibsSpelling()
    {
        JsonElement data = Invoke("PEPM[Common Variable:Oxidation on M]K\n", "peptidoform", "convert",
            "--from", "MZLIB", "--to", "unimod");

        Assert.That(data.GetProperty("source_format").GetString(), Is.EqualTo("mzLib"));
        Assert.That(data.GetProperty("target_format").GetString(), Is.EqualTo("Unimod"));
        Assert.That(data.GetProperty("columns").GetProperty("output")[0].GetString(), Is.EqualTo("PEPM[UNIMOD:35]K"));
    }

    [Test]
    public void Threads_DoNotChangeTheRows()
    {
        string stdin = string.Concat(Enumerable.Range(0, 40).Select(_ => Probe));
        string one = Invoke(stdin, "peptidoform", "convert", "--threads", "1").GetRawText();
        string four = Invoke(stdin, "peptidoform", "convert", "--threads", "4").GetRawText();
        string all = Invoke(stdin, "peptidoform", "convert", "--threads", "-1").GetRawText();

        Assert.That(four, Is.EqualTo(one));
        Assert.That(all, Is.EqualTo(one));
    }

    [TestCase("PEPTIDE\n", new[] { "--to", "Foo" }, "no target format named 'Foo'")]
    [TestCase("PEPTIDE\n", new[] { "--from", "Unimod" }, "no source format named 'Unimod'")]
    [TestCase("PEPTIDE\n", new[] { "--to" }, "--to was given without a value")]
    [TestCase("PEPTIDE\n", new[] { "--mode", "Strict" }, "SequenceConversionHandlingMode")]
    [TestCase("PEPTIDE\n", new[] { "--mode", "1" }, "SequenceConversionHandlingMode")]
    [TestCase("PEPTIDE\n", new[] { "--mode" }, "--mode was given without a value")]
    [TestCase("PEPTIDE\n", new[] { "--threads", "0" }, "--threads")]
    [TestCase("\n  \n", new string[0], "No sequences on stdin")]
    public void BadCalls_AreUsageErrors(string stdin, string[] options, string expected)
    {
        JsonElement error = InvokeExpectingError(stdin, new[] { "peptidoform", "convert" }.Concat(options).ToArray());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
    }

    // ---- plumbing -------------------------------------------------------------------------------

    private static string?[] Strings(JsonElement array) => array.EnumerateArray().Select(e => e.GetString()).ToArray();

    private static JsonElement Invoke(string stdin, params string[] args) => Unwrap(Envelope(stdin, args), expectOk: true);

    private static JsonElement InvokeExpectingError(string stdin, params string[] args) => Unwrap(Envelope(stdin, args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk), $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    /// <summary>Dispatches a verb with <paramref name="stdin"/> as the console, shaped as <c>Program.Main</c> would.</summary>
    private static JsonElement Envelope(string stdin, string[] args)
    {
        TextReader previousIn = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.SerializeToElement(new { ok = true, data }, Program.JsonOptions);
        }
        catch (Program.UsageException usage)
        {
            return JsonSerializer.SerializeToElement(new { ok = false, error = new { type = "usage", message = usage.Message } },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetIn(previousIn);
        }
    }
}
