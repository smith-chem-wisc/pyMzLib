using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>proteins parse-accession</c>, the projection of mzLib's
/// <c>VariantApplication.ParseAccession</c> (#1382) and <c>ProteinAccession.Parse</c>.
/// </summary>
/// <remarks>
/// The probe covers each case the spec's caveats name: UniProt entries, isoforms and 10-character
/// accessions, RefSeq versions with their underscore kept, variant suffixes (several, a deletion, a
/// stop gain), and every input mzLib leaves unrecognized rather than repairs. Each line must reach
/// mzLib exactly as read, so a leading space is still there in <c>input</c>.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ProteinsAccessionTests
{
    // (input, entry_accession, isoform, version, namespace, parent_accession, applied_variants)
    private static readonly (string Input, string Entry, int? Isoform, int? Version, string Namespace, string Parent, string? Variants)[] Probe =
    {
        ("P04637", "P04637", null, null, "uniprot", "P04637", null),
        ("P04637-2", "P04637", 2, null, "uniprot", "P04637-2", null),
        ("A0A023GPI8", "A0A023GPI8", null, null, "uniprot", "A0A023GPI8", null),
        ("P04637_S70N_A80T", "P04637", null, null, "uniprot", "P04637", "S70N_A80T"),
        ("P04637-2_Q5*", "P04637", 2, null, "uniprot", "P04637-2", "Q5*"),
        ("P04637_A80", "P04637", null, null, "uniprot", "P04637", "A80"),
        ("NP_000537", "NP_000537", null, null, "refseq", "NP_000537", null),
        ("NP_000537.3_R72P", "NP_000537", null, 3, "refseq", "NP_000537.3", "R72P"),
        ("P04637_2", "P04637_2", null, null, "unrecognized", "P04637_2", null),
        ("DECOY_P04637", "DECOY_P04637", null, null, "unrecognized", "DECOY_P04637", null),
        ("sp|P04637|P53_HUMAN", "sp|P04637|P53_HUMAN", null, null, "unrecognized", "sp|P04637|P53_HUMAN", null),
        ("p04637", "p04637", null, null, "unrecognized", "p04637", null),
        (" P04637", " P04637", null, null, "unrecognized", " P04637", null),
        ("P04637|P02768", "P04637|P02768", null, null, "unrecognized", "P04637|P02768", null),
        ("NM_000546.6", "NM_000546.6", null, null, "unrecognized", "NM_000546.6", null),
        ("ENSP00000269305", "ENSP00000269305", null, null, "unrecognized", "ENSP00000269305", null),
        ("P04637-99999999999", "P04637", null, null, "uniprot", "P04637-99999999999", null),
    };

    private static string ProbeStdin => string.Concat(Probe.Select(p => p.Input + "\n"));

    [Test]
    public void EveryProbeAccession_ReadsAsMzLibParsesIt()
    {
        JsonElement data = Invoke(ProbeStdin, "proteins", "parse-accession");
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(new[]
            {
                "input", "entry_accession", "isoform", "version", "namespace", "parent_accession",
                "applied_variants", "has_applied_variants",
            }));
            Assert.That(Strings(columns.GetProperty("input")), Is.EqualTo(Probe.Select(p => p.Input)));
            Assert.That(Strings(columns.GetProperty("entry_accession")), Is.EqualTo(Probe.Select(p => p.Entry)));
            Assert.That(Ints(columns.GetProperty("isoform")), Is.EqualTo(Probe.Select(p => p.Isoform)));
            Assert.That(Ints(columns.GetProperty("version")), Is.EqualTo(Probe.Select(p => p.Version)));
            Assert.That(Strings(columns.GetProperty("namespace")), Is.EqualTo(Probe.Select(p => p.Namespace)));
            Assert.That(Strings(columns.GetProperty("parent_accession")), Is.EqualTo(Probe.Select(p => p.Parent)));
            Assert.That(Strings(columns.GetProperty("applied_variants")), Is.EqualTo(Probe.Select(p => p.Variants)));
            Assert.That(columns.GetProperty("has_applied_variants").EnumerateArray().Select(e => e.GetBoolean()),
                Is.EqualTo(Probe.Select(p => p.Variants is not null)));
        });
    }

    [Test]
    public void Envelope_CountsEveryNamespaceAndTheVariantAccessions()
    {
        JsonElement data = Invoke(ProbeStdin, "proteins", "parse-accession");
        JsonElement counts = data.GetProperty("namespace_counts");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(Probe.Length));
            Assert.That(counts.GetProperty("uniprot").GetInt32(), Is.EqualTo(7));
            Assert.That(counts.GetProperty("refseq").GetInt32(), Is.EqualTo(2));
            Assert.That(counts.GetProperty("unrecognized").GetInt32(), Is.EqualTo(8));
            Assert.That(data.GetProperty("with_applied_variants_count").GetInt32(), Is.EqualTo(4));
            Assert.That(data.GetProperty("caveats").GetArrayLength(), Is.EqualTo(7));
        });
    }

    [Test]
    public void EveryNamespaceIsCounted_EvenAtZero()
    {
        JsonElement counts = Invoke("P04637\n", "proteins", "parse-accession").GetProperty("namespace_counts");

        Assert.That(counts.EnumerateObject().Select(p => (p.Name, p.Value.GetInt32())),
            Is.EquivalentTo(new[] { ("unrecognized", 0), ("uniprot", 1), ("refseq", 0) }));
    }

    [Test]
    public void Duplicates_AreKept_AndBlankLinesSkipped_InInputOrder()
    {
        JsonElement data = Invoke("P04637\n\n   \nNP_000537.3\nP04637\n", "proteins", "parse-accession");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(3));
            Assert.That(Strings(data.GetProperty("columns").GetProperty("input")),
                Is.EqualTo(new[] { "P04637", "NP_000537.3", "P04637" }));
        });
    }

    [Test]
    public void ATrailingSpace_ReachesMzLibAndIsNotRepaired()
    {
        JsonElement columns = Invoke("P04637 \n", "proteins", "parse-accession").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("input")[0].GetString(), Is.EqualTo("P04637 "));
            Assert.That(columns.GetProperty("namespace")[0].GetString(), Is.EqualTo("unrecognized"));
        });
    }

    [TestCase("")]
    [TestCase("\n \n")]
    public void NoAccessions_IsAUsageError(string stdin)
    {
        JsonElement error = InvokeExpectingError(stdin, "proteins", "parse-accession");

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("No accessions on stdin"));
    }

    // ---- plumbing -------------------------------------------------------------------------------

    private static string?[] Strings(JsonElement array) =>
        array.EnumerateArray().Select(e => e.ValueKind == JsonValueKind.Null ? null : e.GetString()).ToArray();

    private static int?[] Ints(JsonElement array) =>
        array.EnumerateArray().Select(e => e.ValueKind == JsonValueKind.Null ? (int?)null : e.GetInt32()).ToArray();

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
