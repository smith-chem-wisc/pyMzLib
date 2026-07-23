using System.Diagnostics.CodeAnalysis;
using System.Text;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the Peptidoform workflow, driven from a local UniProt XML rather than the network.
/// </summary>
/// <remarks>
/// Everything worth checking here — the annotation census, digestion, modification combinatorics,
/// the isoform cap — sits downstream of the download, so none of it should need EBI or UniProt to
/// be reachable. <see cref="Peptidoform.UniProtXmlSource"/> exists for that, the same way
/// <see cref="Program.PrideClientFactory"/> does.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PeptidoformTests
{
    private string _tempDirectory = string.Empty;
    private Func<string, Task<string>> _originalSource = null!;

    [SetUp]
    public void CreateTempDirectory()
    {
        _originalSource = Peptidoform.UniProtXmlSource;
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void Cleanup()
    {
        Peptidoform.UniProtXmlSource = _originalSource;
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    /// <summary>A minimal but structurally real UniProt entry: two modified residues, one glycosylation site.</summary>
    private const string MiniEntryXml = """
        <?xml version="1.0" encoding="UTF-8"?>
        <uniprot xmlns="http://uniprot.org/uniprot">
          <entry dataset="Swiss-Prot">
            <accession>P00001</accession>
            <name>TEST_HUMAN</name>
            <protein><recommendedName><fullName>Test protein</fullName></recommendedName></protein>
            <organism><name type="scientific">Homo sapiens</name></organism>
            <feature type="modified residue" description="Phosphoserine"><location><position position="3"/></location></feature>
            <feature type="modified residue" description="Phosphothreonine"><location><position position="12"/></location></feature>
            <feature type="glycosylation site" description="N-linked (GlcNAc...) asparagine"><location><position position="7"/></location></feature>
            <sequence length="24" mass="2600">MASRENKTLIQTGDKVWERSAMKR</sequence>
          </entry>
        </uniprot>
        """;

    private void UseXml(string xml)
    {
        string path = Path.Combine(_tempDirectory, "entry.xml");
        File.WriteAllText(path, xml, Encoding.UTF8);
        Peptidoform.UniProtXmlSource = _ => Task.FromResult(path);
    }

    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        object data = await Program.DispatchAsync(args);
        return JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(data, Program.JsonOptions));
    }

    // ---- the census ----------------------------------------------------------

    [Test]
    public async Task TheCensusReportsAnnotatedFeaturesIncludingTheOnesNotUsed()
    {
        // The whole reason the census exists: a correct count produced by an invisible rule is
        // still a trap. Two modified residues are usable; the glycosylation site has no defined
        // mass and is not — and the caller should be able to learn that rather than infer it.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("uniprot_annotated_features").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("annotated_modifications_loaded").GetInt32(), Is.EqualTo(2));
        });
    }

    [Test]
    public async Task TheCensusNamesEachFeatureTypeAndWhetherItWasLoaded()
    {
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");
        var byType = data.GetProperty("uniprot_features_by_type").EnumerateArray()
            .ToDictionary(e => e.GetProperty("type").GetString()!, e => e.GetProperty("loaded").GetBoolean());

        Assert.Multiple(() =>
        {
            Assert.That(byType["modified residue"], Is.True);
            Assert.That(byType["glycosylation site"], Is.False);
        });
    }

    [Test]
    public async Task SitesAndModificationsAreReportedSeparately()
    {
        // They are different numbers whenever a residue carries alternatives, and conflating them
        // once made a histone look as though 93 annotations had been silently dropped.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");

        Assert.That(data.TryGetProperty("annotated_modification_sites", out _), Is.True);
        Assert.That(data.TryGetProperty("annotated_modifications_loaded", out _), Is.True);
    }

    // ---- the workflow --------------------------------------------------------

    [Test]
    public async Task DigestionAndFragmentationProduceFragmentsForEveryPeptide()
    {
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", "--min-length", "4");

        Assert.That(data.GetProperty("peptide_count").GetInt32(), Is.GreaterThan(0));
        foreach (JsonElement peptide in data.GetProperty("peptides").EnumerateArray())
        {
            Assert.That(peptide.GetProperty("fragments").GetArrayLength(), Is.GreaterThan(0),
                $"no fragments for {peptide.GetProperty("base_sequence").GetString()}");
        }
    }

    [Test]
    public async Task ModificationsChangeTheResult()
    {
        UseXml(MiniEntryXml);
        JsonElement with = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001", "--min-length", "4");

        UseXml(MiniEntryXml);
        JsonElement without = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", "--min-length", "4", "--no-modifications");

        Assert.Multiple(() =>
        {
            Assert.That(with.GetProperty("modifications_applied").GetBoolean(), Is.True);
            Assert.That(without.GetProperty("modifications_applied").GetBoolean(), Is.False);
            Assert.That(with.GetProperty("peptide_count").GetInt32(),
                Is.GreaterThanOrEqualTo(without.GetProperty("peptide_count").GetInt32()));
        });
    }

    [Test]
    public async Task TheDefaultProteaseAppliesTheProlineRule()
    {
        // mzLib's "trypsin|P" is the Keil rule and "trypsin" is not — the reverse of the
        // MaxQuant/Mascot convention. The default here must be the one a mass spectrometrist means.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");

        Assert.That(data.GetProperty("protease").GetString(), Is.EqualTo("trypsin|P"));
    }

    [Test]
    public async Task TheIsoformCapIsReportedSoTruncationIsVisible()
    {
        // A truncated Peptidoform list and a short one look identical from outside.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", "--min-length", "4", "--max-isoforms", "1");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("max_modification_isoforms").GetInt32(), Is.EqualTo(1));
            Assert.That(data.TryGetProperty("peptides_at_isoform_cap", out _), Is.True);
        });
    }

    // ---- input validation ----------------------------------------------------

    [Test]
    public void AnUnknownProteaseNamesSomeAlternatives()
    {
        UseXml(MiniEntryXml);

        var ex = Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", "--protease", "banana"));

        Assert.That(ex!.Message, Does.Contain("trypsin"));
    }

    [TestCase("--dissociation", "banana")]
    [TestCase("--terminus", "sideways")]
    public void UnknownEnumValuesAreUsageErrorsListingWhatIsValid(string option, string value)
    {
        UseXml(MiniEntryXml);

        var ex = Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", option, value));

        Assert.That(ex!.Message, Does.Contain("Known"));
    }

    [Test]
    public void AMissingAccessionIsAUsageError()
    {
        UseXml(MiniEntryXml);
        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync("peptidoform", "fragments"));
    }

    // ---- UniProt availability -------------------------------------------------
    //
    // UniProt is a second external service and gets the same treatment as PRIDE: a permanent
    // caller error must never be dressed as an outage, and a genuine outage must never be dressed
    // as a caller error. Getting this backwards is what made the live suites skip real bugs.

    [TestCase(System.Net.HttpStatusCode.NotFound, Description = "well-formed but unknown accession")]
    [TestCase(System.Net.HttpStatusCode.BadRequest, Description = "malformed accession")]
    public void APermanentAccessionProblemIsAUsageErrorNotAnOutage(System.Net.HttpStatusCode status)
    {
        var ex = Assert.Throws<Program.UsageException>(() =>
            Peptidoform.ThrowIfUniProtRejected(status, "whatever", "P99999999", "https://x/y"));

        Assert.That(ex!.Message, Does.Contain("P99999999"));
    }

    [TestCase(500)]
    [TestCase(502)]
    [TestCase(503)]
    [TestCase(429)]
    [TestCase(408)]
    public void AUniProtOutageClassifiesAsServiceUnavailable(int status)
    {
        // Phrased so ClassifyError can read the code back out of the message — which is the
        // contract that lets both test suites skip rather than fail when UniProt is down.
        var thrown = Assert.Throws<HttpRequestException>(() => Peptidoform.ThrowIfUniProtRejected(
            (System.Net.HttpStatusCode)status, "Service Unavailable", "P02768", "https://x/y"));

        Assert.That(Program.ClassifyError(thrown!), Is.EqualTo(Program.ServiceUnavailableType));
    }

    [TestCase(401)]
    [TestCase(403)]
    [TestCase(418)]
    public void AUniProtClientErrorStaysOurProblem(int status)
    {
        var thrown = Assert.Throws<HttpRequestException>(() => Peptidoform.ThrowIfUniProtRejected(
            (System.Net.HttpStatusCode)status, "Nope", "P02768", "https://x/y"));

        Assert.That(Program.ClassifyError(thrown!), Is.EqualTo(nameof(HttpRequestException)));
    }

    [TestCase(200)]
    [TestCase(204)]
    public void ASuccessfulResponsePassesThrough(int status)
    {
        Assert.DoesNotThrow(() => Peptidoform.ThrowIfUniProtRejected(
            (System.Net.HttpStatusCode)status, "OK", "P02768", "https://x/y"));
    }

    [Test]
    public void AnUnreachableUniProtIsAnOutageNotAFailure()
    {
        // A refused connection or DNS failure never reaches ThrowIfUniProtRejected at all — it
        // surfaces from HttpClient with no status and a wrapped cause, and must still be an outage.
        var refused = new HttpRequestException(
            "Connection refused (rest.uniprot.org:443)",
            new System.Net.Sockets.SocketException(10061));

        Assert.That(Program.ClassifyError(refused), Is.EqualTo(Program.ServiceUnavailableType));
    }

    [Test]
    public void AnEntryWithNoSequenceIsReportedRatherThanReturningNothing()
    {
        UseXml("""
            <?xml version="1.0" encoding="UTF-8"?>
            <uniprot xmlns="http://uniprot.org/uniprot"></uniprot>
            """);

        Assert.ThrowsAsync<Program.UsageException>(async () => await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001"));
    }
}
