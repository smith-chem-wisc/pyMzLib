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
    private Func<string, Task<(string, bool)>> _originalSource = null!;

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

    /// <summary>
    /// An entry carrying a signal-peptide proteolysis product (residues 1-18). The only tryptic
    /// site is K26, so the peptide ending at the signal boundary (1-18) exists ONLY because mzLib
    /// digests at the proteolysis-product boundary — a protease-only digest of the bare sequence
    /// would run 1-26 straight through. That makes it the perfect probe for pyMzLib#8.
    /// </summary>
    private const string SignalPeptideEntryXml = """
        <?xml version="1.0" encoding="UTF-8"?>
        <uniprot xmlns="http://uniprot.org/uniprot">
          <entry dataset="Swiss-Prot">
            <accession>P00002</accession>
            <name>SIGNAL_HUMAN</name>
            <protein><recommendedName><fullName>Signal test protein</fullName></recommendedName></protein>
            <organism><name type="scientific">Homo sapiens</name></organism>
            <feature type="signal peptide"><location><begin position="1"/><end position="18"/></location></feature>
            <feature type="modified residue" description="Phosphoserine"><location><position position="3"/></location></feature>
            <sequence length="30" mass="3000">MASAENLTLQIQTGDISALGDIVATKTQVW</sequence>
          </entry>
        </uniprot>
        """;

    /// <summary>
    /// Residues 1-18 of <see cref="SignalPeptideEntryXml"/> — the exact span of its signal-peptide
    /// feature, so this is DERIVED from the fixture, not a golden value copied from output. If a
    /// digest ever stops producing it, the correct response is to check the fix or the fixture's
    /// single-tryptic-site premise, never to re-copy whatever the code now emits.
    /// </summary>
    private const string BoundaryPeptide = "MASAENLTLQIQTGDISA";

    private void UseXml(string xml)
    {
        string path = Path.Combine(_tempDirectory, "entry.xml");
        File.WriteAllText(path, xml, Encoding.UTF8);
        // CallerDeletes: false - the test owns this fixture, the workflow must not delete it.
        Peptidoform.UniProtXmlSource = _ => Task.FromResult((path, false));
    }

    /// <summary>
    /// The distinct digestion triples — (start residue, end residue, base sequence) — a fragments
    /// result reports. Position, not just letters: truncation products control peptide BOUNDARIES,
    /// so a dropped boundary peptide whose sequence coincided with another's would slip past a
    /// sequence-only comparison but not this one.
    /// </summary>
    private static HashSet<(int start, int end, string bare)> DigestionTriples(JsonElement result) =>
        result.GetProperty("peptides").EnumerateArray()
            .Select(p => (p.GetProperty("one_based_start").GetInt32(),
                          p.GetProperty("one_based_end").GetInt32(),
                          p.GetProperty("base_sequence").GetString()!))
            .ToHashSet();

    /// <summary>Each reported peptide's (base sequence, full sequence) pair.</summary>
    private static IEnumerable<(string bare, string full)> BareAndFullSequences(JsonElement result) =>
        result.GetProperty("peptides").EnumerateArray()
            .Select(p => (p.GetProperty("base_sequence").GetString()!,
                          p.GetProperty("full_sequence").GetString()!));

    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        object data = await Program.DispatchAsync(args);
        return JsonSerializer.Deserialize<JsonElement>(JsonSerializer.Serialize(data, Program.JsonOptions));
    }

    // ---- the census ----------------------------------------------------------

    [Test]
    public async Task AnInjectedFixtureIsNotDeletedByTheWorkflow()
    {
        // The workflow used to delete whatever path the source returned, destroying a fixture a
        // test reused across two calls. Ownership is now explicit; prove the file survives.
        UseXml(MiniEntryXml);
        string fixture = Path.Combine(_tempDirectory, "entry.xml");

        await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");

        Assert.That(File.Exists(fixture), Is.True, "the workflow deleted a fixture it did not own");
    }

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
        // once made a histone look as though 93 annotations had been silently dropped. Checking
        // only that both keys exist would pass even if they were the same field twice; the mini
        // entry has two modified residues at two positions, so both counts must equal 2 and be
        // reported independently.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync("peptidoform", "fragments", "--accession", "P00001");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("annotated_modification_sites").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("annotated_modifications_loaded").GetInt32(), Is.EqualTo(2));
        });
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
    public async Task NoModificationsPreservesProteolysisProductsSoTheControlDiffersOnlyByMods()
    {
        // pyMzLib#8: the --no-modifications control rebuilt the protein through a Protein constructor
        // that dropped ProteolysisProducts, so mzLib stopped digesting at the signal-peptide boundary
        // and the peptide LIST changed, not just its modifications. A control is only a control if a
        // single variable moves. With the products carried across, the two runs must digest to the
        // same BOUNDARIES and differ only in modifications.
        //
        // The fixture's single-tryptic-site property (K26 is the ONLY K/R) is load-bearing: it is why
        // the 1-18 boundary peptide can arise solely from the truncation-product boundary. A K/R added
        // inside residues 1-18 would let ordinary cleavage produce it and silently defang this test —
        // which is why the boundary peptide is asserted explicitly below, so such a fixture edit
        // breaks a visible assertion rather than quietly weakening coverage. --max-mods is pinned so
        // the fixture's Phosphoserine is placeable regardless of future default changes.
        UseXml(SignalPeptideEntryXml);
        JsonElement with = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00002", "--min-length", "7", "--max-mods", "2");

        UseXml(SignalPeptideEntryXml);
        JsonElement without = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00002", "--min-length", "7", "--max-mods", "2",
            "--no-modifications");

        HashSet<(int start, int end, string bare)> withTriples = DigestionTriples(with);
        HashSet<(int start, int end, string bare)> withoutTriples = DigestionTriples(without);

        Assert.Multiple(() =>
        {
            // Non-vacuous, and pinned by position: the 1-18 span is digested only at the truncation
            // boundary (verified — the pre-fix control's set is missing precisely this triple).
            Assert.That(withTriples, Is.Not.Empty, "the fixture must yield peptides for this test to mean anything");
            Assert.That(withTriples, Does.Contain((1, 18, BoundaryPeptide)),
                "the annotated run must contain the 1-18 signal-peptide boundary peptide");

            // The fix: the control keeps those boundary peptides — same (start, end, base sequence)
            // set — instead of silently dropping them.
            Assert.That(withoutTriples, Does.Contain((1, 18, BoundaryPeptide)),
                "the --no-modifications control must keep the 1-18 boundary peptide, not drop it");
            Assert.That(withoutTriples, Is.EquivalentTo(withTriples),
                "annotated and control must digest to the same (start, end, base sequence) boundaries");

            // The OTHER half of the invariant: the control must actually be STRIPPED of modifications,
            // not merely equal by base sequence. Without this, a regression that made --no-modifications
            // return the annotated protein unchanged would collapse both runs together and pass every
            // assertion above. The annotated run must place at least one mod (Phosphoserine at S3, which
            // sits inside the 1-18 peptide); the control must carry none.
            Assert.That(BareAndFullSequences(with).Any(p => p.full != p.bare), Is.True,
                "the annotated run must apply at least one modification (full_sequence != base_sequence)");
            Assert.That(BareAndFullSequences(without).All(p => p.full == p.bare), Is.True,
                "the --no-modifications control must carry no modifications (full_sequence == base_sequence)");
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
        // A truncated peptidoform list and a short one look identical from outside. Asserting only
        // that the field exists would pass even if truncation were never detected, so this forces
        // the cap to bind - max-isoforms 1 on a peptide with two possible modifications - and
        // requires a non-zero count, a peptide actually reported at the cap.
        UseXml(MiniEntryXml);

        JsonElement data = await InvokeAsync(
            "peptidoform", "fragments", "--accession", "P00001", "--min-length", "4",
            "--max-mods", "2", "--max-isoforms", "1");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("max_modification_isoforms").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("peptides_at_isoform_cap").GetInt32(), Is.GreaterThan(0),
                "the cap was forced to bind but no peptide was reported at it");
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
