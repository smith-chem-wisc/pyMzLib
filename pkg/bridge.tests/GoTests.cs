using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.Json;
using Proteomics;
using Readers;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.GeneOntology;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>proteins annotate-go</c> and <c>proteins update-go</c>, on real data.
/// </summary>
/// <remarks>
/// <para>
/// <b>PXD036557_AllQuantifiedProteinGroups.tsv</b> is mzLib's own fixture: six rows of a real
/// MetaMorpheus 1.1.11 search (mzLib #1347), with one two-member group (H2A.Z, <c>P0C0S5|Q71UI9</c>),
/// one contaminant (bovine albumin) and one decoy.
/// <b>pxd036557_proteins.xml</b> is the five human entries those groups name, cut whole and unedited out
/// of UniProt's reviewed human proteome (UP000005640, downloaded 2026-09-18).
/// <b>go-pxd036557.obo</b> is GO release 2026-07-26 with every stanza kept whose id mzLib's annotator
/// emitted when run against the full 48,340-term go.obo and the full 20,416-entry proteome: 412 terms.
/// The closure is mzLib's, not this repository's, and the trimmed pair reproduces all 563 rows of the
/// full run exactly, apart from the two file hashes (checked when the fixtures were made).
/// <b>organelle_map.tsv</b> is a nine-anchor category map in mzLib's format.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class GoTests
{
    private string _temp = "";

    [SetUp]
    public void SetUp()
    {
        _temp = Path.Combine(Path.GetTempPath(), $"pymzlib-go-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_temp);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_temp))
            Directory.Delete(_temp, recursive: true);
    }

    private static string[] Annotate(params string[] extra) =>
        new[] { "proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", Fixture("pxd036557_proteins.xml"), "--go-obo", Fixture("go-pxd036557.obo") }.Concat(extra).ToArray();

    // ---- the table ----------------------------------------------------------------------------

    [Test]
    public void ColumnNames_AreMzLibsOwnSchemaHeaders()
    {
        JsonElement data = Invoke(Annotate("--limit", "0"));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("column_names").EnumerateArray().Select(n => n.GetString()),
                Is.EqualTo(GoAnnotationTsv.Schema.Select(c => c.Header)));
            Assert.That(Proteins.GoAnnotationColumns, Is.EqualTo(GoAnnotationTsv.Schema.Select(c => c.Header)));
        });
    }

    [Test]
    public void RealSearch_EveryNonDecoyGroupGetsRows_AndNoMemberIsPrivileged()
    {
        JsonElement data = Invoke(Annotate());
        JsonElement c = data.GetProperty("columns");
        string[] groups = Strings(c, "protein_group");
        string[] status = Strings(c, "annotation_status");
        int[] nWith = Ints(c, "n_with");
        int[] nMembers = Ints(c, "n_members");

        int[] h2az = Enumerable.Range(0, groups.Length).Where(i => groups[i] == "P0C0S5|Q71UI9").ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("table_row_count").GetInt32(), Is.EqualTo(6));
            Assert.That(data.GetProperty("decoy_group_count").GetInt32(), Is.EqualTo(1), "DECOY_P62750 is skipped");
            Assert.That(data.GetProperty("group_count").GetInt32(), Is.EqualTo(5));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(563));
            Assert.That(groups.Distinct(), Is.EqualTo(new[] { "P68363", "P05141", "P0C0S5|Q71UI9", "P02769", "P63104" }),
                "groups keep the table's order");

            // The contaminant gets exactly one term-less row that says why.
            int albumin = Array.IndexOf(groups, "P02769");
            Assert.That(groups.Count(g => g == "P02769"), Is.EqualTo(1));
            Assert.That(status[albumin], Is.EqualTo("contaminant"));
            Assert.That(c.GetProperty("go_id")[albumin].ValueKind, Is.EqualTo(JsonValueKind.Null));

            // The union is every term either histone carries; the consensus is the 57 both carry.
            Assert.That(h2az, Has.Length.EqualTo(106));
            Assert.That(h2az.Count(i => nWith[i] == nMembers[i]), Is.EqualTo(57));
            Assert.That(h2az.Select(i => nMembers[i]), Is.All.EqualTo(2));
        });
    }

    [Test]
    public void ARow_SaysWhoCarriesTheTerm_HowDirectly_AndOnWhatEvidence()
    {
        JsonElement c = Invoke(Annotate()).GetProperty("columns");
        int innerMembrane = Row(c, "P05141", "GO:0005743");
        int euchromatin = Row(c, "P0C0S5|Q71UI9", "GO:0000791");
        int nucleosome = Row(c, "P0C0S5|Q71UI9", "GO:0000786");

        Assert.Multiple(() =>
        {
            Assert.That(c.GetProperty("go_name")[innerMembrane].GetString(), Is.EqualTo("mitochondrial inner membrane"));
            Assert.That(c.GetProperty("aspect")[innerMembrane].GetString(), Is.EqualTo("cellular_component"));
            Assert.That(c.GetProperty("propagated")[innerMembrane].GetBoolean(), Is.False, "UniProt annotates it directly");
            Assert.That(c.GetProperty("inherited")[innerMembrane].GetBoolean(), Is.False);
            Assert.That(StringList(c, "evidence", innerMembrane), Is.EqualTo(new[] { "ECO:0000250" }));

            // Only one of the two histones is annotated to euchromatin.
            Assert.That(StringList(c, "accession_used", euchromatin), Is.EqualTo(new[] { "P0C0S5" }));
            Assert.That(c.GetProperty("n_with")[euchromatin].GetInt32(), Is.EqualTo(1));

            // Both are annotated to the nucleosome, each on its own evidence.
            JsonElement byMember = c.GetProperty("evidence_by_member")[nucleosome];
            Assert.That(byMember.EnumerateObject().Select(p => p.Name), Is.EqualTo(new[] { "P0C0S5", "Q71UI9" }));
            Assert.That(byMember.GetProperty("Q71UI9")[0].GetString(), Is.EqualTo("ECO:0000353"));
        });
    }

    [Test]
    public void TheHeaderIsMzLibsOwn_AndCountsGroupsNotRows()
    {
        JsonElement data = Invoke(Annotate("--limit", "0"));
        JsonElement header = data.GetProperty("header");

        Assert.Multiple(() =>
        {
            Assert.That(header.GetProperty("go_annotation_format").GetString(), Is.EqualTo("1"));
            Assert.That(header.GetProperty("go_release").GetString(), Is.EqualTo("releases/2026-07-26"));
            Assert.That(header.GetProperty("counter_q_value_max").GetString(), Is.EqualTo("0.01"));
            Assert.That(header.GetProperty("status_annotated").GetString(), Is.EqualTo("4"));
            Assert.That(header.GetProperty("status_contaminant").GetString(), Is.EqualTo("1"));
            Assert.That(header.GetProperty("n_multi_member_groups").GetString(), Is.EqualTo("1"));
            Assert.That(header.GetProperty("source_file_sha256").GetString(),
                Is.EqualTo(data.GetProperty("groups_file_sha256").GetString()));
            Assert.That(header.GetProperty("annotation_db_sha256").GetString(),
                Is.EqualTo(data.GetProperty("annotation_database").GetProperty("sha256").GetString()));
            Assert.That(data.GetProperty("go").GetProperty("term_count").GetInt32(), Is.EqualTo(412));
        });
    }

    [Test]
    public void Out_IsByteForByteTheFileMzLibsWriterProduces()
    {
        string outPath = Path.Combine(_temp, "go.tsv");
        JsonElement data = Invoke(Annotate("--out", outPath, "--limit", "0"));

        // The same run, straight through mzLib with no bridge in between.
        var ontology = GeneOntologyGraph.Load(Fixture("go-pxd036557.obo"));
        List<Protein> proteins = ProteinDbLoader.LoadProteinXML(Fixture("pxd036557_proteins.xml"), true,
            DecoyType.None, Array.Empty<Omics.Modifications.Modification>(), false, null, out _, maxThreads: 1,
            maxHeterozygousVariants: 0);
        string dbSha = Proteins.DecompressedSha256(Fixture("pxd036557_proteins.xml"));
        var annotator = new GoGroupAnnotator(ontology, proteins, dbSha);
        var rows = annotator.AnnotateAll(new ProteinGroupFromTsvFile(Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"))
            .ToGoAnnotationGroups());
        var expected = new StringWriter();
        GoAnnotationTsv.Write(expected, rows, ontology, dbSha, data.GetProperty("groups_file_sha256").GetString());

        Assert.Multiple(() =>
        {
            Assert.That(File.ReadAllText(outPath), Is.EqualTo(expected.ToString()));
            Assert.That(data.GetProperty("written").GetProperty("row_count").GetInt32(), Is.EqualTo(563));
            Assert.That(data.GetProperty("written").GetProperty("path").GetString(), Is.EqualTo(outPath));
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(0));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
            Assert.That(File.Exists(outPath + ".partial"), Is.False);
        });
    }

    [Test]
    public void LimitAndOffset_WindowTheWire_NotTheFile()
    {
        JsonElement all = Invoke(Annotate()).GetProperty("columns");
        string outPath = Path.Combine(_temp, "window.TSV");
        JsonElement data = Invoke(Annotate("--offset", "100", "--limit", "3", "--out", outPath));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(3));
            Assert.That(data.GetProperty("offset").GetInt32(), Is.EqualTo(100));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
            Assert.That(Strings(data.GetProperty("columns"), "go_id"), Is.EqualTo(Strings(all, "go_id").Skip(100).Take(3)));
            Assert.That(File.ReadLines(outPath).Count(l => !l.StartsWith("#!")), Is.EqualTo(564), "header line + 563 rows");
        });
    }

    [Test]
    public void APastTheEndOffset_ReturnsNothing_AndSaysSo()
    {
        JsonElement data = Invoke(Annotate("--offset", "10000"));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("returned_count").GetInt32(), Is.EqualTo(0));
            Assert.That(data.GetProperty("truncated").GetBoolean(), Is.True);
        });
    }

    // ---- categories ---------------------------------------------------------------------------

    [Test]
    public void ACategoryMap_PlacesTermsUnderTheConsumersAnchors()
    {
        string categoriesOut = Path.Combine(_temp, "categories.tsv");
        JsonElement data = Invoke(Annotate("--category-map", Fixture("organelle_map.tsv"), "--categories-out", categoriesOut,
            "--limit", "0"));
        JsonElement categories = data.GetProperty("categories");
        JsonElement c = categories.GetProperty("columns");
        string[] goIds = Strings(c, "go_id");
        int inner = Array.IndexOf(goIds, "GO:0005743");

        Assert.Multiple(() =>
        {
            Assert.That(categories.GetProperty("map_name").GetString(), Is.EqualTo("organelle"));
            Assert.That(categories.GetProperty("anchor_count").GetInt32(), Is.EqualTo(9));
            Assert.That(categories.GetProperty("row_count").GetInt32(), Is.EqualTo(30));
            Assert.That(Strings(c, "category")[inner], Is.EqualTo("mitochondrion"));
            Assert.That(Strings(c, "subcategory")[inner], Is.EqualTo("mitochondrion:inner_membrane"));
            Assert.That(c.GetProperty("subcategory")[Array.IndexOf(goIds, "GO:0005739")].ValueKind, Is.EqualTo(JsonValueKind.Null),
                "the anchor of the category itself has no subcategory");
            Assert.That(File.ReadLines(categoriesOut).First(), Is.EqualTo("#!go_category_format 1"));
            Assert.That(File.ReadLines(categoriesOut).Count(l => !l.StartsWith("#!")), Is.EqualTo(31));
            Assert.That(data.GetProperty("categories_written").GetProperty("path").GetString(), Is.EqualTo(categoriesOut));
        });
    }

    [Test]
    public void NoCategoryMap_NoCategoryTable()
    {
        JsonElement data = Invoke(Annotate("--limit", "0"));

        Assert.That(data.GetProperty("categories").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    // ---- a database newer than the ontology ---------------------------------------------------

    [Test]
    public void AnIdTheReleaseLacks_IsRefused_UnlessSkipped_AndThenListed()
    {
        string database = DatabaseCiting("GO:9999999");

        JsonElement error = InvokeExpectingError("proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", database, "--go-obo", Fixture("go-pxd036557.obo"));
        JsonElement data = Invoke("proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", database, "--go-obo", Fixture("go-pxd036557.obo"), "--skip-unknown-go-ids", "--limit", "0");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("InvalidDataException"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("GO:9999999"));
            Assert.That(data.GetProperty("unresolved_go_ids").EnumerateArray().Select(x => x.GetString()), Is.EqualTo(new[] { "GO:9999999" }));
            Assert.That(data.GetProperty("skip_unknown_go_ids").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("header").GetProperty("unresolved_go_ids").GetString(), Is.EqualTo("1"));
            Assert.That(data.GetProperty("caveats")[0].GetString(), Does.StartWith("1 GO id(s) the database cites are absent"));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(563), "the invented term is dropped, nothing else moves");
        });
    }

    [Test]
    public void AGoOboWithNoRelease_SaysSo_AndTheRowsCarryNull()
    {
        string obo = Path.Combine(_temp, "unversioned.obo");
        File.WriteAllLines(obo, File.ReadLines(Fixture("go-pxd036557.obo")).Where(l => !l.StartsWith("data-version:")));

        JsonElement data = Invoke("proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", Fixture("pxd036557_proteins.xml"), "--go-obo", obo, "--limit", "1");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("go").GetProperty("release").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("columns").GetProperty("go_release")[0].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()),
                Has.Some.StartsWith("go.obo has no data-version header"));
        });
    }

    [Test]
    public void AnOutThatCannotBeWritten_FailsWithoutLeavingAPartialFile()
    {
        string outPath = Path.Combine(_temp, "taken.tsv");
        Directory.CreateDirectory(outPath);

        JsonElement error = InvokeExpectingError(Annotate("--out", outPath, "--limit", "0"));

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.Not.EqualTo("usage"));
            Assert.That(File.Exists(outPath + ".partial"), Is.False);
        });
    }

    // ---- usage, before any file is read -------------------------------------------------------

    [TestCase("go.csv")]
    [TestCase("go")]
    [TestCase("go.tsv.gz")]
    public void Out_MustBeTsv(string name)
    {
        JsonElement error = InvokeExpectingError(Annotate("--out", Path.Combine(_temp, name)));

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(".tsv"));
        });
    }

    [Test]
    public void CategoriesOut_MustBeTsv_AndNeedsAMap()
    {
        string noMap = InvokeExpectingError(Annotate("--categories-out", Path.Combine(_temp, "c.tsv"))).GetProperty("message").GetString()!;
        string notTsv = InvokeExpectingError(Annotate("--category-map", Fixture("organelle_map.tsv"),
            "--categories-out", Path.Combine(_temp, "c.txt"))).GetProperty("message").GetString()!;

        Assert.Multiple(() =>
        {
            Assert.That(noMap, Does.Contain("needs --category-map"));
            Assert.That(notTsv, Does.Contain("--categories-out must name a .tsv"));
        });
    }

    [Test]
    public void AMissingGoObo_SaysHowToFetchOne_AndNothingIsDownloaded()
    {
        string obo = Path.Combine(_temp, "go.obo");
        JsonElement error = InvokeExpectingError("proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", Fixture("pxd036557_proteins.xml"), "--go-obo", obo);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("proteins update-go"));
            Assert.That(File.Exists(obo), Is.False);
        });
    }

    [Test]
    public void AFastaDatabase_IsRefused_BecauseItCarriesNoGo()
    {
        JsonElement error = InvokeExpectingError("proteins", "annotate-go", "--groups", Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"),
            "--database", Fixture("human_extra.fasta"), "--go-obo", Fixture("go-pxd036557.obo"));

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("FASTA, which carries no GO terms"));
    }

    [Test]
    public void AnOutputOverAnInput_IsRefused()
    {
        string groups = Path.Combine(_temp, "groups.tsv");
        File.Copy(Fixture("PXD036557_AllQuantifiedProteinGroups.tsv"), groups);

        JsonElement error = InvokeExpectingError("proteins", "annotate-go", "--groups", groups,
            "--database", Fixture("pxd036557_proteins.xml"), "--go-obo", Fixture("go-pxd036557.obo"), "--out", groups);
        JsonElement same = InvokeExpectingError(Annotate("--category-map", Fixture("organelle_map.tsv"),
            "--out", Path.Combine(_temp, "x.tsv"), "--categories-out", Path.Combine(_temp, "X.tsv")));

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("names an input file"));
            Assert.That(same.GetProperty("message").GetString(), Does.Contain("name the same file"));
        });
    }

    [TestCase("--groups", "missing.tsv", "Protein-group table not found")]
    [TestCase("--database", "missing.xml", "Annotation database not found")]
    [TestCase("--category-map", "missing.tsv", "Category map not found")]
    public void AMissingInput_IsAUsageError(string option, string name, string message)
    {
        var args = Annotate().ToList();
        int at = args.IndexOf(option);
        if (at >= 0)
            args[at + 1] = Path.Combine(_temp, name);
        else
            args.AddRange(new[] { option, Path.Combine(_temp, name) });

        JsonElement error = InvokeExpectingError(args.ToArray());

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.StartWith(message));
        });
    }

    [TestCase("--limit", "-1")]
    [TestCase("--offset", "-1")]
    public void ANegativeWindow_IsAUsageError(string option, string value)
    {
        JsonElement error = InvokeExpectingError(Annotate(option, value));

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("zero or greater"));
    }

    [Test]
    public void AnOutFolderThatDoesNotExist_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Annotate("--out", Path.Combine(_temp, "nowhere", "go.tsv")));

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("does not exist"));
    }

    [Test]
    public void UpdateGo_IntoAFolderThatDoesNotExist_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError("proteins", "update-go", "--go-obo", Path.Combine(_temp, "nowhere", "go.obo"));

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
    }

    // ---- helpers ------------------------------------------------------------------------------

    /// <summary>The fixture database with one GO reference the ontology release does not have.</summary>
    private string DatabaseCiting(string goId)
    {
        string xml = File.ReadAllText(Fixture("pxd036557_proteins.xml"));
        int at = xml.IndexOf("<dbReference type=\"GO\"", StringComparison.Ordinal);
        string injected = $"<dbReference type=\"GO\" id=\"{goId}\"><property type=\"term\" value=\"C:invented\"/>" +
                          "<property type=\"evidence\" value=\"ECO:0000314\"/></dbReference>\n";
        string path = Path.Combine(_temp, "newer.xml");
        File.WriteAllText(path, xml.Insert(at, injected), new UTF8Encoding(false));
        return path;
    }

    private static int Row(JsonElement columns, string group, string goId)
    {
        string[] groups = Strings(columns, "protein_group");
        string[] ids = Strings(columns, "go_id");
        return Enumerable.Range(0, groups.Length).Single(i => groups[i] == group && ids[i] == goId);
    }

    private static string Fixture(string name, [CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "python", "tests", "fixtures", "proteins", name));

    private static string[] Strings(JsonElement columns, string name) =>
        columns.GetProperty(name).EnumerateArray().Select(x => x.ValueKind == JsonValueKind.Null ? null! : x.GetString()!).ToArray();

    private static int[] Ints(JsonElement columns, string name) =>
        columns.GetProperty(name).EnumerateArray().Select(x => x.GetInt32()).ToArray();

    private static string[] StringList(JsonElement columns, string name, int row) =>
        columns.GetProperty(name)[row].EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) => Unwrap(Envelope(args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk), $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    private static JsonElement Envelope(string[] args)
    {
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
        catch (Exception exception)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = Program.ClassifyError(exception), message = Program.Unwrap(exception).Message } },
                Program.JsonOptions);
        }
    }
}

/// <summary>
/// <c>proteins update-go</c> against GO's live PURL. Skipped, never failed, when it is unreachable.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("GeneOntology")]
[ExcludeFromCodeCoverage]
public class GoLiveCanaryTests
{
    [Test]
    public Task UpdateGo_FetchesARelease_AndASecondFetchChangesNothing() =>
        ExternalServiceTestHelper.RunAsync("Gene Ontology PURL", async () =>
        {
            string folder = Path.Combine(Path.GetTempPath(), $"pymzlib-go-live-{Guid.NewGuid():N}");
            Directory.CreateDirectory(folder);
            try
            {
                string obo = Path.Combine(folder, "go.obo");
                JsonElement first = await InvokeAsync("proteins", "update-go", "--go-obo", obo);
                JsonElement second = await InvokeAsync("proteins", "update-go", "--go-obo", obo);

                Assert.Multiple(() =>
                {
                    Assert.That(first.GetProperty("existed_before").GetBoolean(), Is.False);
                    Assert.That(first.GetProperty("go").GetProperty("term_count").GetInt32(), Is.GreaterThan(40000));
                    Assert.That(first.GetProperty("go").GetProperty("release").GetString(), Does.StartWith("releases/"));
                    Assert.That(second.GetProperty("existed_before").GetBoolean(), Is.True);
                    // GO can publish between the two calls; then the second must say so and keep a backup.
                    if (second.GetProperty("changed").GetBoolean())
                        Assert.That(Directory.GetFiles(folder, "go.obo.*"), Is.Not.Empty);
                    else
                        Assert.That(second.GetProperty("previous_sha256").GetString(),
                            Is.EqualTo(second.GetProperty("go").GetProperty("sha256").GetString()));
                });
            }
            finally
            {
                Directory.Delete(folder, recursive: true);
            }
        });

    private static async Task<JsonElement> InvokeAsync(params string[] args)
    {
        try
        {
            object data = await Program.DispatchAsync(args);
            return JsonSerializer.SerializeToElement(data, Program.JsonOptions);
        }
        catch (Exception exception) when (Program.ClassifyError(exception) == Program.ServiceUnavailableType)
        {
            throw new ExternalServiceUnavailableException(Program.Unwrap(exception).Message);
        }
    }
}
