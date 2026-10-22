using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Text.Json;
using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>genes orthology</c>, the projection of mzLib's <see cref="OrthologySnapshot.Build"/>,
/// <see cref="OrthologySnapshot.PairStatus"/> and <see cref="OrthologySnapshotWriter.WriteAsync"/>
/// (mzLib #1381).
/// </summary>
/// <remarks>
/// The inputs are mzLib's own test world (Test/FileReadingTests/TestOrthologySnapshot.cs at the pin),
/// committed in <c>pkg/python/tests/fixtures/orthology/</c> with each dump where Ensembl publishes it,
/// <c>homologies/&lt;species&gt;/Compara.116.protein_default.homologies.tsv.gz</c>. It holds both traps
/// the builder exists for: every human-mouse row is in the MOUSE dump, and H2~M2 and M2~R2 are
/// orthologs while H2~R2 is not called. The expected values are mzLib's test's own.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class GenesOrthologyTests
{
    private const string Human = "homo_sapiens", Mouse = "mus_musculus", Rat = "rattus_norvegicus";
    private const string Trees = "vertebrates.GeneTree_content.default.e116.txt.gz";
    private const string DumpName = "Compara.116.protein_default.homologies.tsv.gz";

    private string _scratch = "";

    [SetUp]
    public void SetUp()
    {
        _scratch = Path.Combine(Path.GetTempPath(), "genes-orthology-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_scratch);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_scratch))
            Directory.Delete(_scratch, recursive: true);
    }

    // ---- the snapshot ---------------------------------------------------------------------------

    [Test]
    public void PairStatus_SaysHowFarTheInferenceGot_AsMzLibsTestDoes()
    {
        JsonElement status = Invoke(Stdin(), Args("--pair-status", $"{Human},{Mouse}")).GetProperty("pair_status");
        JsonElement columns = status.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(status.GetProperty("from_species").GetString(), Is.EqualTo(Human));
            Assert.That(status.GetProperty("to_species").GetString(), Is.EqualTo(Mouse));
            Assert.That(status.GetProperty("row_count").GetInt32(), Is.EqualTo(4));
            Assert.That(Strings(status.GetProperty("column_names")), Is.EqualTo(new[] { "gene_id", "gene_name", "gene_biotype", "status" }));
            Assert.That(Strings(columns.GetProperty("gene_id")), Is.EqualTo(new[]
            {
                "ENSG00000000001", "ENSG00000000002", "ENSG00000000003", "ENSG00000000004",
            }), "every gene of the human gene set, in ordinal order");
            Assert.That(Strings(columns.GetProperty("status")), Is.EqualTo(new[]
            {
                "has_ortholog", "has_ortholog", "tree_lacks_target_species", "not_in_any_tree",
            }), "mzLib's PairStatus_SaysHowFarTheInferenceGot");
            Assert.That(Strings(columns.GetProperty("gene_name")), Is.All.Null, "the GTFs give no gene_name");
            Assert.That(Strings(columns.GetProperty("gene_biotype")), Is.All.EqualTo("protein_coding"));
            Assert.That(status.GetProperty("status_counts").EnumerateObject().Select(p => (p.Name, p.Value.GetInt32())), Is.EqualTo(new[]
            {
                ("has_ortholog", 2), ("no_edge_in_shared_tree", 0), ("tree_lacks_target_species", 1),
                ("not_in_any_tree", 1), ("not_in_gene_set", 0),
            }), "every status in mzLib's order, zeros included");
        });
    }

    [Test]
    public void PairStatus_NoEdgeInSharedTree_IsHumanToRatsH2()
    {
        // ENSGT2 holds R2, and no H2-R2 ortholog was called (mzLib's test: StatusOf(H2, Rat)).
        JsonElement columns = Invoke(Stdin(), Args("--pair-status", $"{Human},{Rat}"))
            .GetProperty("pair_status").GetProperty("columns");

        Assert.That(Strings(columns.GetProperty("status")), Is.EqualTo(new[]
        {
            "has_ortholog", "no_edge_in_shared_tree", "tree_lacks_target_species", "not_in_any_tree",
        }));
    }

    [Test]
    public void Envelope_ReadsEveryDump_BecauseAPairsRowsSitInEitherOne()
    {
        JsonElement data = Invoke(Stdin(), Args());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(Strings(data.GetProperty("species")), Is.EqualTo(new[] { Human, Mouse, Rat }));
            Assert.That(data.GetProperty("format").GetString(), Is.EqualTo(OrthologySnapshotWriter.FormatName));
            Assert.That(data.GetProperty("format_version").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("identical_duplicates_dropped").GetInt32(), Is.Zero);
            Assert.That(data.GetProperty("pairs").EnumerateArray().Select(p => (
                    p.GetProperty("species_a").GetString(), p.GetProperty("species_b").GetString(),
                    p.GetProperty("row_count").GetInt32(), p.GetProperty("ortholog_count").GetInt32(),
                    p.GetProperty("paralog_count").GetInt32(), p.GetProperty("homoeolog_count").GetInt32())),
                Is.EqualTo(new[]
                {
                    (Human, Human, 0, 0, 0, 0), (Human, Mouse, 3, 3, 0, 0), (Human, Rat, 1, 1, 0, 0),
                    (Mouse, Mouse, 1, 0, 1, 0), (Mouse, Rat, 2, 2, 0, 0), (Rat, Rat, 0, 0, 0, 0),
                }), "the three human-mouse rows all came from the mouse dump");
            Assert.That(data.GetProperty("pair_status").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void Envelope_NestedObjects_HaveTheSpecsKeys()
    {
        JsonElement data = Invoke(Stdin(), Args("--pair-status", $"{Mouse},{Human}", "--out", Path.Combine(_scratch, "s")));

        Assert.Multiple(() =>
        {
            Assert.That(Keys(data), Is.EqualTo(new[]
            {
                "release", "species", "format", "format_version", "snapshot_id", "identical_duplicates_dropped",
                "gene_sets", "gene_trees", "dumps", "pairs", "written", "pair_status", "caveats",
            }));
            Assert.That(Keys(data.GetProperty("gene_sets")), Is.EqualTo(new[] { Human, Mouse, Rat }));
            Assert.That(Keys(data.GetProperty("gene_sets").GetProperty(Rat)), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "genome_build", "genebuild_last_updated", "gene_count",
            }), "genes resolve's gene_set fields");
            Assert.That(Keys(data.GetProperty("gene_trees")), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "collection", "member_count", "restricted_to_gene_sets",
            }));
            Assert.That(Keys(data.GetProperty("dumps")), Is.EqualTo(new[] { Human, Mouse, Rat }));
            Assert.That(Keys(data.GetProperty("dumps").GetProperty(Human)), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "rows_read", "rows_kept",
            }));
            Assert.That(Keys(data.GetProperty("pairs")[0]), Is.EqualTo(new[]
            {
                "species_a", "species_b", "row_count", "ortholog_count", "paralog_count", "homoeolog_count",
            }));
            Assert.That(Keys(data.GetProperty("written")), Is.EqualTo(new[] { "directory", "files" }));
            Assert.That(Keys(data.GetProperty("written").GetProperty("files")[0]), Is.EqualTo(new[] { "path", "rows", "bytes", "sha256" }));
            Assert.That(Keys(data.GetProperty("pair_status")), Is.EqualTo(new[]
            {
                "from_species", "to_species", "row_count", "status_counts", "column_names", "columns",
            }));
        });
    }

    [Test]
    public void Inputs_CarryMzLibsProvenance()
    {
        JsonElement data = Invoke(Stdin(), Args());
        JsonElement human = data.GetProperty("gene_sets").GetProperty(Human);
        JsonElement trees = data.GetProperty("gene_trees");
        JsonElement mouseDump = data.GetProperty("dumps").GetProperty(Mouse);

        Assert.Multiple(() =>
        {
            Assert.That(human.GetProperty("source_file_name").GetString(), Is.EqualTo("Homo_sapiens.GRCh38.116.gtf"));
            Assert.That(human.GetProperty("sha256").GetString(), Is.EqualTo(Sha256(Fixture("Homo_sapiens.GRCh38.116.gtf"))));
            Assert.That(human.GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(human.GetProperty("genome_build").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(human.GetProperty("gene_count").GetInt32(), Is.EqualTo(4));

            Assert.That(trees.GetProperty("source_file_name").GetString(), Is.EqualTo(Trees));
            Assert.That(trees.GetProperty("sha256").GetString(), Is.EqualTo(Sha256(Fixture(Trees))), "of the compressed bytes");
            Assert.That(trees.GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(trees.GetProperty("collection").GetString(), Is.EqualTo("vertebrates"));
            Assert.That(trees.GetProperty("member_count").GetInt32(), Is.EqualTo(9), "unrestricted: the zebrafish gene too");
            Assert.That(trees.GetProperty("restricted_to_gene_sets").ValueKind, Is.EqualTo(JsonValueKind.Null));

            Assert.That(mouseDump.GetProperty("source_file_name").GetString(), Is.EqualTo(DumpName));
            Assert.That(mouseDump.GetProperty("sha256").GetString(), Is.EqualTo(Sha256(Fixture($"homologies/{Mouse}/{DumpName}"))));
            Assert.That(mouseDump.GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(mouseDump.GetProperty("rows_read").GetInt32(), Is.EqualTo(5));
            Assert.That(mouseDump.GetProperty("rows_kept").GetInt32(), Is.EqualTo(5));
        });
    }

    [Test]
    public void SnapshotId_IsMzLibs_ForTheSameInputs()
    {
        string id = Invoke(Stdin(), Args()).GetProperty("snapshot_id").GetString()!;

        var sets = new Dictionary<string, EnsemblGeneSet>
        {
            [Human] = EnsemblGeneSet.LoadGtf(Fixture("Homo_sapiens.GRCh38.116.gtf")),
            [Mouse] = EnsemblGeneSet.LoadGtf(Fixture("Mus_musculus.GRCm39.116.gtf")),
            [Rat] = EnsemblGeneSet.LoadGtf(Fixture("Rattus_norvegicus.GRCr8.116.gtf")),
        };
        var dumps = new[] { Human, Mouse, Rat }
            .Select(s => ComparaHomologyDump.Load(Fixture($"homologies/{s}/{DumpName}"), sets.Keys)).ToList();
        var direct = OrthologySnapshot.Build("116", sets, ComparaGeneTreeContent.Load(Fixture(Trees)), dumps);

        Assert.That(id, Is.EqualTo(OrthologySnapshotWriter.SnapshotId(direct)));
        Assert.That(id, Does.Match("^[0-9a-f]{64}$"));
    }

    [Test]
    public void DumpOrderOnStdin_DoesNotChangeTheSnapshot()
    {
        JsonElement given = Invoke(Stdin(), Args());
        JsonElement reversed = Invoke(Stdin(dumpOrder: new[] { Rat, Human, Mouse }), Args());

        Assert.That(reversed.GetProperty("snapshot_id").GetString(), Is.EqualTo(given.GetProperty("snapshot_id").GetString()));
        Assert.That(Keys(reversed.GetProperty("dumps")), Is.EqualTo(new[] { Human, Mouse, Rat }), "in species order, as mzLib holds them");
    }

    // ---- --restrict-gene-trees ------------------------------------------------------------------

    [Test]
    public void RestrictGeneTrees_KeepsOnlyTheGivenGeneSets_AndChangesNothingElse()
    {
        string unrestrictedDir = Path.Combine(_scratch, "all"), restrictedDir = Path.Combine(_scratch, "restricted");
        JsonElement all = Invoke(Stdin(), Args("--out", unrestrictedDir, "--pair-status", $"{Human},{Rat}"));
        JsonElement restricted = Invoke(Stdin(), Args("--restrict-gene-trees", "--out", restrictedDir, "--pair-status", $"{Human},{Rat}"));
        JsonElement trees = restricted.GetProperty("gene_trees");

        Assert.Multiple(() =>
        {
            Assert.That(trees.GetProperty("member_count").GetInt32(), Is.EqualTo(8), "the zebrafish gene is not kept");
            Assert.That(Strings(trees.GetProperty("restricted_to_gene_sets")), Is.EqualTo(new[] { Human, Mouse, Rat }
                .Select(s => restricted.GetProperty("gene_sets").GetProperty(s).GetProperty("sha256").GetString())
                .OrderBy(s => s, StringComparer.Ordinal)));
            Assert.That(restricted.GetProperty("snapshot_id").GetString(), Is.EqualTo(all.GetProperty("snapshot_id").GetString()));
            Assert.That(Shas(restricted), Is.EqualTo(Shas(all)), "every written byte is the same");
            Assert.That(File.ReadAllBytes(Path.Combine(restrictedDir, "manifest.json")),
                Is.EqualTo(File.ReadAllBytes(Path.Combine(unrestrictedDir, "manifest.json"))));
            Assert.That(restricted.GetProperty("pair_status").GetRawText(), Is.EqualTo(all.GetProperty("pair_status").GetRawText()));
            Assert.That(Strings(all.GetProperty("caveats")), Has.Some.Contains("read unrestricted: 9 genes were held in memory, 8 of them"));
            Assert.That(Strings(restricted.GetProperty("caveats")), Has.None.Contains("unrestricted"));
        });
    }

    [Test]
    public void RestrictGeneTrees_WithAValue_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(), Args("--restrict-gene-trees", "true"));

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("--restrict-gene-trees takes no value"));
    }

    // ---- --out ----------------------------------------------------------------------------------

    [Test]
    public void Out_WritesMzLibsLayout_AndReportsEveryFileButTheManifest()
    {
        string directory = Path.Combine(_scratch, "compara-116");
        JsonElement data = Invoke(Stdin(), Args("--out", directory));
        JsonElement written = data.GetProperty("written");
        var files = written.GetProperty("files").EnumerateArray().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(written.GetProperty("directory").GetString(), Is.EqualTo(Path.GetFullPath(directory)));
            Assert.That(files.Select(f => f.GetProperty("path").GetString()), Is.EqualTo(new[]
            {
                "genes/homo_sapiens.parquet", "genes/mus_musculus.parquet", "genes/rattus_norvegicus.parquet",
                "members/homo_sapiens.parquet", "members/mus_musculus.parquet", "members/rattus_norvegicus.parquet",
                "pairs/homo_sapiens__homo_sapiens.parquet", "pairs/homo_sapiens__mus_musculus.parquet",
                "pairs/homo_sapiens__rattus_norvegicus.parquet", "pairs/mus_musculus__mus_musculus.parquet",
                "pairs/mus_musculus__rattus_norvegicus.parquet", "pairs/rattus_norvegicus__rattus_norvegicus.parquet",
                "views.sql",
            }), "mzLib's Write_LayoutAndManifest");
            Assert.That(Rows(files, "members/homo_sapiens.parquet"), Is.EqualTo(3), "H4 is in no tree");
            Assert.That(Rows(files, "pairs/homo_sapiens__mus_musculus.parquet"), Is.EqualTo(3));
            Assert.That(Rows(files, "pairs/homo_sapiens__homo_sapiens.parquet"), Is.Zero, "an empty pair is still written");
            Assert.That(files.Single(f => f.GetProperty("path").GetString() == "views.sql").GetProperty("rows").ValueKind,
                Is.EqualTo(JsonValueKind.Null), "views.sql is not a table");
            foreach (JsonElement f in files)
            {
                string path = Path.Combine(directory, f.GetProperty("path").GetString()!);
                Assert.That(f.GetProperty("sha256").GetString(), Is.EqualTo(Sha256(path)), path);
                Assert.That(f.GetProperty("bytes").GetInt64(), Is.EqualTo(new FileInfo(path).Length), path);
            }

            using JsonDocument manifest = JsonDocument.Parse(File.ReadAllText(Path.Combine(directory, "manifest.json")));
            Assert.That(manifest.RootElement.GetProperty("snapshot_id").GetString(), Is.EqualTo(data.GetProperty("snapshot_id").GetString()));
            Assert.That(manifest.RootElement.GetProperty("files").GetArrayLength(), Is.EqualTo(files.Count));
            Assert.That(manifest.RootElement.GetProperty("gene_set_sha256").GetProperty(Rat).GetString(),
                Is.EqualTo(data.GetProperty("gene_sets").GetProperty(Rat).GetProperty("sha256").GetString()));
        });
    }

    [Test]
    public void Out_TheSameInputsGiveTheSameBytes()
    {
        JsonElement first = Invoke(Stdin(), Args("--out", Path.Combine(_scratch, "a")));
        JsonElement second = Invoke(Stdin(), Args("--out", Path.Combine(_scratch, "b")));

        Assert.That(Shas(second), Is.EqualTo(Shas(first)));
        Assert.That(File.ReadAllBytes(Path.Combine(_scratch, "b", "manifest.json")),
            Is.EqualTo(File.ReadAllBytes(Path.Combine(_scratch, "a", "manifest.json"))));
    }

    [Test]
    public void Out_AnEmptyExistingDirectory_IsWrittenInto()
    {
        string directory = Path.Combine(_scratch, "empty");
        Directory.CreateDirectory(directory);

        Assert.That(Invoke(Stdin(), Args("--out", directory)).GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Object));
    }

    [Test]
    public void Out_ADirectoryThatIsNotEmpty_IsRefusedBeforeAnyInputIsRead()
    {
        string directory = Path.Combine(_scratch, "used");
        Directory.CreateDirectory(directory);
        File.WriteAllText(Path.Combine(directory, "old.txt"), "x");
        // A gene-trees file mzLib would refuse: reaching it would be a correctness failure, not usage.
        string badTrees = Write("vertebrates.GeneTree_content.default.e116.txt", "ENSGT1\tP1\tG1\n");

        JsonElement error = InvokeExpectingError(Stdin(), Args("--out", directory), trees: badTrees);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("is not empty"));
        Assert.That(Directory.EnumerateFileSystemEntries(directory), Has.Exactly(1).Items, "nothing was written");
    }

    [Test]
    public void Out_AFile_IsAUsageError()
    {
        string file = Write("snapshot.tsv", "x");

        JsonElement error = InvokeExpectingError(Stdin(), Args("--out", file));

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("is a file"));
    }

    // ---- caveats --------------------------------------------------------------------------------

    [Test]
    public void Caveats_EnsemblsSharedDumpName_IsNamed()
    {
        string?[] caveats = Strings(Invoke(Stdin(), Args()).GetProperty("caveats"));

        Assert.That(caveats, Has.Some.Contains($"share the file name {DumpName}"));
        Assert.That(caveats, Has.None.Contains("no release in the file name"));
    }

    [Test]
    public void Caveats_AFileNameWithoutARelease_IsNotHeldToIt()
    {
        string renamed = Copy(Fixture($"homologies/{Human}/{DumpName}"), "human.homologies.tsv.gz");

        JsonElement data = Invoke(Stdin(dumpPath: s => s == Human ? renamed : null), Args());

        Assert.That(data.GetProperty("dumps").GetProperty(Human).GetProperty("release").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Strings(data.GetProperty("caveats")), Has.Some.Contains(
            "human.homologies.tsv.gz (the homo_sapiens dump) carries no release in the file name, so mzLib could not hold it to --release 116"));
    }

    [Test]
    public void Caveats_AGtfNamedForItsOwnRelease_IsAccepted_AndSaidSo()
    {
        string gtf = Copy(Fixture("Rattus_norvegicus.GRCr8.116.gtf"), "Rattus_norvegicus.GRCr8.63.gtf");

        JsonElement data = Invoke(Stdin(gtfPath: s => s == Rat ? gtf : null), Args());

        Assert.That(data.GetProperty("gene_sets").GetProperty(Rat).GetProperty("release").GetString(), Is.EqualTo("63"));
        Assert.That(Strings(data.GetProperty("caveats")), Has.Some.Contains("is named for release 63, not 116"));
    }

    // ---- a compact gene table ---------------------------------------------------------------------

    [Test]
    public void AGeneSetTable_KeepsTheGtfsProvenance_SoTheSnapshotIsTheSame()
    {
        string table = Path.Combine(_scratch, "mouse.genes.tsv");
        EnsemblGeneSetWriter.Write(table, EnsemblGeneSet.LoadGtf(Fixture("Mus_musculus.GRCm39.116.gtf")));

        JsonElement fromGtf = Invoke(Stdin(), Args());
        JsonElement fromTable = Invoke(Stdin(kind: s => s == Mouse ? "gene-set" : "gtf", gtfPath: s => s == Mouse ? table : null), Args());

        Assert.That(fromTable.GetProperty("snapshot_id").GetString(), Is.EqualTo(fromGtf.GetProperty("snapshot_id").GetString()));
        Assert.That(fromTable.GetProperty("gene_sets").GetProperty(Mouse).GetRawText(),
            Is.EqualTo(fromGtf.GetProperty("gene_sets").GetProperty(Mouse).GetRawText()));
    }

    // ---- mzLib's refusals -----------------------------------------------------------------------

    [Test]
    public void AMissingDump_IsMzLibsRefusal_AsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(dumpOrder: new[] { Human, Rat }), Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo(
            "No homology dump for mus_musculus. Each pair's relationships are split between its two species' dumps, " +
            "so every species' dump is required."));
    }

    [Test]
    public void AnotherRelease_IsMzLibsRefusal_AsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(), ArgsFor("117", null));

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain($"{Trees} names release 116, not 117"));
    }

    [Test]
    public void TwoDumpsOfOneSpecies_IsMzLibsRefusal_AsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(dumpOrder: new[] { Human, Mouse, Mouse, Rat }), Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("Two dumps of mus_musculus"));
    }

    [Test]
    public void ADumpOfASpeciesNotGiven_IsMzLibsRefusal_AsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(species: new[] { Human, Mouse }), Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("which is not one of the species"));
    }

    [Test]
    public void AMalformedDump_IsACorrectnessFailure_UnderMzLibsType()
    {
        string bad = Write("Compara.116.protein_default.homologies.tsv", "gene_stable_id\tspecies\n");

        JsonElement error = InvokeExpectingError(Stdin(dumpPath: s => s == Rat ? bad : null), Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("InvalidDataException"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("unexpected columns"));
    }

    [Test]
    public void AHomologyNamingAGeneOutsideItsGeneSet_IsACorrectnessFailure()
    {
        // The human gene set without ENSG00000000002, which two mouse rows name: a gene set and dumps
        // from different releases look exactly like this.
        string gtf = Write("Homo_sapiens.GRCh38.116.gtf", string.Concat(
            File.ReadAllLines(Fixture("Homo_sapiens.GRCh38.116.gtf")).Where(l => !l.Contains("ENSG00000000002")).Select(l => l + "\n")));

        JsonElement error = InvokeExpectingError(Stdin(gtfPath: s => s == Human ? gtf : null), Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("InvalidDataException"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("names ENSG00000000002, which is not in the homo_sapiens gene set"));
    }

    // ---- the bridge's own input -----------------------------------------------------------------

    private static readonly object[] MalformedStdin =
    {
        new object[] { "homo_sapiens\tgtf\tx.gtf\n", "No such line was found" },
        new object[] { "--\nx.tsv.gz\n", "No species lines before the '--' line" },
        new object[] { "homo_sapiens\tgtf\tx.gtf\n--\n", "No homology dumps after the '--' line" },
        new object[] { "homo_sapiens\tx.gtf\n--\nx.tsv.gz\n", "is not three tab-separated fields" },
        new object[] { "homo_sapiens\tfasta\tx.gtf\n--\nx.tsv.gz\n", "names its gene set as 'fasta'" },
        new object[] { "homo_sapiens\tgtf\tx.gtf\nhomo_sapiens\tgtf\ty.gtf\n--\nx.tsv.gz\n", "is listed twice" },
        new object[] { "homo_sapiens\tgtf\tabsent.gtf\n--\nx.tsv.gz\n", "GTF of homo_sapiens not found" },
    };

    [TestCaseSource(nameof(MalformedStdin))]
    public void MalformedStdin_IsAUsageError(string stdin, string message)
    {
        JsonElement error = InvokeExpectingError(stdin, Args());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(message));
    }

    [TestCase("homo_sapiens", "takes two species separated by a comma")]
    [TestCase("homo_sapiens,", "takes two species separated by a comma")]
    [TestCase("homo_sapiens,danio_rerio", "names 'danio_rerio', which is not one of the species on stdin")]
    [TestCase("mus_musculus,mus_musculus", "names one species twice")]
    public void PairStatus_NotTwoGivenSpecies_IsRefusedBeforeAnyInputIsRead(string value, string message)
    {
        string badTrees = Write("vertebrates.GeneTree_content.default.e116.txt", "ENSGT1\tP1\tG1\n");

        JsonElement error = InvokeExpectingError(Stdin(), Args("--pair-status", value), trees: badTrees);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(message));
    }

    [TestCase("--release")]
    [TestCase("--gene-trees")]
    public void ARequiredOption_LeftOut_IsAUsageError(string option)
    {
        var args = Args().ToList();
        int at = args.IndexOf(option);
        args.RemoveRange(at, 2);

        JsonElement error = InvokeExpectingError(Stdin(), args.ToArray());

        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo($"Missing required option {option}."));
    }

    [TestCase("--out")]
    [TestCase("--pair-status")]
    public void AValuedOption_WithoutAValue_IsAUsageError(string option)
    {
        JsonElement error = InvokeExpectingError(Stdin(), Args(option));

        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo($"Option {option} was given without a value."));
    }

    [Test]
    public void AMissingGeneTreesFile_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(), Args(), trees: Path.Combine(_scratch, "absent.txt.gz"));

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("Gene-trees file not found"));
    }

    // ---- plumbing -------------------------------------------------------------------------------

    private static readonly Dictionary<string, string> GtfOf = new()
    {
        [Human] = "Homo_sapiens.GRCh38.116.gtf",
        [Mouse] = "Mus_musculus.GRCm39.116.gtf",
        [Rat] = "Rattus_norvegicus.GRCr8.116.gtf",
    };

    /// <summary>The stdin of a call: species lines, '--', dump paths. Each override returns null to keep the fixture.</summary>
    private static string Stdin(string[]? species = null, string[]? dumpOrder = null, Func<string, string>? kind = null,
        Func<string, string?>? gtfPath = null, Func<string, string?>? dumpPath = null)
    {
        species ??= new[] { Human, Mouse, Rat };
        dumpOrder ??= new[] { Human, Mouse, Rat };
        var lines = species.Select(s =>
            $"{s}\t{kind?.Invoke(s) ?? "gtf"}\t{gtfPath?.Invoke(s) ?? Fixture(GtfOf[s])}").ToList();
        lines.Add("--");
        lines.AddRange(dumpOrder.Select(s => dumpPath?.Invoke(s) ?? Fixture($"homologies/{s}/{DumpName}")));
        return string.Join("\n", lines) + "\n";
    }

    private string[] Args(params string[] extra) => ArgsFor("116", null, extra);

    private string[] ArgsFor(string release, string? trees, params string[] extra) =>
        new[] { "genes", "orthology", "--release", release, "--gene-trees", trees ?? Fixture(Trees) }.Concat(extra).ToArray();

    private JsonElement InvokeExpectingError(string stdin, string[] args, string trees) =>
        InvokeExpectingError(stdin, ReplaceTrees(args, trees));

    private static string[] ReplaceTrees(string[] args, string trees)
    {
        string[] copy = args.ToArray();
        copy[Array.IndexOf(copy, "--gene-trees") + 1] = trees;
        return copy;
    }

    private string Write(string name, string text)
    {
        string path = Path.Combine(_scratch, name);
        File.WriteAllText(path, text);
        return path;
    }

    private string Copy(string source, string name)
    {
        string path = Path.Combine(_scratch, name);
        File.Copy(source, path);
        return path;
    }

    /// <summary>A file of pkg/python/tests/fixtures/orthology, committed beside these tests.</summary>
    private static string Fixture(string name, [CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "python", "tests", "fixtures", "orthology", name));

    private static string Sha256(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

    private static string?[] Strings(JsonElement array) =>
        array.EnumerateArray().Select(e => e.ValueKind == JsonValueKind.Null ? null : e.GetString()).ToArray();

    private static string[] Keys(JsonElement element) => element.EnumerateObject().Select(p => p.Name).ToArray();

    private static long? Rows(List<JsonElement> files, string path)
    {
        JsonElement rows = files.Single(f => f.GetProperty("path").GetString() == path).GetProperty("rows");
        return rows.ValueKind == JsonValueKind.Null ? null : rows.GetInt64();
    }

    private static string[] Shas(JsonElement data) =>
        data.GetProperty("written").GetProperty("files").EnumerateArray()
            .Select(f => f.GetProperty("path").GetString() + " " + f.GetProperty("sha256").GetString()).ToArray();

    private static JsonElement Invoke(string stdin, string[] args)
    {
        JsonElement envelope = JsonDocument.Parse(Raw(stdin, args)).RootElement;
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.True, $"Expected success, got: {envelope}");
        return envelope.GetProperty("data").Clone();
    }

    private static JsonElement InvokeExpectingError(string stdin, string[] args)
    {
        JsonElement envelope = JsonDocument.Parse(Raw(stdin, args)).RootElement;
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.False, $"Expected a failure, got: {envelope}");
        return envelope.GetProperty("error").Clone();
    }

    /// <summary>Runs a verb as Main would and returns the envelope text exactly as it would be written.</summary>
    private static string Raw(string stdin, string[] args)
    {
        TextReader previousIn = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.Serialize(new Program.Envelope { Ok = true, Data = data }, Program.JsonOptions);
        }
        catch (Exception exception)
        {
            string type = Program.Unwrap(exception) is Program.UsageException ? "usage" : Program.ClassifyError(exception);
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
