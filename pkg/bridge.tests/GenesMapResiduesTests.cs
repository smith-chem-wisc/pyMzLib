using System.Diagnostics.CodeAnalysis;
using System.IO.Compression;
using System.Runtime.CompilerServices;
using System.Text.Json;
using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>genes map-residues</c>, the projection of mzLib's <see cref="ComparaGeneTreeAlignment"/>
/// (<see cref="ComparaGeneTreeAlignment.Load"/> and <see cref="ComparaGeneTreeAlignment.MapResidue"/>,
/// mzLib #1465).
/// </summary>
/// <remarks>
/// The inputs are a real cut of Ensembl release 116, committed in
/// <c>pkg/python/tests/fixtures/map_residues/</c> (PROVENANCE.txt there): the TP53 and GAPDH gene trees'
/// alignments with the human, mouse and rat canonical rows verbatim, their gene-tree rows, and the
/// genes' GTF rows. Every expected row is mzLib's own answer: the tests load the same files with mzLib
/// and call <see cref="ComparaGeneTreeAlignment.MapResidue"/> directly.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class GenesMapResiduesTests
{
    private const string Human = "homo_sapiens", Mouse = "mus_musculus", Rat = "rattus_norvegicus";
    private const string AlignmentName = "Compara.116.protein_default.aa.fasta.gz";
    private const string TreesName = "vertebrates.GeneTree_content.default.e116.txt.gz";
    private const string Tp53 = "ENSP00000269305", Trp53 = "ENSMUSP00000104298", RatTp53 = "ENSRNOP00000074031";
    private const string MouseGapdh = "ENSMUSP00000113213", Tp53Isoform = "ENSP00000352610", Trp53Isoform = "ENSMUSP00000005371";

    private static readonly (string Source, int Position, string Target)[] Queries =
    {
        (Tp53, 15, Trp53), (Tp53, 20, Trp53), (Tp53, 46, Trp53), (Tp53, 392, Trp53), (Trp53, 18, Tp53),
        (Tp53, 15, RatTp53), (Tp53, 15, MouseGapdh), (Tp53Isoform, 15, Trp53), (Tp53, 15, Trp53Isoform),
    };

    private string _scratch = "";

    [SetUp]
    public void SetUp()
    {
        _scratch = Path.Combine(Path.GetTempPath(), "genes-map-residues-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_scratch);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_scratch))
            Directory.Delete(_scratch, recursive: true);
    }

    // ---- the answers are mzLib's ----------------------------------------------------------------

    [Test]
    public void EveryRow_IsMzLibsOwnMapResidue_InInputOrder()
    {
        JsonElement data = Invoke(Stdin(Queries), Args());
        ComparaGeneTreeAlignment mzlib = LoadWithMzLib();
        JsonElement columns = data.GetProperty("columns");

        for (int i = 0; i < Queries.Length; i++)
        {
            var (source, position, target) = Queries[i];
            ComparaColumnMapping expected = mzlib.MapResidue(source, position, target);
            Assert.Multiple(() =>
            {
                Assert.That(columns.GetProperty("source_protein_id")[i].GetString(), Is.EqualTo(source));
                Assert.That(columns.GetProperty("position")[i].GetInt32(), Is.EqualTo(position));
                Assert.That(columns.GetProperty("target_protein_id")[i].GetString(), Is.EqualTo(target));
                Assert.That(columns.GetProperty("outcome")[i].GetString(), Is.EqualTo(expected.Outcome.ToString()));
                Assert.That(NullableInt(columns.GetProperty("column")[i]), Is.EqualTo(expected.Column));
                Assert.That(NullableInt(columns.GetProperty("target_position")[i]), Is.EqualTo(expected.TargetPosition));
                Assert.That(NullableString(columns.GetProperty("target_residue")[i]), Is.EqualTo(expected.TargetResidue?.ToString()));
            });
        }
    }

    [Test]
    public void TheTp53Queries_GiveEveryOutcome_AndAGapIsNeverTheNearestResidue()
    {
        JsonElement data = Invoke(Stdin(Queries), Args());
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(columns.GetProperty("outcome")), Is.EqualTo(new[]
            {
                "Aligned", "Aligned", "GapInTarget", "Aligned", "Aligned", "Aligned",
                "DifferentAlignments", "SourceNotInAlignment", "TargetNotInAlignment",
            }));
            // Human S15, S20 and S392 are mouse S18, S23 and S389; S46 faces a gap in mouse.
            Assert.That(columns.GetProperty("target_position").EnumerateArray().Select(NullableInt),
                Is.EqualTo(new int?[] { 18, 23, null, 389, 15, 15, null, null, null }));
            Assert.That(columns.GetProperty("column").EnumerateArray().Select(NullableInt),
                Is.EqualTo(new int?[] { 476, 547, 863, 3036, 476, 476, null, null, null }));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(9));
            Assert.That(data.GetProperty("outcome_counts").EnumerateObject().Select(p => (p.Name, p.Value.GetInt32())), Is.EqualTo(new[]
            {
                ("Aligned", 5), ("GapInTarget", 1), ("DifferentAlignments", 1), ("SourceNotInAlignment", 1), ("TargetNotInAlignment", 1),
            }), "every outcome in mzLib's declaration order");
        });
    }

    [Test]
    public void Envelope_RecordsBothFiles_AndTheGeneSets_WithTheSpecsKeys()
    {
        JsonElement data = Invoke(Stdin(Queries), Args("--proteins"));
        ComparaGeneTreeAlignment mzlib = LoadWithMzLib();
        JsonElement alignment = data.GetProperty("alignment");
        JsonElement trees = data.GetProperty("gene_trees");

        Assert.Multiple(() =>
        {
            Assert.That(Keys(data), Is.EqualTo(new[]
            {
                "alignment", "gene_trees", "gene_sets", "row_count", "outcome_counts", "column_names", "columns",
                "proteins_included", "proteins", "caveats",
            }));
            Assert.That(Keys(alignment), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "collection", "alignment_count", "protein_count",
                "restricted_to_gene_tree_content",
            }));
            Assert.That(Keys(trees), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "collection", "member_count", "restricted_to_gene_sets",
            }), "genes orthology's gene_trees fields");
            Assert.That(Keys(data.GetProperty("gene_sets")), Is.EqualTo(new[] { Human, Mouse, Rat }));
            Assert.That(Keys(data.GetProperty("gene_sets").GetProperty(Human)), Is.EqualTo(new[]
            {
                "source_file_name", "sha256", "release", "genome_build", "genebuild_last_updated", "gene_count",
            }), "genes resolve's gene_set fields");
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(Proteins.ResidueMappingColumns));

            Assert.That(alignment.GetProperty("source_file_name").GetString(), Is.EqualTo(mzlib.SourceFileName));
            Assert.That(alignment.GetProperty("sha256").GetString(), Is.EqualTo(mzlib.SourceSha256));
            Assert.That(alignment.GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(alignment.GetProperty("collection").GetString(), Is.EqualTo("protein_default"));
            Assert.That(alignment.GetProperty("alignment_count").GetInt32(), Is.EqualTo(mzlib.AlignmentCount).And.EqualTo(2));
            Assert.That(alignment.GetProperty("protein_count").GetInt32(), Is.EqualTo(mzlib.Count).And.EqualTo(5));
            Assert.That(alignment.GetProperty("restricted_to_gene_tree_content").GetString(),
                Is.EqualTo(trees.GetProperty("sha256").GetString()));
            Assert.That(trees.GetProperty("member_count").GetInt32(), Is.EqualTo(5));
            Assert.That(trees.GetProperty("restricted_to_gene_sets").GetArrayLength(), Is.EqualTo(3));
            Assert.That(data.GetProperty("gene_sets").GetProperty(Human).GetProperty("genome_build").GetString(), Is.EqualTo("GRCh38.p14"));
            Assert.That(data.GetProperty("proteins_included").GetBoolean(), Is.True);
            Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).Contains("2 queries name a protein the alignment did not keep"));
        });
    }

    [Test]
    public void Proteins_AreMzLibsKeptProteins_InItsOrder()
    {
        JsonElement proteins = Invoke(Stdin(Queries), Args("--proteins")).GetProperty("proteins");
        List<ComparaAlignedProtein> expected = LoadWithMzLib().Proteins.ToList();
        JsonElement columns = proteins.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Keys(proteins), Is.EqualTo(new[] { "row_count", "column_names", "columns" }));
            Assert.That(Strings(proteins.GetProperty("column_names")), Is.EqualTo(Proteins.AlignedProteinColumns));
            Assert.That(proteins.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));
            Assert.That(Strings(columns.GetProperty("protein_id")), Is.EqualTo(expected.Select(p => p.ProteinId)));
            Assert.That(columns.GetProperty("alignment_index").EnumerateArray().Select(e => e.GetInt32()), Is.EqualTo(expected.Select(p => p.AlignmentIndex)));
            Assert.That(Strings(columns.GetProperty("tree_id")), Is.EqualTo(expected.Select(p => p.TreeId)));
            Assert.That(Strings(columns.GetProperty("sequence")), Is.EqualTo(expected.Select(p => p.Sequence)));
            Assert.That(columns.GetProperty("sequence")[Array.IndexOf(Strings(columns.GetProperty("protein_id")), Tp53)].GetString(), Has.Length.EqualTo(393));
        });
    }

    [Test]
    public void WithoutProteins_TheBlockIsNull_AndNoCaveatWhenEveryProteinWasKept()
    {
        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, Trp53) }), Args());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("proteins_included").GetBoolean(), Is.False);
            Assert.That(data.GetProperty("proteins").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("caveats").GetArrayLength(), Is.Zero);
        });
    }

    // ---- the gene sets only restrict -----------------------------------------------------------

    [Test]
    public void WithoutSpeciesLines_MzLibsDefaultKeepsEveryGene_AndACaveatSaysSo()
    {
        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, Trp53) }, species: Array.Empty<string>()), Args());

        Assert.Multiple(() =>
        {
            Assert.That(Keys(data.GetProperty("gene_sets")), Is.Empty);
            Assert.That(data.GetProperty("gene_trees").GetProperty("restricted_to_gene_sets").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("columns").GetProperty("target_position")[0].GetInt32(), Is.EqualTo(18));
            Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).StartsWith("No gene sets were given, so mzLib kept every gene in " + TreesName + " (5)"));
        });
    }

    [Test]
    public void ARatOnlyGeneSet_KeepsOnlyRatProteins_SoHumanIsNotInTheAlignment()
    {
        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, RatTp53) }, species: new[] { Rat }), Args());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("alignment").GetProperty("protein_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("columns").GetProperty("outcome")[0].GetString(), Is.EqualTo("SourceNotInAlignment"));
            Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).StartsWith("1 query names a protein"));
        });
    }

    [Test]
    public void AGeneSetTable_IsReadWithMzLibsReader()
    {
        string table = Path.Combine(_scratch, "rat.genes.tsv");
        EnsemblGeneSetWriter.Write(table, EnsemblGeneSet.LoadGtf(Fixture("Rattus_norvegicus.GRCr8.116.gtf")));
        string stdin = $"{Rat}\tgene-set\t{table}\n--\n{RatTp53}\t15\t{Tp53}\n";

        JsonElement data = Invoke(stdin, Args());
        Assert.That(data.GetProperty("gene_sets").GetProperty(Rat).GetProperty("source_file_name").GetString(),
            Is.EqualTo("Rattus_norvegicus.GRCr8.116.gtf"), "the table keeps the GTF's provenance");
    }

    // ---- release caveats -------------------------------------------------------------------------

    [Test]
    public void FilesNamedForDifferentReleases_GetACaveat_ButAnAnswer()
    {
        string trees = Copy(Fixture(TreesName), "vertebrates.GeneTree_content.default.e115.txt.gz");

        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, Trp53) }), ArgsWith(trees: trees));
        Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).StartsWith(
            "The alignment is named for release 116 and the gene trees for release 115."));
    }

    [Test]
    public void AnAlignmentNamedForNoRelease_GetsACaveat()
    {
        string alignment = Copy(Fixture(AlignmentName), "tp53_tree.aa.fasta.gz");

        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, Trp53) }), ArgsWith(alignment: alignment));
        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("alignment").GetProperty("release").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("alignment").GetProperty("collection").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).StartsWith("tp53_tree.aa.fasta.gz carries no release"));
        });
    }

    [Test]
    public void GeneTreesNamedForNoRelease_GetsACaveat()
    {
        string trees = Copy(Fixture(TreesName), "trees.txt.gz");

        JsonElement data = Invoke(Stdin(new[] { (Tp53, 15, Trp53) }), ArgsWith(trees: trees));
        Assert.That(Strings(data.GetProperty("caveats")), Has.Exactly(1).StartsWith("trees.txt.gz carries no release"));
    }

    // ---- mzLib's refusals ------------------------------------------------------------------------

    [TestCase(394)]
    [TestCase(0)]
    [TestCase(-1)]
    public void APositionOutsideTheSource_IsMzLibsRefusal_AndStopsTheCall(int position)
    {
        JsonElement error = InvokeExpectingError(Stdin(new[] { (Tp53, 15, Trp53), (Tp53, position, Trp53) }), Args());

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            string message = error.GetProperty("message").GetString()!;
            Assert.That(message, Does.StartWith($"Query 2 ({Tp53}, position {position}, {Trp53}): {Tp53} has 393 residues."));
            Assert.That(message, Does.Not.Contain("(Parameter"), "mzLib's message without .NET's suffix");
            Assert.That(message, Does.Contain("no query of this call was answered"));
        });
    }

    [Test]
    public void ARaggedAlignment_IsMzLibsInvalidDataException()
    {
        string bad = WriteGz(AlignmentName, ">ENSDARP00000000001\nPEPT\n>ENSDARP00000000002\nPEP\n//\n");

        JsonElement error = InvokeExpectingError(Stdin(new[] { (Tp53, 15, Trp53) }), ArgsWith(alignment: bad));
        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("InvalidDataException"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("is 3 columns wide, the alignment 4"));
        });
    }

    [Test]
    public void AMalformedGeneTreeRow_IsMzLibsInvalidDataException()
    {
        string bad = WriteGz(TreesName, "ENSGT1\tENSP00000269305\tENSG00000141510\n");

        JsonElement error = InvokeExpectingError(Stdin(new[] { (Tp53, 15, Trp53) }), ArgsWith(trees: bad));
        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("InvalidDataException"));
    }

    // ---- the bridge's own usage errors ----------------------------------------------------------

    [Test]
    public void AnOptionWithoutAValue_IsAUsageError()
    {
        Assert.Multiple(() =>
        {
            Assert.That(InvokeExpectingError(Stdin(Queries), new[] { "genes", "map-residues", "--gene-trees", Fixture(TreesName), "--alignment" })
                .GetProperty("message").GetString(), Is.EqualTo("Option --alignment was given without a value."));
            Assert.That(InvokeExpectingError(Stdin(Queries), new[] { "genes", "map-residues", "--alignment", Fixture(AlignmentName), "--gene-trees" })
                .GetProperty("message").GetString(), Is.EqualTo("Option --gene-trees was given without a value."));
        });
    }

    [Test]
    public void AMissingRequiredOption_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(Queries), new[] { "genes", "map-residues", "--gene-trees", Fixture(TreesName) });
        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo("Missing required option --alignment."));
    }

    [Test]
    public void AValueAfterProteins_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(Queries), Args("--proteins", "true"));
        Assert.That(error.GetProperty("message").GetString(), Does.StartWith("--proteins takes no value (got 'true')"));
    }

    [TestCase("ENSP00000269305\t15\n", "is not three tab-separated fields")]
    [TestCase("ENSP00000269305\t15\tENSMUSP00000104298\textra\n", "is not three tab-separated fields")]
    [TestCase("ENSP00000269305\t\tENSMUSP00000104298\n", "is not three tab-separated fields")]
    [TestCase("ENSP00000269305\tS15\tENSMUSP00000104298\n", "the position 'S15' is not an integer")]
    [TestCase("ENSP00000269305\t15.0\tENSMUSP00000104298\n", "the position '15.0' is not an integer")]
    [TestCase("", "No queries after the '--' line")]
    public void AMalformedQuery_IsAUsageError(string queries, string message)
    {
        string stdin = SpeciesLines(new[] { Human, Mouse, Rat }) + "--\n" + queries;
        JsonElement error = InvokeExpectingError(stdin, Args());
        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(message));
        });
    }

    [Test]
    public void NoSeparatorLine_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError($"{Tp53}\t15\t{Trp53}\n", Args());
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("No such line was found"));
    }

    [TestCase("homo_sapiens\tgtf\n", "is not three tab-separated fields")]
    [TestCase("homo_sapiens\tfasta\tx.fa\n", "names its gene set as 'fasta'")]
    public void AMalformedSpeciesLine_IsAUsageError(string speciesLine, string message)
    {
        JsonElement error = InvokeExpectingError(speciesLine + $"--\n{Tp53}\t15\t{Trp53}\n", Args());
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(message));
    }

    [Test]
    public void ASpeciesListedTwice_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError(Stdin(Queries, species: new[] { Human, Human }), Args());
        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo("Species 'homo_sapiens' is listed twice; give one gene set per species."));
    }

    [Test]
    public void MissingFiles_AreUsageErrors_BeforeAnyInputIsRead()
    {
        string absent = Path.Combine(_scratch, "absent.gz");
        Assert.Multiple(() =>
        {
            Assert.That(InvokeExpectingError(Stdin(Queries), ArgsWith(alignment: absent)).GetProperty("message").GetString(),
                Is.EqualTo($"Compara gene-tree alignment not found: '{absent}'."));
            Assert.That(InvokeExpectingError(Stdin(Queries), ArgsWith(trees: absent)).GetProperty("message").GetString(),
                Is.EqualTo($"Gene-trees file not found: '{absent}'."));
            Assert.That(InvokeExpectingError($"{Human}\tgtf\t{absent}\n--\n{Tp53}\t15\t{Trp53}\n", Args()).GetProperty("message").GetString(),
                Is.EqualTo($"GTF of {Human} not found: '{absent}'."));
            Assert.That(InvokeExpectingError($"{Human}\tgene-set\t{absent}\n--\n{Tp53}\t15\t{Trp53}\n", Args()).GetProperty("message").GetString(),
                Is.EqualTo($"Gene-set table of {Human} not found: '{absent}'."));
        });
    }

    // ---- plumbing -------------------------------------------------------------------------------

    private static readonly Dictionary<string, string> GtfOf = new()
    {
        [Human] = "Homo_sapiens.GRCh38.116.gtf",
        [Mouse] = "Mus_musculus.GRCm39.116.gtf",
        [Rat] = "Rattus_norvegicus.GRCr8.116.gtf",
    };

    private static string SpeciesLines(IEnumerable<string> species) =>
        string.Concat(species.Select(s => $"{s}\tgtf\t{Fixture(GtfOf[s])}\n"));

    private static string Stdin(IEnumerable<(string Source, int Position, string Target)> queries, string[]? species = null) =>
        SpeciesLines(species ?? new[] { Human, Mouse, Rat }) + "--\n" +
        string.Concat(queries.Select(q => $"{q.Source}\t{q.Position}\t{q.Target}\n"));

    private static string[] Args(params string[] extra) => ArgsWith(null, null).Concat(extra).ToArray();

    private static string[] ArgsWith(string? alignment = null, string? trees = null) =>
        new[] { "genes", "map-residues", "--alignment", alignment ?? Fixture(AlignmentName), "--gene-trees", trees ?? Fixture(TreesName) };

    /// <summary>mzLib's own reading of the same files, restricted to the same three gene sets.</summary>
    private static ComparaGeneTreeAlignment LoadWithMzLib()
    {
        var sets = GtfOf.Values.Select(g => EnsemblGeneSet.LoadGtf(Fixture(g))).ToList();
        return ComparaGeneTreeAlignment.Load(Fixture(AlignmentName), ComparaGeneTreeContent.Load(Fixture(TreesName), sets));
    }

    private string Copy(string source, string name)
    {
        string path = Path.Combine(_scratch, name);
        File.Copy(source, path);
        return path;
    }

    private string WriteGz(string name, string text)
    {
        string path = Path.Combine(_scratch, name);
        using (var file = File.Create(path))
        using (var gz = new GZipStream(file, CompressionLevel.Optimal))
        using (var writer = new StreamWriter(gz))
            writer.Write(text);
        return path;
    }

    /// <summary>A file of pkg/python/tests/fixtures/map_residues, committed beside these tests.</summary>
    private static string Fixture(string name, [CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "python", "tests", "fixtures", "map_residues", name));

    private static int? NullableInt(JsonElement e) => e.ValueKind == JsonValueKind.Null ? null : e.GetInt32();

    private static string? NullableString(JsonElement e) => e.ValueKind == JsonValueKind.Null ? null : e.GetString();

    private static string?[] Strings(JsonElement array) => array.EnumerateArray().Select(NullableString).ToArray();

    private static string[] Keys(JsonElement element) => element.EnumerateObject().Select(p => p.Name).ToArray();

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
