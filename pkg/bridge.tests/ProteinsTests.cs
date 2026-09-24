using System.Diagnostics.CodeAnalysis;
using System.IO.Compression;
using System.Runtime.CompilerServices;
using System.Text.Json;
using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>proteins read</c>, <c>genes resolve</c> and <c>proteins classify-peptides</c>, against
/// the small databases in <c>pkg/python/tests/fixtures/proteins/</c>.
/// </summary>
/// <remarks>
/// <para>
/// The fixtures are chosen so each outcome the verbs can report has a real example.
/// <b>human_subset.xml</b> is two real UniProt entries: GAPDH (P04406), whose Ensembl links carry
/// versioned gene ids, and PSCA (O43653), whose gene the tiny test GTF leaves out.
/// <b>human_extra.fasta</b> holds human albumin and three DYNC1I2 isoforms; <b>contaminants.fasta</b>
/// holds bovine albumin, which shares YLYEIAR with human albumin. <b>mouse_aifm1.fasta</b> is a mouse
/// entry, for the FASTA <c>OX=</c> taxonomy id.
/// </para>
/// <para>
/// The GTF and xref table are hand-made: three genes in Ensembl's own GTF layout, and four rows in
/// Ensembl's uniprot.tsv layout. mzLib at the pin ships no GTF test input; its own resolver tests build
/// theirs inline the same way.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ProteinsTests
{
    private string _temp = "";

    [SetUp]
    public void SetUp()
    {
        _temp = Path.Combine(Path.GetTempPath(), $"pymzlib-proteins-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_temp);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_temp))
            Directory.Delete(_temp, recursive: true);
    }

    // ---- proteins read ------------------------------------------------------------------------

    [Test]
    public void Read_UniProtXml_ReportsOrganismTaxonGenesAndMass()
    {
        JsonElement data = Invoke("proteins", "read", "--path", Fixture("human_subset.xml"));
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(columns, "accession"), Is.EqualTo(new[] { "P04406", "O43653" }));
            Assert.That(Strings(columns, "organism"), Is.All.EqualTo("Homo sapiens"));
            Assert.That(Strings(columns, "ncbi_taxonomy_id"), Is.All.EqualTo("9606"));
            Assert.That(Strings(columns, "primary_gene_name"), Is.EqualTo(new[] { "GAPDH", "PSCA" }));
            Assert.That(columns.GetProperty("length")[0].GetInt32(), Is.EqualTo(335));
            Assert.That(columns.GetProperty("monoisotopic_mass")[0].GetDouble(), Is.EqualTo(36030.3979).Within(0.001));
            Assert.That(columns.GetProperty("ensembl_gene_ids")[0][0].GetString(), Is.EqualTo("ENSG00000111640"));
            Assert.That(data.GetProperty("go_terms").ValueKind, Is.EqualTo(JsonValueKind.Null), "not asked for");
            Assert.That(data.GetProperty("tables").EnumerateArray().Select(t => t.GetString()), Is.EqualTo(new[] { "proteins" }));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(2));
        });
    }

    [Test]
    public void Read_Fasta_TakesTheTaxonFromOx_AndSaysGoAndEnsemblAreAbsent()
    {
        JsonElement data = Invoke("proteins", "read", "--path", Fixture("mouse_aifm1.fasta"),
            "--tables", "proteins,go_terms,ensembl_genes");
        JsonElement file = data.GetProperty("files")[0];

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data.GetProperty("columns"), "ncbi_taxonomy_id"), Is.EqualTo(new[] { "10090" }));
            Assert.That(Strings(data.GetProperty("columns"), "organism"), Is.EqualTo(new[] { "Mus musculus" }));
            Assert.That(file.GetProperty("file_type").GetString(), Is.EqualTo("Fasta"));
            Assert.That(file.GetProperty("absent_fields").EnumerateArray().Select(a => a.GetString()),
                Is.SupersetOf(new[] { "go_terms", "ensembl_genes" }));
            Assert.That(data.GetProperty("go_terms").GetProperty("row_count").GetInt32(), Is.Zero);
        });
    }

    [Test]
    public void Read_FastaWithoutOx_HasANullTaxon_NotAGuessedOne()
    {
        string fasta = Write("no_ox.fasta", ">sp|P62805|H4_HUMAN Histone H4 OS=Homo sapiens GN=H4C1 PE=1 SV=2\nMSGRGKGGKGLGKGGAKRHRK\n");

        JsonElement columns = Invoke("proteins", "read", "--path", fasta).GetProperty("columns");

        Assert.That(columns.GetProperty("ncbi_taxonomy_id")[0].ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void Read_GoTerms_AreOneRowPerTermWithAspectAndSortedEvidence()
    {
        JsonElement go = Invoke("proteins", "read", "--path", Fixture("human_subset.xml"), "--tables", "go_terms")
            .GetProperty("go_terms");
        JsonElement columns = go.GetProperty("columns");
        string[] aspects = Strings(columns, "aspect");
        int cytoplasm = Array.IndexOf(Strings(columns, "go_id"), "GO:0005737");

        Assert.Multiple(() =>
        {
            Assert.That(go.GetProperty("row_count").GetInt32(), Is.EqualTo(aspects.Length).And.GreaterThan(0));
            Assert.That(aspects, Is.All.AnyOf("BiologicalProcess", "CellularComponent", "MolecularFunction", "Unknown"));
            Assert.That(Strings(columns, "term_name")[cytoplasm], Is.EqualTo("cytoplasm"), "the C: prefix is removed");
            Assert.That(aspects[cytoplasm], Is.EqualTo("CellularComponent"));
            foreach (JsonElement codes in columns.GetProperty("evidence_codes").EnumerateArray())
            {
                string[] values = codes.EnumerateArray().Select(c => c.GetString()!).ToArray();
                Assert.That(values, Is.Ordered.Using((IComparer<string>)StringComparer.Ordinal));
            }
        });
    }

    [Test]
    public void Read_EnsemblGenes_KeepTheVersionedIdBesideTheStableOne()
    {
        JsonElement ensembl = Invoke("proteins", "read", "--path", Fixture("human_subset.xml"), "--tables", "ensembl_genes")
            .GetProperty("ensembl_genes").GetProperty("columns");

        int gapdh = Array.IndexOf(Strings(ensembl, "versioned_gene_id"), "ENSG00000111640.15");
        int psca = Array.IndexOf(Strings(ensembl, "gene_id"), "ENSG00000167653");

        Assert.Multiple(() =>
        {
            Assert.That(gapdh, Is.GreaterThanOrEqualTo(0));
            Assert.That(ensembl.GetProperty("gene_version")[gapdh].GetInt32(), Is.EqualTo(15));
            Assert.That(ensembl.GetProperty("gene_version")[psca].ValueKind, Is.EqualTo(JsonValueKind.Null),
                "an unversioned id has no version; absent is not 0");
        });
    }

    [Test]
    public void Read_AccessionFilter_KeepsOnlyThoseAndNamesTheMisses()
    {
        JsonElement data = InvokeWithStdin("O43653\nP99999\n", "proteins", "read", "--path", Fixture("human_subset.xml"),
            "--accessions-stdin", "--tables", "proteins,go_terms");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(data.GetProperty("columns"), "accession"), Is.EqualTo(new[] { "O43653" }));
            Assert.That(Strings(data.GetProperty("go_terms").GetProperty("columns"), "accession"), Is.All.EqualTo("O43653"));
            Assert.That(data.GetProperty("accession_filter_count").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("accessions_not_found").EnumerateArray().Select(a => a.GetString()),
                Is.EqualTo(new[] { "P99999" }));
        });
    }

    [Test]
    public void Read_ManyDatabases_AreInInputOrder_AndIdenticalAtAnyThreadCount()
    {
        string stdin = string.Join('\n', Fixture("human_subset.xml"), Fixture("human_extra.fasta"),
            Fixture("contaminants.fasta") + "\tcontaminant", Fixture("mouse_aifm1.fasta"));

        string one = InvokeWithStdin(stdin, "proteins", "read", "--paths-stdin", "--threads", "1",
            "--tables", "proteins,go_terms,ensembl_genes").GetRawText();
        JsonElement four = InvokeWithStdin(stdin, "proteins", "read", "--paths-stdin", "--threads", "4",
            "--tables", "proteins,go_terms,ensembl_genes");

        Assert.Multiple(() =>
        {
            Assert.That(four.GetRawText(), Is.EqualTo(one), "output must not depend on --threads");
            Assert.That(four.GetProperty("columns").GetProperty("source_index").EnumerateArray().Select(i => i.GetInt32()),
                Is.Ordered);
            Assert.That(four.GetProperty("files")[2].GetProperty("contaminant").GetBoolean(), Is.True);
            Assert.That(Bools(four.GetProperty("columns"), "is_contaminant").Count(b => b), Is.EqualTo(1));
            Assert.That(four.GetProperty("file_count").GetInt32(), Is.EqualTo(4));
        });
    }

    [Test]
    public void Read_AllCoresIsAccepted()
    {
        JsonElement data = Invoke("proteins", "read", "--path", Fixture("mouse_aifm1.fasta"), "--threads", "-1");

        Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(1));
    }

    [Test]
    public void Read_PathsAndAccessionsShareStdin_SeparatedByADashDashLine()
    {
        string stdin = $"{Fixture("human_subset.xml")}\n{Fixture("human_extra.fasta")}\n--\nP02768\nP04406\n";

        JsonElement data = InvokeWithStdin(stdin, "proteins", "read", "--paths-stdin", "--accessions-stdin");

        Assert.That(Strings(data.GetProperty("columns"), "accession"), Is.EqualTo(new[] { "P04406", "P02768" }));
    }

    [Test]
    public void Read_WithSequences_AddsTheSequenceColumn()
    {
        JsonElement data = Invoke("proteins", "read", "--path", Fixture("mouse_aifm1.fasta"), "--sequences");

        Assert.That(Strings(data.GetProperty("columns"), "sequence")[0], Does.StartWith("MFRCGGLAGAF"));
    }

    [Test]
    public void Read_ASequenceWithX_HasANullMass_AndACaveatSaysWhy()
    {
        string fasta = Write("x.fasta", ">sp|Q00001|TEST_HUMAN Test OS=Homo sapiens OX=9606 GN=TST PE=1 SV=1\nPEPXTIDEK\n");

        JsonElement data = Invoke("proteins", "read", "--path", fasta);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").GetProperty("monoisotopic_mass")[0].ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("caveats")[0].GetString(), Does.Contain("no monoisotopic_mass"));
        });
    }

    [Test]
    public void Read_AGzippedDatabase_Reads()
    {
        string gz = Path.Combine(_temp, "mouse.fasta.gz");
        using (FileStream output = File.Create(gz))
        using (var zip = new GZipStream(output, CompressionMode.Compress))
            File.OpenRead(Fixture("mouse_aifm1.fasta")).CopyTo(zip);

        JsonElement data = Invoke("proteins", "read", "--path", gz);

        Assert.That(Strings(data.GetProperty("columns"), "accession"), Is.EqualTo(new[] { "Q9Z0X1" }));
    }

    [Test]
    public void Read_AnAppliedVariant_IsDisclosedAsACaveat()
    {
        string seqVar = MzLibDatabase("SeqVar.xml");

        JsonElement file = Invoke("proteins", "read", "--path", seqVar).GetProperty("files")[0];

        Assert.That(file.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()),
            Has.Some.Contains("sequence variants mzLib applied"));
    }

    [Test]
    public void Read_OnErrorSkip_RecordsTheFailureAndReadsTheRest()
    {
        string missing = Path.Combine(_temp, "absent.fasta");
        string stdin = $"{Fixture("mouse_aifm1.fasta")}\n{missing}\n";

        JsonElement data = InvokeWithStdin(stdin, "proteins", "read", "--paths-stdin", "--on-error", "skip");
        JsonElement failed = data.GetProperty("files")[1];

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("read_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("failed_count").GetInt32(), Is.EqualTo(1));
            Assert.That(failed.GetProperty("error").GetProperty("kind").GetString(), Is.EqualTo("usage"));
            Assert.That(failed.GetProperty("error").GetProperty("message").GetString(), Does.Contain("absent.fasta"));
        });
    }

    [Test]
    public void Read_OnErrorFail_ReportsTheFirstFailingFileInInputOrder()
    {
        string first = Path.Combine(_temp, "first-missing.fasta");
        string second = Path.Combine(_temp, "second-missing.fasta");

        JsonElement error = InvokeWithStdinExpectingError($"{first}\n{second}\n", "proteins", "read", "--paths-stdin",
            "--threads", "4");

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("first-missing"));
    }

    [Test]
    public void Read_AFileThatCannotParse_UnderSkip_IsACorrectnessFailure()
    {
        string broken = Write("broken.xml", "<uniprot><entry><accession>P1</accession><sequence>MPEP");

        JsonElement data = InvokeWithStdin(broken + "\n", "proteins", "read", "--paths-stdin", "--on-error", "skip");

        Assert.That(data.GetProperty("files")[0].GetProperty("error").GetProperty("kind").GetString(), Is.EqualTo("correctness"));
    }

    [TestCase(new[] { "proteins", "read" }, "", "Missing required option --path")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--paths-stdin" }, "", "not both")]
    [TestCase(new[] { "proteins", "read", "--paths-stdin", "--contaminant" }, "a.fasta", "--contaminant applies to --path")]
    [TestCase(new[] { "proteins", "read", "--paths-stdin" }, "", "no database paths")]
    [TestCase(new[] { "proteins", "read", "--paths-stdin" }, "a.fasta\na.fasta", "listed twice")]
    [TestCase(new[] { "proteins", "read", "--paths-stdin" }, "a.fasta\tdecoy", "only roles")]
    [TestCase(new[] { "proteins", "read", "--paths-stdin", "--accessions-stdin" }, "a.fasta\nP04406", "only '--'")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--accessions-stdin" }, "", "no accessions")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--tables", "proteins,peptides" }, "", "--tables")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--threads", "0" }, "", "--threads")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--threads", "-2" }, "", "--threads")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--on-error", "ignore" }, "", "'fail' or 'skip'")]
    [TestCase(new[] { "proteins", "read", "--path", "a.fasta", "--threads" }, "", "without a value")]
    [TestCase(new[] { "proteins", "read", "--path", "proteins.csv" }, "", "not a protein database")]
    [TestCase(new[] { "proteins", "read", "--path", "absent.fasta" }, "", "not found")]
    public void Read_BadInput_IsAUsageError(string[] args, string stdin, string expected)
    {
        JsonElement error = InvokeWithStdinExpectingError(stdin, args);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    // ---- genes resolve ------------------------------------------------------------------------

    [Test]
    public void Resolve_ColumnsAreMzLibsOwnTsvSchema()
    {
        JsonElement data = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"));
        string[] names = data.GetProperty("column_names").EnumerateArray().Select(n => n.GetString()!).ToArray();

        Assert.That(names.Skip(2), Is.EqualTo(GeneResolutionTsv.Schema.Select(c => c.Header)));
    }

    [Test]
    public void Resolve_ReportsEachOutcomeAndThePinningHashes()
    {
        string stdin = string.Join('\n', Fixture("human_subset.xml"), Fixture("human_extra.fasta"),
            Fixture("contaminants.fasta") + "\tcontaminant");

        JsonElement data = InvokeWithStdin(stdin, "genes", "resolve", "--paths-stdin",
            "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"), "--xref", Fixture("Homo_sapiens.GRCh38.116.uniprot.tsv"));
        JsonElement columns = data.GetProperty("columns");
        string[] accessions = Strings(columns, "accession");
        string[] outcomes = Strings(columns, "outcome");
        int gapdh = Array.IndexOf(accessions, "P04406");

        Assert.Multiple(() =>
        {
            Assert.That(outcomes[gapdh], Is.EqualTo("resolved"));
            Assert.That(Strings(columns, "gene_symbol")[gapdh], Is.EqualTo("GAPDH"));
            Assert.That(Strings(columns, "versioned_gene_id")[gapdh], Is.EqualTo("ENSG00000111640.15"));
            Assert.That(columns.GetProperty("ensembl_xref_agrees")[gapdh].GetBoolean(), Is.True);
            Assert.That(outcomes[Array.IndexOf(accessions, "O43653")], Is.EqualTo("off_primary_only"));
            Assert.That(outcomes[Array.IndexOf(accessions, "P02768")], Is.EqualTo("not_in_source"));
            Assert.That(outcomes[Array.IndexOf(accessions, "P02769")], Is.EqualTo("contaminant_not_mapped"));
            Assert.That(columns.GetProperty("isoform")[Array.IndexOf(accessions, "Q13409-2")].GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("gene_set").GetProperty("release").GetString(), Is.EqualTo("116"));
            Assert.That(data.GetProperty("gene_set").GetProperty("genome_build").GetString(), Is.EqualTo("GRCh38.p14"));
            Assert.That(data.GetProperty("outcome_counts").GetProperty("not_in_source").GetInt32(), Is.EqualTo(4));
            Assert.That(data.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()),
                Has.Some.Contains("came from a FASTA"));
            Assert.That(Strings(columns, "search_database_sha256").Distinct().Count(), Is.EqualTo(3), "one hash per database");
        });
    }

    [Test]
    public void Resolve_AGeneOnlyEnsemblLinks_IsItsOwnRow_AndChangesNoOutcome()
    {
        string xref = Write("Homo_sapiens.GRCh38.116.uniprot.tsv",
            string.Join('\t', EnsemblXrefTable.Columns) + "\n" +
            "ENSG00000000003\tENST1\tENSP1\tP04406\tUniprot/SWISSPROT\tINFERRED_PAIR\t100\t100\t-\n");

        JsonElement columns = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"),
            "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"), "--xref", xref).GetProperty("columns");
        string[] sources = Strings(columns, "source");
        int extra = Array.IndexOf(sources, "ensembl_xref");

        Assert.Multiple(() =>
        {
            Assert.That(extra, Is.GreaterThanOrEqualTo(0));
            Assert.That(Strings(columns, "gene_symbol")[extra], Is.EqualTo("TSPAN6"));
            Assert.That(Strings(columns, "outcome")[extra], Is.EqualTo("resolved"), "the search database's outcome");
            Assert.That(Strings(columns, "ensembl_xref_info_type")[extra], Is.EqualTo("INFERRED_PAIR"));
        });
    }

    [Test]
    public void Resolve_WithoutXref_AgreementIsNull_AndACaveatSaysItIsUnknownNotFalse()
    {
        JsonElement data = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"));

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").GetProperty("ensembl_xref_agrees").EnumerateArray()
                .Select(v => v.ValueKind), Is.All.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("xref").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()), Has.Some.Contains("unknown, not false"));
        });
    }

    [Test]
    public void Resolve_AGtfWithoutReleaseOrBuild_SaysSo()
    {
        string gtf = Write("genes.gtf", "1\thavana\tgene\t1\t9\t.\t+\t.\tgene_id \"ENSG00000111640\"; gene_biotype \"protein_coding\";\n");

        JsonElement data = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gtf", gtf);
        string[] caveats = data.GetProperty("caveats").EnumerateArray().Select(c => c.GetString()!).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("gene_set").GetProperty("release").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(caveats, Has.Some.Contains("no Ensembl release number"));
            Assert.That(caveats, Has.Some.Contains("genome-build"));
        });
    }

    [Test]
    public void Resolve_TheDatabaseHash_IsOfTheDecompressedBytes()
    {
        string gz = Path.Combine(_temp, "human_subset.xml.gz");
        using (FileStream output = File.Create(gz))
        using (var zip = new GZipStream(output, CompressionMode.Compress))
            File.OpenRead(Fixture("human_subset.xml")).CopyTo(zip);

        string plain = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"))
            .GetProperty("files")[0].GetProperty("search_database_sha256").GetString()!;
        string zipped = Invoke("genes", "resolve", "--path", gz, "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"))
            .GetProperty("files")[0].GetProperty("search_database_sha256").GetString()!;

        Assert.That(zipped, Is.EqualTo(plain).And.Length.EqualTo(64));
    }

    [Test]
    public void Resolve_ACompactGeneSet_KeysRowsExactlyAsItsGtf()
    {
        EnsemblGeneSet set = EnsemblGeneSet.LoadGtf(Fixture("Homo_sapiens.GRCh38.116.gtf"));
        string table = Path.Combine(_temp, "genes.tsv");
        EnsemblGeneSetWriter.Write(table, set);

        JsonElement fromGtf = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"));
        JsonElement fromTable = Invoke("genes", "resolve", "--path", Fixture("human_subset.xml"), "--gene-set", table);

        Assert.That(fromTable.GetProperty("columns").GetRawText(), Is.EqualTo(fromGtf.GetProperty("columns").GetRawText()));
    }

    [TestCase(new[] { "genes", "resolve", "--path", "a.xml" }, "exactly one of --gtf")]
    [TestCase(new[] { "genes", "resolve", "--path", "a.xml", "--gtf", "a.gtf", "--gene-set", "b.tsv" }, "exactly one of --gtf")]
    [TestCase(new[] { "genes", "resolve", "--path", "a.xml", "--gtf", "absent.gtf" }, "GTF not found")]
    [TestCase(new[] { "genes", "resolve", "--path", "a.xml", "--gene-set", "absent.tsv" }, "Gene set not found")]
    [TestCase(new[] { "genes", "resolve", "--path", "a.xml", "--gtf" }, "without a value")]
    public void Resolve_BadInput_IsAUsageError(string[] args, string expected)
    {
        JsonElement error = InvokeExpectingError(args);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    [Test]
    public void Resolve_AMissingXref_IsAUsageError()
    {
        JsonElement error = InvokeExpectingError("genes", "resolve", "--path", Fixture("human_subset.xml"),
            "--gtf", Fixture("Homo_sapiens.GRCh38.116.gtf"), "--xref", Path.Combine(_temp, "absent.tsv"));

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("xref table not found"));
    }

    // ---- proteins classify-peptides -----------------------------------------------------------

    private string ClassificationStdin(params string[] peptides) =>
        string.Join('\n', Fixture("human_subset.xml"), Fixture("human_extra.fasta"),
            Fixture("contaminants.fasta") + "\tcontaminant", "--") + "\n" + string.Join('\n', peptides);

    [Test]
    public void Classify_ReportsEachSharingClass()
    {
        JsonElement data = InvokeWithStdin(
            ClassificationStdin("VGVNGFGR", "ALSEQINIFFDYSGR", "YLYEIAR", "PEPTIDEK"),
            "proteins", "classify-peptides", "--paths-stdin");
        JsonElement columns = data.GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(columns, "sharing"),
                Is.EqualTo(new[] { "Unique", "SharedWithinGene", "SharedAcrossGenes", "NotInDatabase" }));
            Assert.That(columns.GetProperty("accessions")[2].EnumerateArray().Select(a => a.GetString()),
                Is.EqualTo(new[] { "P02768", "P02769" }), "human and bovine albumin");
            Assert.That(columns.GetProperty("shared_gene_keys")[1].EnumerateArray().Select(a => a.GetString()),
                Has.Member("gene:Homo sapiens:DYNC1I2"));
            Assert.That(data.GetProperty("i_and_l_equivalent").GetBoolean(), Is.True);
            Assert.That(data.GetProperty("target_protein_count").GetInt32(), Is.EqualTo(7));
            Assert.That(data.GetProperty("sharing_counts").GetProperty("Unique").GetInt32(), Is.EqualTo(1));
        });
    }

    [Test]
    public void Classify_TreatsIAndLAsOneResidue_AndEchoesThePeptideAsGiven()
    {
        // GAPDH has LVINGNPITIFQER. Every I written as L must still find it.
        JsonElement columns = InvokeWithStdin("LVLNGNPLTLFQER\n", "proteins", "classify-peptides", "--path", Fixture("human_subset.xml")).GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(Strings(columns, "peptide"), Is.EqualTo(new[] { "LVLNGNPLTLFQER" }));
            Assert.That(Strings(columns, "sharing"), Is.EqualTo(new[] { "Unique" }));
        });
    }

    [Test]
    public void Classify_IgnoresDecoysAndCountsThem()
    {
        string fasta = Write("with_decoy.fasta",
            ">sp|P00001|A_HUMAN A OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\nMPEPTIDEKAAA\n" +
            ">DECOY_P00002 decoy\nMPEPTIDEKCCC\n");

        JsonElement data = InvokeWithStdin("PEPTIDEK\n", "proteins", "classify-peptides", "--path", fasta);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("decoy_proteins_ignored").GetInt32(), Is.EqualTo(1));
            Assert.That(Strings(data.GetProperty("columns"), "sharing"), Is.EqualTo(new[] { "Unique" }));
        });
    }

    [TestCase(new[] { "proteins", "classify-peptides", "--path", "a.fasta" }, "", "No peptides")]
    [TestCase(new[] { "proteins", "classify-peptides", "--paths-stdin" }, "a.fasta\n--\n", "No peptides")]
    [TestCase(new[] { "proteins", "classify-peptides", "--path", "a.fasta", "--on-error", "skip" }, "PEPTIDEK", "no --on-error skip")]
    public void Classify_BadInput_IsAUsageError(string[] args, string stdin, string expected)
    {
        JsonElement error = InvokeWithStdinExpectingError(stdin, args);

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
        });
    }

    [TestCase("peptidek")]
    [TestCase("PEPT[Phospho]IDEK")]
    public void Classify_AModifiedOrLowerCasePeptide_IsAUsageErrorNamingIt(string peptide)
    {
        JsonElement error = InvokeWithStdinExpectingError(peptide + "\n", "proteins", "classify-peptides",
            "--path", Fixture("human_subset.xml"));

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain(peptide).And.Contain("base sequence"));
        });
    }

    [Test]
    public void Classify_ExplicitFailIsAccepted()
    {
        JsonElement data = InvokeWithStdin("VGVNGFGR\n", "proteins", "classify-peptides", "--path", Fixture("human_subset.xml"),
            "--on-error", "fail");

        Assert.That(data.GetProperty("peptide_count").GetInt32(), Is.EqualTo(1));
    }

    // ---- helpers ------------------------------------------------------------------------------

    /// <summary>A database from pkg/python/tests/fixtures/proteins, which is committed beside these tests.</summary>
    private static string Fixture(string name, [CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "python", "tests", "fixtures", "proteins", name));

    /// <summary>A database from the pinned mzLib worktree, or an ignored test without one.</summary>
    private static string MzLibDatabase(string name, [CallerFilePath] string thisFile = "")
    {
        string root = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));
        string path = Path.Combine(root, "code", "mzLib", "mzLib", "Test", "DatabaseTests", name);
        if (!File.Exists(path))
            Assert.Ignore($"mzLib database fixture not present in the worktree: {path}");
        return path;
    }

    private string Write(string name, string content)
    {
        string path = Path.Combine(_temp, name);
        File.WriteAllText(path, content);
        return path;
    }

    private static string[] Strings(JsonElement columns, string name) =>
        columns.GetProperty(name).EnumerateArray().Select(x => x.ValueKind == JsonValueKind.Null ? null! : x.GetString()!).ToArray();

    private static bool[] Bools(JsonElement columns, string name) =>
        columns.GetProperty(name).EnumerateArray().Select(x => x.GetBoolean()).ToArray();

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(null, args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) => Unwrap(Envelope(null, args), expectOk: false);

    private static JsonElement InvokeWithStdin(string stdin, params string[] args) => Unwrap(Envelope(stdin, args), expectOk: true);

    private static JsonElement InvokeWithStdinExpectingError(string stdin, params string[] args) =>
        Unwrap(Envelope(stdin, args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk), $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    /// <summary>Dispatches a verb with <paramref name="stdin"/> as the console, shaped as <c>Program.Main</c> would.</summary>
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
            return JsonSerializer.SerializeToElement(new { ok = false, error = new { type = "usage", message = usage.Message } },
                Program.JsonOptions);
        }
        catch (Exception exception)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = Program.ClassifyError(exception), message = Program.Unwrap(exception).Message } },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetIn(previousIn);
        }
    }
}
