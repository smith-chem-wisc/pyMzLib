using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Omics.SequenceAlignment;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>proteins align</c>, the projection of mzLib's <c>PairwiseAligner</c> (#1428).
/// </summary>
/// <remarks>
/// The short pairs are mzLib's own test cases (Test/Omics/SequenceAlignment/TestPairwiseAligner.cs),
/// so the scores here are mzLib's documented answers. The real pairs are the DYNC1I2 isoforms and the
/// human and bovine albumins in pyMzLib's committed FASTA fixtures, which the proteins guide aligns.
/// What is pinned: an option left out is mzLib's default and is echoed; each residue maps to its
/// counterpart or to null, never to a neighbour; mzLib's refusals are usage errors naming the pair.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ProteinsAlignTests
{
    // mzLib's Align_AnInternalDeletion_LeavesTheDeletedResiduesFacingAGap: QRQ (source 9-11) is deleted.
    private const string Deletion = "MKTAYIAKQRQISFVKSHFSRQ\tMKTAYIAKISFVKSHFSRQ\n";

    // mzLib's Align_AnExtraNTerminus_IsFreeByDefault_AndChargedWhenAsked.
    private const string Overhang = "MLLLPEPTIDE\tPEPTIDE\n";

    [Test]
    public void ADeletion_LeavesItsResiduesUnmappedAndMapsAcrossTheGap()
    {
        JsonElement data = Invoke(Deletion, "proteins", "align");
        JsonElement columns = data.GetProperty("columns");
        JsonElement map = data.GetProperty("residue_map").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(1));
            Assert.That(columns.GetProperty("score")[0].GetInt32(), Is.EqualTo(39 + 55 - 14)); // mzLib's own figure
            Assert.That(columns.GetProperty("identical")[0].GetInt32(), Is.EqualTo(19));
            Assert.That(columns.GetProperty("substituted")[0].GetInt32(), Is.Zero);
            Assert.That(columns.GetProperty("source_residues_facing_gap")[0].GetInt32(), Is.EqualTo(3));
            Assert.That(columns.GetProperty("target_residues_facing_gap")[0].GetInt32(), Is.Zero);
            Assert.That(columns.GetProperty("source_length")[0].GetInt32(), Is.EqualTo(22));
            Assert.That(columns.GetProperty("target_length")[0].GetInt32(), Is.EqualTo(19));
            Assert.That(columns.GetProperty("aligned_source")[0].GetString(), Is.EqualTo("MKTAYIAKQRQISFVKSHFSRQ"));
            Assert.That(columns.GetProperty("aligned_target")[0].GetString(), Is.EqualTo("MKTAYIAK---ISFVKSHFSRQ"));

            // One row per residue of each sequence: 22 source rows, then 19 target rows.
            Assert.That(data.GetProperty("residue_map").GetProperty("row_count").GetInt32(), Is.EqualTo(41));
            Assert.That(Strings(map.GetProperty("sequence")),
                Is.EqualTo(Enumerable.Repeat("source", 22).Concat(Enumerable.Repeat("target", 19))));
            Assert.That(Ints(map.GetProperty("position")),
                Is.EqualTo(Enumerable.Range(1, 22).Concat(Enumerable.Range(1, 19)).Cast<int?>()));

            // Source 9-11 (Q, R, Q) face the gap: null, never the nearest residue.
            Assert.That(Ints(map.GetProperty("mapped_position")).Skip(8).Take(3), Is.All.Null);
            Assert.That(Strings(map.GetProperty("mapped_residue")).Skip(8).Take(3), Is.All.Null);
            Assert.That(Strings(map.GetProperty("residue")).Skip(8).Take(3), Is.EqualTo(new[] { "Q", "R", "Q" }));
            // Source 12 (I) is target 9, and target 9 is source 12: both directions are mzLib's TryMap.
            Assert.That(map.GetProperty("mapped_position")[11].GetInt32(), Is.EqualTo(9));
            Assert.That(map.GetProperty("mapped_residue")[11].GetString(), Is.EqualTo("I"));
            Assert.That(map.GetProperty("mapped_position")[22 + 8].GetInt32(), Is.EqualTo(12));
            Assert.That(Ints(map.GetProperty("pair_index")), Is.All.EqualTo(0));
        });
    }

    [Test]
    public void OptionsLeftOut_AreMzLibsDefaults_AndTheEnvelopeSaysWhichWereUsed()
    {
        JsonElement data = Invoke(Overhang, "proteins", "align");
        var mzLib = new PairwiseAligner();

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("aligner_id").GetString(), Is.EqualTo(mzLib.Id));
            Assert.That(data.GetProperty("aligner_id").GetString(),
                Is.EqualTo("global-affine;BLOSUM62;open=11;extend=1;end-gaps=free"));
            Assert.That(data.GetProperty("matrix").GetString(), Is.EqualTo(mzLib.Matrix.Name));
            Assert.That(data.GetProperty("gap_open").GetInt32(), Is.EqualTo(mzLib.GapOpen));
            Assert.That(data.GetProperty("gap_extend").GetInt32(), Is.EqualTo(mzLib.GapExtend));
            Assert.That(data.GetProperty("end_gaps").GetString(), Is.EqualTo("free"));
            Assert.That(data.GetProperty("max_cells").GetInt64(), Is.EqualTo(mzLib.MaxCells));
            Assert.That(data.GetProperty("max_cells").GetInt64(), Is.EqualTo(200_000_000));
            // Free end gaps: the extra MLLL costs nothing (mzLib's score, 39).
            Assert.That(data.GetProperty("columns").GetProperty("score")[0].GetInt32(), Is.EqualTo(39));
            Assert.That(data.GetProperty("columns").GetProperty("aligned_target")[0].GetString(), Is.EqualTo("----PEPTIDE"));
            Assert.That(data.GetProperty("caveats").GetArrayLength(), Is.EqualTo(7));
        });
    }

    [Test]
    public void EndGapsCharged_ChargesTheOverhang()
    {
        JsonElement data = Invoke(Overhang, "proteins", "align", "--end-gaps", "charged");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("end_gaps").GetString(), Is.EqualTo("charged"));
            Assert.That(data.GetProperty("aligner_id").GetString(), Does.EndWith("end-gaps=charged"));
            // mzLib's figure: 39 less one 4-residue gap, 11 + 4 x 1.
            Assert.That(data.GetProperty("columns").GetProperty("score")[0].GetInt32(), Is.EqualTo(39 - (11 + 4)));
        });
    }

    [Test]
    public void GapCosts_ReachTheAligner()
    {
        JsonElement data = Invoke(Overhang, "proteins", "align", "--end-gaps", "CHARGED", "--gap-open", "5",
            "--gap-extend", "2", "--max-cells", "1000");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("aligner_id").GetString(),
                Is.EqualTo("global-affine;BLOSUM62;open=5;extend=2;end-gaps=charged"));
            Assert.That(data.GetProperty("gap_open").GetInt32(), Is.EqualTo(5));
            Assert.That(data.GetProperty("gap_extend").GetInt32(), Is.EqualTo(2));
            Assert.That(data.GetProperty("max_cells").GetInt64(), Is.EqualTo(1000));
            Assert.That(data.GetProperty("columns").GetProperty("score")[0].GetInt32(), Is.EqualTo(39 - (5 + 4 * 2)));
        });
    }

    [Test]
    public void ManyPairs_AreRowsInInputOrder_BlankLinesSkipped()
    {
        JsonElement data = Invoke(Deletion + "\n  \n" + Overhang + Deletion, "proteins", "align");
        JsonElement columns = data.GetProperty("columns");
        JsonElement map = data.GetProperty("residue_map");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(3));
            Assert.That(Ints(columns.GetProperty("pair_index")), Is.EqualTo(new int?[] { 0, 1, 2 }));
            Assert.That(Ints(columns.GetProperty("score")), Is.EqualTo(new int?[] { 80, 39, 80 }));
            Assert.That(map.GetProperty("row_count").GetInt32(), Is.EqualTo(41 + 18 + 41));
            Assert.That(Ints(map.GetProperty("columns").GetProperty("pair_index")),
                Is.EqualTo(Enumerable.Repeat(0, 41).Concat(Enumerable.Repeat(1, 18)).Concat(Enumerable.Repeat(2, 41)).Cast<int?>()));
        });
    }

    [Test]
    public void NoResidueMap_DropsOnlyTheMap()
    {
        JsonElement with = Invoke(Deletion, "proteins", "align");
        JsonElement without = Invoke(Deletion, "proteins", "align", "--no-residue-map");

        Assert.Multiple(() =>
        {
            Assert.That(without.GetProperty("residue_map").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(without.GetProperty("columns").GetRawText(), Is.EqualTo(with.GetProperty("columns").GetRawText()));
        });
    }

    [Test]
    public void IdenticalSequences_AlignWithoutTheTable_SoMaxCellsDoesNotLimitThem()
    {
        // mzLib skips the dynamic programme for two identical sequences free of X and *, so a pair
        // that would need 64 cells aligns under --max-cells 1, and one residue different does not.
        JsonElement same = Invoke("PEPTIDE\tpeptide\n", "proteins", "align", "--max-cells", "1");
        JsonElement error = InvokeExpectingError("PEPTIDE\tPEPTIDE\nAAAA\tAAAC\n", "proteins", "align", "--max-cells", "1");

        Assert.Multiple(() =>
        {
            Assert.That(same.GetProperty("columns").GetProperty("identical")[0].GetInt32(), Is.EqualTo(7));
            // Case is compared insensitively and kept as given.
            Assert.That(same.GetProperty("columns").GetProperty("aligned_target")[0].GetString(), Is.EqualTo("peptide"));
            Assert.That(same.GetProperty("residue_map").GetProperty("columns").GetProperty("mapped_residue")[0].GetString(),
                Is.EqualTo("p"));
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(),
                Does.StartWith("Pair 1 was refused by mzLib's PairwiseAligner:").And.Contain("MaxCells"));
        });
    }

    [Test]
    public void AnUnnamedCharacter_IsScoredAsX()
    {
        // The caveat: mzLib refuses only '-' and '.', so U, a digit or a space score as BLOSUM62's X.
        int Score(string pair) => Invoke(pair + "\n", "proteins", "align").GetProperty("columns").GetProperty("score")[0].GetInt32();
        int asX = Score("PEPTXDE\tPEPTXDE");

        Assert.Multiple(() =>
        {
            Assert.That(Score("PEPTUDE\tPEPTXDE"), Is.EqualTo(asX));
            Assert.That(Score("PEPT1DE\tPEPTXDE"), Is.EqualTo(asX));
            Assert.That(Score("PEPT DE\tPEPTXDE"), Is.EqualTo(asX));
        });
    }

    [Test]
    public void AnEmptySequence_FacesEveryResidueOfTheOtherWithAGap()
    {
        JsonElement data = Invoke("\tPEP\n", "proteins", "align");
        JsonElement map = data.GetProperty("residue_map").GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").GetProperty("aligned_source")[0].GetString(), Is.EqualTo("---"));
            Assert.That(data.GetProperty("columns").GetProperty("score")[0].GetInt32(), Is.Zero);
            Assert.That(Strings(map.GetProperty("sequence")), Is.EqualTo(new[] { "target", "target", "target" }));
            Assert.That(Ints(map.GetProperty("mapped_position")), Is.All.Null);
        });
    }

    [Test]
    public void RealIsoforms_TheExonIsoform2BLacksFacesAGap()
    {
        // DYNC1I2 (Q13409) and its isoform 2B (Q13409-2), from the committed human_extra.fasta.
        Dictionary<string, string> fasta = ReadFasta("human_extra.fasta");
        JsonElement data = Invoke($"{fasta["Q13409"]}\t{fasta["Q13409-2"]}\n", "proteins", "align");
        JsonElement columns = data.GetProperty("columns");
        JsonElement map = data.GetProperty("residue_map").GetProperty("columns");
        int?[] toTarget = Ints(map.GetProperty("mapped_position")).Take(638).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("source_length")[0].GetInt32(), Is.EqualTo(638));
            Assert.That(columns.GetProperty("target_length")[0].GetInt32(), Is.EqualTo(632));
            Assert.That(columns.GetProperty("identical")[0].GetInt32(), Is.EqualTo(632));
            Assert.That(columns.GetProperty("source_residues_facing_gap")[0].GetInt32(), Is.EqualTo(6));
            Assert.That(columns.GetProperty("aligned_target")[0].GetString(), Does.Contain("PESPI------VPPPM"));
            // Before the missing six residues positions agree; after them they are six apart.
            Assert.That(toTarget[74], Is.EqualTo(75));
            Assert.That(toTarget.Skip(75).Take(6), Is.All.Null);
            Assert.That(toTarget[81], Is.EqualTo(76));
            Assert.That(toTarget[637], Is.EqualTo(632));
        });
    }

    [Test]
    public void RealOrthologs_HumanAndBovineAlbuminDifferInResidueAndLength()
    {
        Dictionary<string, string> human = ReadFasta("human_extra.fasta");
        Dictionary<string, string> bovine = ReadFasta("contaminants.fasta");
        JsonElement columns = Invoke($"{human["P02768"]}\t{bovine["P02769"]}\n", "proteins", "align", "--no-residue-map")
            .GetProperty("columns");

        Assert.Multiple(() =>
        {
            Assert.That(columns.GetProperty("identical")[0].GetInt32(), Is.EqualTo(465));
            Assert.That(columns.GetProperty("substituted")[0].GetInt32(), Is.EqualTo(142));
            Assert.That(columns.GetProperty("source_residues_facing_gap")[0].GetInt32(), Is.EqualTo(2));
            Assert.That(columns.GetProperty("score")[0].GetInt32(), Is.EqualTo(2565));
        });
    }

    [TestCase("", new string[0], "No sequence pairs on stdin")]
    [TestCase("\n \n", new string[0], "No sequence pairs on stdin")]
    [TestCase("PEPTIDE\n", new string[0], "Pair 0 has 1 tab-separated field(s)")]
    [TestCase("PEPTIDE\tPEPTIDE\nA\tB\tC\n", new string[0], "Pair 1 has 3 tab-separated field(s)")]
    [TestCase("PEP-TIDE\tPEPTIDE\n", new string[0], "gap character")]
    [TestCase("PEPTIDE\tPEP.TIDE\n", new string[0], "gap character")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--gap-open", "-1" }, "--gap-open was refused by mzLib's PairwiseAligner")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--gap-extend", "0" }, "--gap-extend was refused by mzLib's PairwiseAligner")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--max-cells", "0" }, "--max-cells was refused by mzLib's PairwiseAligner")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--max-cells", "3000000000" }, "--max-cells was refused by mzLib's PairwiseAligner")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--max-cells", "lots" }, "--max-cells must be an integer")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--gap-open", "1.5" }, "--gap-open must be an integer")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--gap-open" }, "--gap-open was given without a value")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--end-gaps" }, "--end-gaps was given without a value")]
    [TestCase("PEPTIDE\tPEPTIDE\n", new[] { "--end-gaps", "both" }, "must be 'free' or 'charged'")]
    public void BadCalls_AreUsageErrors(string stdin, string[] options, string expected)
    {
        JsonElement error = InvokeExpectingError(stdin, new[] { "proteins", "align" }.Concat(options).ToArray());

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(expected));
    }

    // ---- plumbing -------------------------------------------------------------------------------

    private static Dictionary<string, string> ReadFasta(string name, [CallerFilePath] string thisFile = "")
    {
        string path = Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", "python", "tests",
            "fixtures", "proteins", name));
        var sequences = new Dictionary<string, string>();
        string? accession = null;
        foreach (string line in File.ReadLines(path))
        {
            if (line.StartsWith('>'))
            {
                accession = line.Split('|')[1];
                sequences[accession] = "";
            }
            else if (accession is not null)
            {
                sequences[accession] += line.Trim();
            }
        }
        return sequences;
    }

    private static string?[] Strings(JsonElement array) => array.EnumerateArray().Select(e => e.ValueKind == JsonValueKind.Null ? null : e.GetString()).ToArray();

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
