using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats entrapment-fdp</c>. Its numbers must be mzLib's EntrapmentFdp.Sweep, which mzLib holds to
/// the entrapment project's reference script (results/fdp/fdp_analysis.py at commit 92f3108) on a
/// hand-built search in its own test (EntrapmentFdpTests.cs). These tests read that search and the
/// script's answers out of mzLib's test source, and hold the WIRE to them.
/// </summary>
public partial class StatisticsTests
{
    private static JsonElement RunEntrapmentFdp(params string[] options) =>
        Run(Statistics.Entrapment, "", new[] { "stats", "entrapment-fdp" }.Concat(options).ToArray());

    private static string Invariant(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    /// <summary>
    /// The hand-built search of mzLib's EntrapmentFdpTests: targets T1..T7 as (score, q), entrapments
    /// as (score, q, partner) with the partner a target's name or null, and fdp_analysis.py's output at
    /// r = 1 as (q, N_T, N_E, lower bound, combined, paired). Parsed from the test source, so the bridge
    /// test cannot drift from the numbers mzLib tests against.
    /// </summary>
    private static (Dictionary<string, (double Score, double Q)> Targets, List<(double Score, double Q, string? Partner)> Entrapments,
        List<double[]> Reference) EntrapmentReference()
    {
        string source = File.ReadAllText(MzLibTest("StatisticalModels", "EntrapmentFdpTests.cs"));
        var targets = Regex.Matches(source, @"ScoredIdentification (T\d) = new\(([\d.]+), ([\d.]+)\);")
            .ToDictionary(m => m.Groups[1].Value, m => (Parse(m.Groups[2].Value), Parse(m.Groups[3].Value)));
        string block = source[source.IndexOf("List<ScoredEntrapment> Entrapments =>", StringComparison.Ordinal)..];
        block = block[..block.IndexOf("];", StringComparison.Ordinal)];
        var entrapments = Regex.Matches(block, @"new\(([\d.]+), ([\d.]+), (T\d|null)\)")
            .Select(m => (Parse(m.Groups[1].Value), Parse(m.Groups[2].Value), m.Groups[3].Value == "null" ? null : m.Groups[3].Value))
            .ToList();
        string table = source[source.IndexOf("ReferenceR1 =", StringComparison.Ordinal)..];
        table = table[..table.IndexOf("];", StringComparison.Ordinal)];
        var reference = Regex.Matches(table, @"\(([\d.]+), (\d+), (\d+), ([\d.]+), ([\d.]+), ([\d.]+)\)")
            .Select(m => Enumerable.Range(1, 6).Select(i => Parse(m.Groups[i].Value)).ToArray())
            .ToList();
        Assert.That(targets, Has.Count.EqualTo(7));
        Assert.That(entrapments, Has.Count.EqualTo(4));
        Assert.That(reference, Has.Count.EqualTo(7));
        return (targets, entrapments, reference);
    }

    private static (string Targets, string Entrapments, string Reference) EntrapmentReferenceTables()
    {
        var (targets, entrapments, reference) = EntrapmentReference();
        var t = new StringBuilder("score\tq_value\n");
        foreach (var (score, q) in targets.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => kv.Value))
            t.Append(Invariant(score)).Append('\t').Append(Invariant(q)).Append('\n');
        var e = new StringBuilder("score\tq_value\tpartner_score\tpartner_q_value\n");
        foreach (var (score, q, partner) in entrapments)
        {
            e.Append(Invariant(score)).Append('\t').Append(Invariant(q)).Append('\t');
            e.Append(partner is null ? "\t" : Invariant(targets[partner].Score) + "\t" + Invariant(targets[partner].Q)).Append('\n');
        }
        var r = new StringBuilder("q_value_threshold\ttarget_count\tentrapment_count\tlower_bound\tcombined\tpaired\n");
        foreach (double[] row in reference)
            r.Append(string.Join('\t', row.Select(Invariant))).Append('\n');
        return (t.ToString(), e.ToString(), r.ToString());
    }

    private (string Targets, string Entrapments) WriteEntrapmentReference()
    {
        var (targets, entrapments, _) = EntrapmentReferenceTables();
        return (WriteTemp("targets.tsv", targets), WriteTemp("entrapments.tsv", entrapments));
    }

    // ---- agreement with the reference script -----------------------------------------------------

    [Test]
    public void EntrapmentFdp_ReproducesTheReferenceScript_ThroughTheWire()
    {
        var (_, _, reference) = EntrapmentReference();
        var (targets, entrapments) = WriteEntrapmentReference();
        string thresholds = string.Join(',', reference.Select(r => Invariant(r[0])));
        JsonElement root = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "1", "--thresholds", thresholds);

        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(reference.Count));
        Assert.That(Doubles(root, "q_value_threshold"), Is.EqualTo(reference.Select(r => r[0])));
        Assert.That(Doubles(root, "target_count"), Is.EqualTo(reference.Select(r => r[1])));
        Assert.That(Doubles(root, "entrapment_count"), Is.EqualTo(reference.Select(r => r[2])));
        // The script prints five decimals; mzLib's test compares at 1e-5.
        for (int i = 0; i < reference.Count; i++)
        {
            Assert.That(Doubles(root, "lower_bound")[i], Is.EqualTo(reference[i][3]).Within(1e-5), $"q {reference[i][0]}");
            Assert.That(Doubles(root, "combined")[i], Is.EqualTo(reference[i][4]).Within(1e-5), $"q {reference[i][0]}");
            Assert.That(Doubles(root, "paired")[i], Is.EqualTo(reference[i][5]).Within(1e-5), $"q {reference[i][0]}");
        }
        Assert.That(root.GetProperty("target_candidate_count").GetInt32(), Is.EqualTo(7));
        Assert.That(root.GetProperty("entrapment_candidate_count").GetInt32(), Is.EqualTo(4));
        Assert.That(root.GetProperty("partnered_count").GetInt32(), Is.EqualTo(3));
        Assert.That(root.GetProperty("foreign_count").GetInt32(), Is.EqualTo(0));
        Assert.That(root.GetProperty("thresholds_given").GetBoolean(), Is.True);
        Assert.That(root.GetProperty("empty_threshold_count").GetInt32(), Is.EqualTo(0));
    }

    [Test]
    public void EntrapmentFdp_IsMzLibsSweep_CellForCell_AtMzLibsDefaultThresholds()
    {
        var (t, e, _) = EntrapmentReference();
        var (targets, entrapments) = WriteEntrapmentReference();
        JsonElement root = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "1");

        var direct = StatisticalModels.EntrapmentFdp.Sweep(
            t.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => new ScoredIdentification(kv.Value.Score, kv.Value.Q)).ToList(),
            e.Select(x => new ScoredEntrapment(x.Score, x.Q,
                x.Partner is null ? null : new ScoredIdentification(t[x.Partner].Score, t[x.Partner].Q))).ToList(),
            1.0, StatisticalModels.EntrapmentFdp.DefaultQValueThresholds);

        Assert.That(root.GetProperty("thresholds").EnumerateArray().Select(x => x.GetDouble()),
            Is.EqualTo(StatisticalModels.EntrapmentFdp.DefaultQValueThresholds), "absent, the thresholds are mzLib's default, echoed");
        Assert.That(root.GetProperty("thresholds_given").GetBoolean(), Is.False);
        Assert.That(Doubles(root, "q_value_threshold"), Is.EqualTo(direct.Select(p => p.QValueThreshold)));
        Assert.That(Doubles(root, "target_count"), Is.EqualTo(direct.Select(p => (double)p.TargetCount)));
        Assert.That(Doubles(root, "entrapment_count"), Is.EqualTo(direct.Select(p => (double)p.EntrapmentCount)));
        Assert.That(Doubles(root, "lower_bound"), Is.EqualTo(direct.Select(p => p.LowerBound)));
        Assert.That(Doubles(root, "combined"), Is.EqualTo(direct.Select(p => p.Combined)));
        Assert.That(Doubles(root, "paired"), Is.EqualTo(direct.Select(p => p.Paired!.Value)));
        int empty = StatisticalModels.EntrapmentFdp.DefaultQValueThresholds.Count - direct.Count;
        Assert.That(root.GetProperty("empty_threshold_count").GetInt32(), Is.EqualTo(empty));
        Assert.That(empty, Is.EqualTo(0), "T1 is discovered at the first default threshold, 0.001");
    }

    [Test]
    public void EntrapmentFdp_AtRatioTwo_ScalesCombined_AndReportsNoPaired()
    {
        var (targets, entrapments) = WriteEntrapmentReference();
        JsonElement root = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "2", "--thresholds", "0.01");
        Assert.That(Doubles(root, "lower_bound")[0], Is.EqualTo(3.0 / 7).Within(1e-12));
        Assert.That(Doubles(root, "combined")[0], Is.EqualTo(3 * 1.5 / 7).Within(1e-12), "fdp_analysis.py at r = 2 gives 0.64286");
        Assert.That(Strings(root, "paired")[0], Is.Null);
        Assert.That(Caveats(root), Has.Some.EqualTo("paired is null at every threshold. Here the ratio is 2 and 0 of 4 entrapments are foreign."));
        // 0.9998 is not 1: mzLib reports paired only for r exactly 1.
        JsonElement almost = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "0.9998", "--thresholds", "0.01");
        Assert.That(Strings(almost, "paired")[0], Is.Null);
        JsonElement one = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "1.0", "--thresholds", "0.01");
        Assert.That(Strings(one, "paired")[0], Is.Not.Null);
    }

    // ---- foreign entrapments, partners, empty thresholds -----------------------------------------

    [Test]
    public void EntrapmentFdp_ForeignColumn_IsMzLibsIsForeign_WithItsDefaultWhenBlank()
    {
        string targets = WriteTemp("t.tsv", "q_value\tscore\n0.001\t100\n0.002\t90\n");
        string allForeign = WriteTemp("ef.tsv", "score\tq_value\tpartner_score\tpartner_q_value\tforeign\n85\t0.003\t\t\t1\n75\t0.004\t\t\ttrue\n");
        JsonElement foreign = RunEntrapmentFdp("--targets", targets, "--entrapments", allForeign, "--ratio", "1", "--thresholds", "0.005");
        Assert.That(Strings(foreign, "paired")[0], Is.Null, "every entrapment foreign: paired does not apply");
        Assert.That(foreign.GetProperty("foreign_count").GetInt32(), Is.EqualTo(2));

        string oneNative = WriteTemp("en.tsv", "score\tq_value\tpartner_score\tpartner_q_value\tforeign\n85\t0.003\t\t\t1\n75\t0.004\t\t\t\n");
        JsonElement native = RunEntrapmentFdp("--targets", targets, "--entrapments", oneNative, "--ratio", "1", "--thresholds", "0.005");
        Assert.That(native.GetProperty("foreign_count").GetInt32(), Is.EqualTo(1), "a blank foreign cell is mzLib's default, false");
        Assert.That(Doubles(native, "paired")[0], Is.EqualTo(2.0 * 2 / 4).Within(1e-12), "no partner anywhere: 2 N_E / (N_T + N_E)");
        Assert.That(Caveats(native), Has.Some.StartsWith("No entrapment has a partner, so paired is 2 N_E / (N_T + N_E)"));
        Assert.That(Doubles(native, "target_count")[0], Is.EqualTo(2), "columns are found by name, in any order");
    }

    [Test]
    public void EntrapmentFdp_IsNotCappedAtOne_AndLeavesOutThresholdsWithNoDiscoveries()
    {
        // Two entrapments that each beat a discovered partner, one target: paired = (2 + 0 + 2*2) / 3 = 2.
        string targets = WriteTemp("t.tsv", "score\tq_value\n50\t0.01\n");
        string entrapments = WriteTemp("e.tsv", "score\tq_value\tpartner_score\tpartner_q_value\n60\t0.005\t50\t0.01\n70\t0.004\t50\t0.01\n");
        JsonElement root = RunEntrapmentFdp("--targets", targets, "--entrapments", entrapments, "--ratio", "1", "--thresholds", "0.001,0.01");
        Assert.That(Doubles(root, "q_value_threshold"), Is.EqualTo(new[] { 0.01 }), "nothing is discovered at 0.001");
        Assert.That(Doubles(root, "paired")[0], Is.EqualTo(2.0).Within(1e-12));
        Assert.That(root.GetProperty("empty_threshold_count").GetInt32(), Is.EqualTo(1));
        Assert.That(Caveats(root), Has.Some.StartsWith("1 of 2 thresholds had no discoveries"));
    }

    // ---- refusals --------------------------------------------------------------------------------

    [Test]
    public void EntrapmentFdp_RefusesBadInput_AsUsage()
    {
        var (targets, entrapments) = WriteEntrapmentReference();
        string Refusal(params string[] options) => Usage(Statistics.Entrapment, "", new[] { "stats", "entrapment-fdp" }.Concat(options).ToArray());
        string[] files = { "--targets", targets, "--entrapments", entrapments };

        Assert.That(Refusal(files), Does.Contain("Missing required option --ratio"));
        Assert.That(Refusal(files.Concat(new[] { "--ratio" }).ToArray()), Does.Contain("without a value"));
        Assert.That(Refusal(files.Concat(new[] { "--ratio", "0" }).ToArray()),
            Does.StartWith("Option --ratio 0: The entrapment-to-target ratio must be a positive finite number."));
        Assert.That(Refusal(files.Concat(new[] { "--ratio", "1", "--thresholds", "0.01,x" }).ToArray()),
            Does.Contain("Option --thresholds entry 2 must be a finite number"));
        Assert.That(Refusal(files.Concat(new[] { "--ratio", "1", "--thresholds" }).ToArray()), Does.Contain("without a value"));
        Assert.That(Refusal("--targets", targets, "--ratio", "1"), Does.Contain("--entrapments"));

        string Bad(string targetsText, string entrapmentsText) => Refusal(
            "--targets", WriteTemp("bt.tsv", targetsText), "--entrapments", WriteTemp("be.tsv", entrapmentsText), "--ratio", "1");
        const string goodT = "score\tq_value\n1\t0.01\n";
        const string goodE = "score\tq_value\tpartner_score\tpartner_q_value\n1\t0.01\t\t\n";
        Assert.That(Bad("score\tqvalue\n1\t0.01\n", goodE), Does.Contain("Missing: q_value.").And.Contain("Not a column this file takes: qvalue."));
        Assert.That(Bad("id\tscore\tq_value\nA\t1\t0.01\n", goodE), Does.Contain("Not a column this file takes: id."));
        Assert.That(Bad("\tscore\tq_value\n0\t1\t0.01\n", goodE), Does.Contain("empty column"), "a pandas index column");
        Assert.That(Bad(goodT, "score\tq_value\tpartner_score\tpartner_q_value\tforiegn\n1\t0.01\t\t\t1\n"),
            Does.Contain("Not a column this file takes: foriegn."), "a misspelt optional column is refused, not defaulted");
        Assert.That(Bad("score\tq_value\n1\t\n", goodE), Does.Contain("Line 2 of the targets file: q_value must be a finite number"));
        Assert.That(Bad(goodT, "score\tq_value\tpartner_score\tpartner_q_value\n1\t0.01\t5\t\n"),
            Does.Contain("gives a partner score but no q_value"));
        Assert.That(Bad(goodT, "score\tq_value\tpartner_score\tpartner_q_value\tforeign\n1\t0.01\t\t\tyes\n"),
            Does.Contain("foreign is 'yes'"));
        Assert.That(Bad(goodT, "score\tq_value\tpartner_score\tpartner_q_value\n1\t0.01\t\n"), Does.Contain("has 3 fields"));
    }

    // ---- the fixtures the Python docs run on -----------------------------------------------------

    /// <summary>
    /// mzLib's hand-built entrapment search (EntrapmentFdpTests.cs) as the verb's two TSVs, each
    /// entrapment's partner written out as its target's score and q-value, and fdp_analysis.py's
    /// output for it at r = 1.
    /// </summary>
    [Test]
    public void EntrapmentFixtures_AreMzLibsReferenceSearch()
    {
        var (targets, entrapments, reference) = EntrapmentReferenceTables();
        AssertFixture("entrapment_reference_targets.tsv", targets);
        AssertFixture("entrapment_reference_entrapments.tsv", entrapments);
        AssertFixture("entrapment_reference_r1.tsv", reference);
    }
}
