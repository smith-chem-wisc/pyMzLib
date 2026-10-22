using System.Globalization;
using System.Text;
using System.Text.Json;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats combine-p</c> and <c>stats correlate</c>. Their numbers must be mzLib's
/// PValueCombination and SpearmanCorrelation, which mzLib holds to R 4.6.1 (pchisq, qnorm,
/// cor.test) and ppcor 1.1 with its own reference cases. These tests hold the WIRE to those same
/// cases, so a parsing slip, a sample re-ordering or a lost sign in the bridge shows as a
/// disagreement with R, not as a silent change.
/// </summary>
public partial class StatisticsTests
{
    // ---- helpers ---------------------------------------------------------------------------------

    private string WriteTemp(string name, string content)
    {
        string path = Path.Combine(_temp, name);
        File.WriteAllText(path, content);
        return path;
    }

    /// <summary>Relative agreement to <see cref="Tolerance"/>; <paramref name="scaleFloor"/> 1 compares a
    /// unit-scale value (a z near 0) absolutely, as mzLib's own reference test does.</summary>
    private static void AssertNear(double actual, double expected, string what, double scaleFloor = 1e-300)
    {
        double scale = Math.Max(Math.Abs(expected), scaleFloor);
        Assert.That(Math.Abs(actual - expected) / scale, Is.LessThan(Tolerance),
            $"{what}: {actual:R} vs reference {expected:R}");
    }

    private static List<string?> Strings(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? null : e.ToString()).ToList();

    private static int[] Ints(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray().Select(e => e.GetInt32()).ToArray();

    private static List<string> Caveats(JsonElement root) =>
        root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();

    private static string R(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    private static JsonElement CombineP(string stdin, string method) =>
        Run(Statistics.CombineP, stdin, "stats", "combine-p", "--method", method);

    // ---- stats combine-p: agreement with R -------------------------------------------------------

    [Test]
    public void CombineP_ReproducesRsFisherAndStouffer_ThroughTheWire()
    {
        var (inHeader, inRows) = ReadReference("combination_inputs.tsv");
        var (outHeader, outRows) = ReadReference("combination_results.tsv");
        int c = Array.IndexOf(inHeader, "case"), p = Array.IndexOf(inHeader, "p"), w = Array.IndexOf(inHeader, "weight");
        string plain = string.Join('\n', inRows.Select(r => $"{r[c]}\t{r[p]}")) + "\n";
        string weightedLines = string.Join('\n', inRows.Select(r => $"{r[c]}\t{r[p]}\t{r[w]}")) + "\n";

        JsonElement fisher = CombineP(plain, "fisher");
        JsonElement stouffer = CombineP(plain, "stouffer");
        JsonElement weighted = CombineP(weightedLines, "stouffer");

        string[] cases = outRows.Select(r => r[Array.IndexOf(outHeader, "case")]).ToArray();
        foreach (JsonElement root in new[] { fisher, stouffer, weighted })
        {
            Assert.That(Strings(root, "feature"), Is.EqualTo(cases), "features come back in first-appearance order");
            Assert.That(Ints(root, "studies"), Is.EqualTo(cases.Select(k => inRows.Count(r => r[c] == k))));
            Assert.That(Ints(root, "omitted"), Is.All.EqualTo(0));
            Assert.That(root.GetProperty("line_count").GetInt32(), Is.EqualTo(inRows.Count));
            Assert.That(root.GetProperty("feature_count").GetInt32(), Is.EqualTo(cases.Length));
            Assert.That(root.GetProperty("signed").GetBoolean(), Is.False);
        }
        Assert.That(fisher.GetProperty("method").GetString(), Is.EqualTo("fisher"));
        Assert.That(stouffer.GetProperty("method").GetString(), Is.EqualTo("stouffer"));
        Assert.That(stouffer.GetProperty("weighted").GetBoolean(), Is.False);
        Assert.That(weighted.GetProperty("weighted").GetBoolean(), Is.True);

        double Out(int i, string column) => Parse(outRows[i][Array.IndexOf(outHeader, column)]);
        for (int i = 0; i < cases.Length; i++)
        {
            AssertNear(Doubles(fisher, "statistic")[i], Out(i, "fisher_x"), $"{cases[i]} Fisher X");
            AssertNear(Doubles(fisher, "p_value")[i], Out(i, "fisher_p"), $"{cases[i]} Fisher p");
            // Z is on a unit scale: where opposite tests cancel to about 1e-16, compare it absolutely.
            AssertNear(Doubles(stouffer, "statistic")[i], Out(i, "stouffer_z"), $"{cases[i]} Stouffer Z", scaleFloor: 1);
            AssertNear(Doubles(stouffer, "p_value")[i], Out(i, "stouffer_p"), $"{cases[i]} Stouffer p");
            AssertNear(Doubles(weighted, "statistic")[i], Out(i, "weighted_z"), $"{cases[i]} weighted Z", scaleFloor: 1);
            AssertNear(Doubles(weighted, "p_value")[i], Out(i, "weighted_p"), $"{cases[i]} weighted p");
        }
    }

    // ---- stats combine-p: signed, omitted, infinite ----------------------------------------------

    [Test]
    public void CombineP_Signed_ConvertsEachTwoSidedPValue_WithMzLibsOneSided()
    {
        // A: the same two-sided p, opposite effects, which must cancel. B: weighted, mixed directions.
        JsonElement unweighted = CombineP("A\t0.04\t\t1.2\nA\t0.04\t\t-0.5\nC\t0.01\t\t0.3\nC\t0.2\t\t-0.1\n", "stouffer");
        CombinedPValue a = PValueCombination.Stouffer(new[] { PValueCombination.OneSided(0.04, 1.2), PValueCombination.OneSided(0.04, -0.5) });
        CombinedPValue c = PValueCombination.Stouffer(new[] { PValueCombination.OneSided(0.01, 0.3), PValueCombination.OneSided(0.2, -0.1) });
        Assert.That(Doubles(unweighted, "statistic"), Is.EqualTo(new[] { a.Statistic, c.Statistic }));
        Assert.That(Doubles(unweighted, "p_value"), Is.EqualTo(new[] { a.PValue, c.PValue }));
        Assert.That(Doubles(unweighted, "p_value")[0], Is.EqualTo(0.5).Within(1e-12), "opposite effects cancel");
        Assert.That(unweighted.GetProperty("signed").GetBoolean(), Is.True);
        Assert.That(unweighted.GetProperty("weighted").GetBoolean(), Is.False);
        Assert.That(Caveats(unweighted), Has.Some.StartsWith("Each p-value was read as two-sided"));
        Assert.That(Caveats(unweighted), Has.Some.StartsWith("For a negative effect whose two-sided p-value is below about 2e-16"));

        JsonElement weighted = CombineP("B\t0.01\t2\t0.3\nB\t0.2\t3\t-0.1\n", "stouffer");
        CombinedPValue b = PValueCombination.Stouffer(
            new[] { PValueCombination.OneSided(0.01, 0.3), PValueCombination.OneSided(0.2, -0.1) }, new[] { 2.0, 3.0 });
        Assert.That(Doubles(weighted, "statistic")[0], Is.EqualTo(b.Statistic));
        Assert.That(Doubles(weighted, "p_value")[0], Is.EqualTo(b.PValue));
        Assert.That(weighted.GetProperty("weighted").GetBoolean(), Is.True);

        JsonElement fisher = CombineP("B\t0.01\t\t0.3\nB\t0.2\t\t-0.1\n", "fisher");
        CombinedPValue f = PValueCombination.Fisher(new[] { PValueCombination.OneSided(0.01, 0.3), PValueCombination.OneSided(0.2, -0.1) });
        Assert.That(Doubles(fisher, "p_value")[0], Is.EqualTo(f.PValue));
        Assert.That(Caveats(fisher), Has.None.StartsWith("For a negative effect"), "the rounding caveat is Stouffer's");
    }

    [Test]
    public void CombineP_OmitsUntestableTests_AndCountsThem()
    {
        JsonElement root = CombineP("A\t0.01\nA\t\nA\tNA\nB\tNaN\nC\t0.5\n", "fisher");
        Assert.That(Ints(root, "studies"), Is.EqualTo(new[] { 1, 0, 1 }));
        Assert.That(Ints(root, "omitted"), Is.EqualTo(new[] { 2, 1, 0 }));
        Assert.That(Doubles(root, "statistic")[0], Is.EqualTo(-2 * Math.Log(0.01)).Within(1e-12));
        Assert.That(Doubles(root, "p_value")[0], Is.EqualTo(0.01).Within(1e-12), "one test: Fisher returns its own p");
        Assert.That(Strings(root, "statistic")[1], Is.Null);
        Assert.That(Strings(root, "p_value")[1], Is.Null);
        Assert.That(Caveats(root), Has.Some.StartsWith("3 tests were omitted"));
        Assert.That(Caveats(root), Has.Some.StartsWith("1 feature has no test left to combine"));

        // Signed: an effect of 0, or none, has no direction, so mzLib's OneSided gives NaN and it is omitted.
        JsonElement signed = CombineP("A\t0.01\t\t0\nA\t0.02\t\t1\nA\t0.3\t\tNA\n", "stouffer");
        Assert.That(Ints(signed, "studies"), Is.EqualTo(new[] { 1 }));
        Assert.That(Ints(signed, "omitted"), Is.EqualTo(new[] { 2 }));
        Assert.That(Doubles(signed, "p_value")[0], Is.EqualTo(0.01).Within(1e-12), "the one test left: 0.02 two-sided, positive");
        Assert.That(Caveats(signed), Has.Some.Contains("or an effect that is 0, blank or NA"));
    }

    [Test]
    public void CombineP_InfiniteAndContradictoryStatistics_CrossAsNull_WithTheirPValue()
    {
        JsonElement fisher = CombineP("A\t0\nA\t0.5\n", "fisher");
        Assert.That(Strings(fisher, "statistic")[0], Is.Null);
        Assert.That(Doubles(fisher, "p_value")[0], Is.EqualTo(0));
        Assert.That(Caveats(fisher), Has.Some.StartsWith("1 feature has an infinite statistic"));

        JsonElement stouffer = CombineP("A\t0\nA\t0.5\nB\t1\nB\t0.3\nC\t0\nC\t1\n", "stouffer");
        Assert.That(Strings(stouffer, "statistic"), Is.EqualTo(new string?[] { null, null, null }));
        Assert.That(Strings(stouffer, "p_value")[0], Is.EqualTo("0"));
        Assert.That(Strings(stouffer, "p_value")[1], Is.EqualTo("1"));
        Assert.That(Strings(stouffer, "p_value")[2], Is.Null, "plus and minus infinity together are undefined");
        Assert.That(Ints(stouffer, "studies"), Is.EqualTo(new[] { 2, 2, 2 }), "mzLib still counts every finite input");
        Assert.That(Caveats(stouffer), Has.Some.StartsWith("1 feature holds both a p-value of 0 and a p-value of 1"));
        Assert.That(Caveats(stouffer), Has.Some.StartsWith("2 features have an infinite statistic"));
    }

    [Test]
    public void CombineP_InterleavedFeatures_ComeBackInFirstAppearanceOrder()
    {
        JsonElement root = CombineP("B\t0.1\nA\t0.2\nB\t0.3\n", "stouffer");
        Assert.That(Strings(root, "feature"), Is.EqualTo(new[] { "B", "A" }));
        Assert.That(Ints(root, "studies"), Is.EqualTo(new[] { 2, 1 }));
        Assert.That(root.GetProperty("column_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(Statistics.CombineColumns));
        Assert.That(Caveats(root), Has.Count.EqualTo(3), "only the static caveats apply");
    }

    [Test]
    public void CombineP_RefusesBadInput_AsUsage()
    {
        string Combine(string stdin, params string[] options) =>
            Usage(Statistics.CombineP, stdin, new[] { "stats", "combine-p" }.Concat(options).ToArray());

        Assert.That(Combine("A\t0.1\n"), Does.Contain("Missing required option --method"));
        Assert.That(Combine("A\t0.1\n", "--method"), Does.Contain("without a value"));
        Assert.That(Combine("A\t0.1\n", "--method", "Fisher"), Does.Contain("'fisher' or 'stouffer'"));
        Assert.That(Combine("", "--method", "fisher"), Does.Contain("No tests"));
        Assert.That(Combine("A\n", "--method", "fisher"), Does.Contain("1 tab-separated field;").And.Contain("expected 2 to 4"));
        Assert.That(Combine("A\t0.1\t1\t1\t1\n", "--method", "stouffer"), Does.Contain("5 tab-separated fields"));
        Assert.That(Combine("\t0.1\n", "--method", "fisher"), Does.Contain("empty feature name"));
        Assert.That(Combine("A\tsmall\n", "--method", "fisher"), Does.Contain("Line 1 p_value is not a number: 'small'"));
        Assert.That(Combine("A\t0.1\t2\nA\t0.2\n", "--method", "stouffer"),
            Does.Contain("Line 2 of stdin has no weight").And.Contain("every line or on none"));
        Assert.That(Combine("A\t0.1\t\t1\nA\t0.2\n", "--method", "stouffer"),
            Does.Contain("Line 2 of stdin has no effect field").And.Contain("every line"));
        Assert.That(Combine("A\t0.1\t2\n", "--method", "fisher"), Does.Contain("Fisher's method takes no weights"));
        Assert.That(Combine("A\t0.1\nB\t1.5\n", "--method", "fisher"),
            Does.Contain("Feature 'B'").And.Contain("outside [0, 1]").And.Not.Contain("(Parameter"));
        Assert.That(Combine("A\t1.5\t\t1\n", "--method", "stouffer"),
            Does.Contain("Feature 'A'").And.Contain("outside [0, 1]").And.Not.Contain("(Parameter"));
        Assert.That(Combine("A\t0.1\t0\n", "--method", "stouffer"),
            Does.Contain("Feature 'A'").And.Contain("finite and positive").And.Not.Contain("(Parameter"));
    }

    [Test]
    public void CombinationName_NamesEveryMzLibMethod_InSnakeCase()
    {
        var names = Enum.GetValues<PValueCombinationMethod>().Select(Statistics.CombinationName).ToList();
        Assert.That(names, Is.EqualTo(new[] { "fisher", "stouffer" }));
    }

    // ---- stats correlate: agreement with R -------------------------------------------------------

    /// <summary>One reference case as the verb's files: x the one feature row, y the trait, z the
    /// covariates. The trait and covariate rows are written REVERSED, so the test also proves the
    /// bridge matches samples by name rather than by position.</summary>
    private (string Responses, string Trait, string? Given) WriteCase(
        string name, IReadOnlyList<double> x, IReadOnlyList<double> y, IReadOnlyList<double[]>? z = null)
    {
        var samples = Enumerable.Range(1, x.Count).Select(i => $"s{i}").ToList();
        string responses = WriteTemp($"{name}_responses.tsv",
            "feature\t" + string.Join('\t', samples) + "\n" + name + "\t" + string.Join('\t', x.Select(R)) + "\n");
        var trait = new StringBuilder("sample\ty\n");
        for (int s = samples.Count - 1; s >= 0; s--)
            trait.Append(samples[s]).Append('\t').Append(R(y[s])).Append('\n');
        string? given = null;
        if (z != null)
        {
            var g = new StringBuilder("sample\t" + string.Join('\t', z.Select((_, j) => $"z{j + 1}")) + "\n");
            for (int s = samples.Count - 1; s >= 0; s--)
                g.Append(samples[s]).Append('\t').Append(string.Join('\t', z.Select(col => R(col[s])))).Append('\n');
            given = WriteTemp($"{name}_given.tsv", g.ToString());
        }
        return (responses, WriteTemp($"{name}_trait.tsv", trait.ToString()), given);
    }

    [Test]
    public void Correlate_ReproducesRsCorTest_ThroughTheWire()
    {
        var (inHeader, inRows) = ReadReference("spearman_inputs.tsv");
        var (outHeader, outRows) = ReadReference("spearman_results.tsv");
        int c = Array.IndexOf(inHeader, "case");
        foreach (string[] o in outRows)
        {
            string name = o[Array.IndexOf(outHeader, "case")];
            var rows = inRows.Where(r => r[c] == name).ToList();
            double[] x = rows.Select(r => Parse(r[Array.IndexOf(inHeader, "x")])).ToArray();
            double[] y = rows.Select(r => Parse(r[Array.IndexOf(inHeader, "y")])).ToArray();
            var (responses, trait, _) = WriteCase(name, x, y);

            JsonElement root = Run(Statistics.Correlate, "", "stats", "correlate", "--responses", responses, "--trait", trait);
            Assert.That(Strings(root, "method")[0], Is.EqualTo(o[Array.IndexOf(outHeader, "method")]), name);
            Assert.That(Ints(root, "n")[0], Is.EqualTo(x.Length), name);
            AssertNear(Doubles(root, "rho")[0], Parse(o[Array.IndexOf(outHeader, "rho")]), $"{name} rho");
            AssertNear(Doubles(root, "p_value")[0], Parse(o[Array.IndexOf(outHeader, "p")]), $"{name} p");
            Assert.That(root.GetProperty("partial").GetBoolean(), Is.False);
        }
    }

    [Test]
    public void Correlate_Given_ReproducesPpcor_ThroughTheWire()
    {
        var (inHeader, inRows) = ReadReference("partial_spearman_inputs.tsv");
        var (outHeader, outRows) = ReadReference("partial_spearman_results.tsv");
        int c = Array.IndexOf(inHeader, "case");
        double In(string[] r, string column) => Parse(r[Array.IndexOf(inHeader, column)]);
        foreach (string[] o in outRows)
        {
            string name = o[Array.IndexOf(outHeader, "case")];
            var rows = inRows.Where(r => r[c] == name).ToList();
            int k = int.Parse(o[Array.IndexOf(outHeader, "k")], CultureInfo.InvariantCulture);
            var z = Enumerable.Range(1, k).Select(j => rows.Select(r => In(r, $"z{j}")).ToArray()).ToList();
            var (responses, trait, given) = WriteCase(name, rows.Select(r => In(r, "x")).ToArray(), rows.Select(r => In(r, "y")).ToArray(), z);

            JsonElement root = Run(Statistics.Correlate, "", "stats", "correlate",
                "--responses", responses, "--trait", trait, "--given", given!);
            Assert.That(Strings(root, "method")[0], Is.EqualTo("asymptotic"), name);
            Assert.That(Ints(root, "n")[0], Is.EqualTo(int.Parse(o[Array.IndexOf(outHeader, "n")], CultureInfo.InvariantCulture)), name);
            AssertNear(Doubles(root, "rho")[0], Parse(o[Array.IndexOf(outHeader, "rho")]), $"{name} rho");
            AssertNear(Doubles(root, "p_value")[0], Parse(o[Array.IndexOf(outHeader, "p")]), $"{name} p");
            Assert.That(root.GetProperty("partial").GetBoolean(), Is.True);
            Assert.That(root.GetProperty("covariate_names").EnumerateArray().Select(e => e.GetString()),
                Is.EqualTo(Enumerable.Range(1, k).Select(j => $"z{j}")));
        }
    }

    // ---- stats correlate: missing values, not estimable, envelope --------------------------------

    private const string SixSampleResponses =
        "id\ts1\ts2\ts3\ts4\ts5\ts6\n" +
        "up\t1\t2\t3\t5\t4\t6\n" +
        "gappy\t3\t\t1\tNA\t2\t9\n" +
        "pair\t1\t\t\t\t\t2\n" +
        "flat\t4\t4\t4\t4\t4\t4\n";

    [Test]
    public void Correlate_OmitsMissingValuesPairwise_AsMzLibDoes()
    {
        string responses = WriteTemp("r.tsv", SixSampleResponses);
        // Rows out of order and one trait missing: s5 is in no feature's correlation.
        string trait = WriteTemp("t.tsv", "sample\tdose\ns6\t60\ns1\t10\ns3\t30\ns2\t20\ns5\tNA\ns4\t40\n");
        JsonElement root = Run(Statistics.Correlate, "", "stats", "correlate", "--responses", responses, "--trait", trait);

        double[] dose = { 10, 20, 30, 40, double.NaN, 60 };
        double[][] rows =
        {
            new double[] { 1, 2, 3, 5, 4, 6 },
            new[] { 3, double.NaN, 1, double.NaN, 2, 9 },
            new[] { 1, double.NaN, double.NaN, double.NaN, double.NaN, 2 },
            new double[] { 4, 4, 4, 4, 4, 4 },
        };
        SpearmanResult[] direct = rows.Select(x => SpearmanCorrelation.Correlate(x, dose)).ToArray();
        Assert.That(Ints(root, "n"), Is.EqualTo(direct.Select(r => r.N)));
        Assert.That(Ints(root, "n"), Is.EqualTo(new[] { 5, 3, 2, 5 }));
        Assert.That(Doubles(root, "rho")[0], Is.EqualTo(direct[0].Rho));
        Assert.That(Doubles(root, "p_value")[0], Is.EqualTo(direct[0].PValue));
        Assert.That(Doubles(root, "p_value")[1], Is.EqualTo(direct[1].PValue));
        Assert.That(Strings(root, "method"), Is.EqualTo(new[] { "exact", "exact", "not_estimable", "not_estimable" }));
        Assert.That(Strings(root, "has_ties"), Is.EqualTo(new[] { "False", "False", "False", "True" }),
            "a constant feature is tied; two samples are too few for mzLib to rank");

        Assert.That(root.GetProperty("trait_name").GetString(), Is.EqualTo("dose"));
        Assert.That(root.GetProperty("trait_missing_count").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("missing_count").GetInt32(), Is.EqualTo(6));
        Assert.That(root.GetProperty("given_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("given_incomplete_count").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("covariate_names").GetArrayLength(), Is.EqualTo(0));
        Assert.That(root.GetProperty("trait_file").GetString(), Is.EqualTo(Path.GetFullPath(trait)));
        Assert.That(root.GetProperty("sample_names").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "s1", "s2", "s3", "s4", "s5", "s6" }));
    }

    [Test]
    public void Correlate_NotEstimable_IsNull_WithACountedCaveat_AndMethodCounts()
    {
        string responses = WriteTemp("r.tsv", SixSampleResponses);
        string trait = WriteTemp("t.tsv", "sample\tdose\ns1\t1\ns2\t2\ns3\t3\ns4\t4\ns5\t5\ns6\t6\n");
        JsonElement root = Run(Statistics.Correlate, "", "stats", "correlate", "--responses", responses, "--trait", trait);

        Assert.That(Strings(root, "rho").Skip(2), Is.All.Null);
        Assert.That(Strings(root, "p_value").Skip(2), Is.All.Null);
        JsonElement counts = root.GetProperty("method_counts");
        Assert.That(counts.EnumerateObject().Select(p => p.Name), Is.EquivalentTo(new[] { "exact", "asymptotic", "not_estimable" }));
        Assert.That(counts.GetProperty("exact").GetInt32(), Is.EqualTo(2));
        Assert.That(counts.GetProperty("asymptotic").GetInt32(), Is.EqualTo(0));
        Assert.That(counts.GetProperty("not_estimable").GetInt32(), Is.EqualTo(2));
        Assert.That(Caveats(root), Has.Some.StartsWith("2 of 4 features are not_estimable: fewer than 3 complete samples"));
        Assert.That(Caveats(root), Has.Count.EqualTo(5), "four static caveats and the not_estimable count");
        Assert.That(Caveats(root), Has.Some.StartsWith("For 10 <= n < 1290 without ties"));
    }

    [Test]
    public void Correlate_Given_AsksMzLibsPartialForm_AndCountsIncompleteSamples()
    {
        string responses = WriteTemp("r.tsv",
            "id\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n" +
            "f1\t1\t3\t2\t5\t4\t7\t6\t8\n" +
            "f2\t8\t1\t6\t\t4\t3\t2\t5\n");
        string trait = WriteTemp("t.tsv", "sample\tage\ns1\t21\ns2\t25\ns3\t33\ns4\t38\ns5\t41\ns6\t52\ns7\t58\ns8\t64\n");
        string given = WriteTemp("g.tsv", "sample\tbatch\tload\ns8\t2\t1.1\ns7\t2\t0.9\ns6\t2\t1.4\ns5\t1\tNA\ns4\t1\t1.0\ns3\t1\t0.7\ns2\t1\t1.2\ns1\t1\t0.8\n");
        JsonElement root = Run(Statistics.Correlate, "", "stats", "correlate",
            "--responses", responses, "--trait", trait, "--given", given);

        double[] age = { 21, 25, 33, 38, 41, 52, 58, 64 };
        IReadOnlyList<double>[] covariates =
        {
            new double[] { 1, 1, 1, 1, 1, 2, 2, 2 },
            new[] { 0.8, 1.2, 0.7, 1.0, double.NaN, 1.4, 0.9, 1.1 },
        };
        SpearmanResult f1 = SpearmanCorrelation.PartialCorrelate(new double[] { 1, 3, 2, 5, 4, 7, 6, 8 }, age, covariates);
        SpearmanResult f2 = SpearmanCorrelation.PartialCorrelate(new[] { 8, 1, 6, double.NaN, 4, 3, 2, 5 }, age, covariates);
        Assert.That(Doubles(root, "rho"), Is.EqualTo(new[] { f1.Rho, f2.Rho }));
        Assert.That(Doubles(root, "p_value"), Is.EqualTo(new[] { f1.PValue, f2.PValue }));
        Assert.That(Ints(root, "n"), Is.EqualTo(new[] { 7, 6 }));
        Assert.That(Strings(root, "has_ties"), Is.EqualTo(new[] { "True", "True" }), "the batch covariate is tied");
        Assert.That(Strings(root, "method"), Is.All.EqualTo("asymptotic"));

        Assert.That(root.GetProperty("partial").GetBoolean(), Is.True);
        Assert.That(root.GetProperty("covariate_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "batch", "load" }));
        Assert.That(root.GetProperty("given_incomplete_count").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("given_file").GetString(), Is.EqualTo(Path.GetFullPath(given)));
        Assert.That(Caveats(root), Has.Count.EqualTo(5));
        Assert.That(Caveats(root), Has.Some.StartsWith("Partial Spearman"));
        Assert.That(Caveats(root), Has.Some.StartsWith("mzLib drops a covariate"));
        Assert.That(Caveats(root), Has.None.StartsWith("For 10 <= n < 1290"), "that caveat is the plain test's");
    }

    [Test]
    public void Correlate_RefusesBadInput_AsUsage()
    {
        string responses = WriteTemp("r.tsv", SixSampleResponses);
        string trait = WriteTemp("t.tsv", "sample\tdose\ns1\t1\ns2\t2\ns3\t3\ns4\t4\ns5\t5\ns6\t6\n");
        string Correlate(params string[] options) =>
            Usage(Statistics.Correlate, "", new[] { "stats", "correlate" }.Concat(options).ToArray());
        string WithTrait(string content) => Correlate("--responses", responses, "--trait", WriteTemp("bad.tsv", content));
        string WithGiven(string content) =>
            Correlate("--responses", responses, "--trait", trait, "--given", WriteTemp("given.tsv", content));

        Assert.That(Correlate("--responses", responses), Does.Contain("Missing required option --trait"));
        Assert.That(Correlate("--trait", trait), Does.Contain("Missing required option --responses"));
        Assert.That(Correlate("--responses", responses, "--trait", Path.Combine(_temp, "absent.tsv")),
            Does.Contain("--trait file was not found"));
        Assert.That(Correlate("--responses", responses, "--trait", trait, "--given"), Does.Contain("--given was given without a value"));
        Assert.That(Correlate("--responses", responses, "--trait", trait, "--given", Path.Combine(_temp, "absent.tsv")),
            Does.Contain("--given file was not found"));
        Assert.That(WithTrait("sample\tdose\tsex\ns1\t1\t0\ns2\t2\t0\ns3\t3\t0\ns4\t4\t1\ns5\t5\t1\ns6\t6\t1\n"),
            Does.Contain("has 2 value columns (dose, sex)").And.Contain("--given"));
        Assert.That(WithTrait("sample\ns1\n"), Does.Contain("needs a sample name column and at least one value column"));
        Assert.That(WithTrait("sample\tdose\ns1\t1\ns2\t2\ns3\t3\ns4\t4\ns5\t5\ns7\t6\n"),
            Does.Contain("Not in the trait file: s6").And.Contain("In the trait file but not in the responses: s7"));
        Assert.That(WithTrait("sample\tdose\ns1\t1\ns2\tlots\ns3\t3\ns4\t4\ns5\t5\ns6\t6\n"),
            Does.Contain("sample 's2', column 'dose' is not a number: 'lots'"));
        Assert.That(WithTrait("sample\tdose\ns1\t1\ns1\t2\ns3\t3\ns4\t4\ns5\t5\ns6\t6\n"),
            Does.Contain("Sample 's1' appears more than once in the trait file"));
        Assert.That(WithTrait("sample\tdose\ns1\t1\t9\n"), Does.Contain("has 3 fields; the header has 2"));
        Assert.That(WithGiven("sample\tbatch\ns1\t1\ns2\t1\ns3\t1\ns4\t2\ns5\t2\n"),
            Does.Contain("given file's samples must be exactly").And.Contain("Not in the given file: s6"));
        Assert.That(WithGiven("sample\tbatch\tbatch\ns1\t1\t1\n"), Does.Contain("column 'batch' appears more than once"));
    }

    [Test]
    public void SpearmanMethodName_NamesEveryMzLibMethod_InSnakeCase()
    {
        var names = Enum.GetValues<SpearmanPValueMethod>().Select(Statistics.SpearmanMethodName).ToList();
        Assert.That(names, Is.EqualTo(new[] { "exact", "asymptotic", "not_estimable" }));
    }

    // ---- the fixtures the Python docs run on -----------------------------------------------------

    /// <summary>
    /// mzLib's p-value combination reference cases and R's results for them, copied verbatim, and
    /// the provenance file that records how R produced them. The Python docstring combines the first
    /// and compares with the second.
    /// </summary>
    [Test]
    public void CombinationFixtures_AreMzLibsReferenceCases()
    {
        AssertFixture("combination_inputs.tsv", File.ReadAllText(Reference("combination_inputs.tsv")));
        AssertFixture("combination_results.tsv", File.ReadAllText(Reference("combination_results.tsv")));
        AssertFixture("combination_PROVENANCE.txt", File.ReadAllText(Reference("PROVENANCE_glm_lmm.txt")));
    }

    /// <summary>
    /// The per-sample tables the correlate examples read, each a column of a committed design:
    /// <list type="bullet">
    /// <item><c>malat_dilution_load_ng.tsv</c>: each run's MALAT1 load in ng, 250 where
    /// <c>malat_250ng</c> is 1, 125 where <c>malat_125ng</c> is 1, and 500 otherwise;</item>
    /// <item><c>limma_reference_age.tsv</c> and <c>limma_reference_sex.tsv</c>: the
    /// <c>age_decades</c> and <c>sex</c> columns of <c>limma_reference_design.tsv</c>, verbatim.</item>
    /// </list>
    /// </summary>
    [Test]
    public void CorrelateFixtures_AreColumnsOfTheCommittedDesigns()
    {
        static List<string[]> Rows(string name) =>
            File.ReadAllLines(Path.Combine(StatsFixtures, name)).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();

        var malat = Rows("malat_dilution_design.tsv");
        Assert.That(malat[0], Is.EqualTo(new[] { "run", "intercept", "malat_250ng", "malat_125ng" }));
        var load = new StringBuilder("run\tmalat_ng\n");
        foreach (string[] row in malat.Skip(1))
            load.Append(row[0]).Append('\t').Append(row[2] == "1" ? "250" : row[3] == "1" ? "125" : "500").Append('\n');
        AssertFixture("malat_dilution_load_ng.tsv", load.ToString());

        var limma = Rows("limma_reference_design.tsv");
        Assert.That(limma[0], Is.EqualTo(new[] { "sample", "Intercept", "age_decades", "sex" }));
        foreach (var (file, column) in new[] { ("limma_reference_age.tsv", 2), ("limma_reference_sex.tsv", 3) })
        {
            var table = new StringBuilder("sample\t" + limma[0][column] + "\n");
            foreach (string[] row in limma.Skip(1))
                table.Append(row[0]).Append('\t').Append(row[column]).Append('\n');
            AssertFixture(file, table.ToString());
        }
    }
}
