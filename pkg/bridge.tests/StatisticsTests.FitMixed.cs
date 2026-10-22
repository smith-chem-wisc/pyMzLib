using System.Text;
using System.Text.Json;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats fit-mixed</c>. Its numbers must be mzLib's MixedModel, which mzLib holds to nlme 3.1.169's
/// lme (REML and ML, random intercept per group) with its own reference fixture
/// (make_glm_lmm_fixtures.R). These tests hold the WIRE to that same fixture, at mzLib's own
/// tolerance (nlme's optimizer limits agreement to about 1e-5).
/// </summary>
public partial class StatisticsTests
{
    private const double Lme = 1e-5;

    private static JsonElement FitMixed(string stdin, params string[] options) =>
        Run(Statistics.FitMixed, stdin, new[] { "stats", "fit-mixed" }.Concat(options).ToArray());

    /// <summary>
    /// mzLib's lme reference input as the verb's three TSVs: feature ids f1.., and the sample names its
    /// response header already has (V1..V29) added to the design and the groups, cells verbatim. With
    /// <paramref name="shuffle"/> the design rows are reversed and the group rows rotated, so a test also
    /// proves both are matched to the responses by sample name.
    /// </summary>
    private static (string Responses, string Design, string Groups) LmmReferenceTables(bool shuffle)
    {
        var (sampleHeader, responseRows) = ReadReference("lmm_responses.tsv");
        var (designHeader, designRows) = ReadReference("lmm_design.tsv");
        var (groupHeader, groupRows) = ReadReference("lmm_groups.tsv");
        var responses = new StringBuilder("feature\t" + string.Join('\t', sampleHeader) + "\n");
        for (int f = 0; f < responseRows.Count; f++)
            responses.Append('f').Append(f + 1).Append('\t').Append(string.Join('\t', responseRows[f])).Append('\n');
        int n = designRows.Count;
        var design = new StringBuilder("sample\t" + string.Join('\t', designHeader) + "\n");
        foreach (int s in shuffle ? Enumerable.Range(0, n).Reverse() : Enumerable.Range(0, n))
            design.Append(sampleHeader[s]).Append('\t').Append(string.Join('\t', designRows[s])).Append('\n');
        var groups = new StringBuilder("sample\t" + groupHeader[0] + "\n");
        foreach (int s in Enumerable.Range(0, n).Select(i => shuffle ? (i + 7) % n : i))
            groups.Append(sampleHeader[s]).Append('\t').Append(groupRows[s][0]).Append('\n');
        return (responses.ToString(), design.ToString(), groups.ToString());
    }

    private (string Responses, string Design, string Groups) WriteLmmReference()
    {
        var (responses, design, groups) = LmmReferenceTables(shuffle: true);
        return (WriteTemp("lmm_responses_with_ids.tsv", responses), WriteTemp("lmm_design_with_samples.tsv", design),
            WriteTemp("lmm_groups_with_samples.tsv", groups));
    }

    // ---- agreement with nlme ---------------------------------------------------------------------

    [TestCase("REML", "reml")]
    [TestCase("ML", "ml")]
    public void FitMixed_ReproducesNlme_ThroughTheWire(string method, string estimator)
    {
        var (responses, design, groups) = WriteLmmReference();
        JsonElement root = FitMixed("Intercept\nage\nwp\n", "--responses", responses, "--design", design, "--groups", groups,
            "--estimator", estimator);
        Assert.That(root.GetProperty("estimator").GetString(), Is.EqualTo(estimator));
        Assert.That(root.GetProperty("grouping").GetString(), Is.EqualTo("group"));
        Assert.That(root.GetProperty("group_names").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "d1", "d2", "d3", "d4", "d5", "d6", "d7" }), "in order of first appearance");

        var (h, all) = ReadReference("lmm_fit.tsv");
        var rows = all.Where(r => r[Array.IndexOf(h, "method")] == method).ToList();
        double Col(string[] r, string c) => Parse(r[Array.IndexOf(h, c)]);
        void Close(double actual, double expected, string what) =>
            Assert.That(Math.Abs(actual - expected) / Math.Max(Math.Abs(expected), 1e-300), Is.LessThan(Lme), $"{what}: {actual:R} vs {expected:R}");
        string[] coefficients = { "Intercept", "age", "wp" };
        List<string?> status = Strings(root, "status");
        int compared = 0;
        foreach (string[] r in rows)
        {
            int f = (int)Col(r, "feature");
            string reference = r[Array.IndexOf(h, "status")];
            if (reference == "r_error") continue;
            if (reference == "not_fitted") { Assert.That(status[f], Is.Not.EqualTo("fitted"), $"feature {f}"); continue; }
            Assert.That(status[f], Is.EqualTo("fitted"), $"feature {f}");
            compared++;
            double sigma2 = Col(r, "sigma2");
            for (int j = 0; j < 3; j++)
            {
                string c = coefficients[j];
                Close(Doubles(root, "estimate", c)[f], Col(r, $"b{j}"), $"{method} feature {f} {c} estimate");
                Close(Doubles(root, "standard_error", c)[f], Col(r, $"se{j}"), $"{method} feature {f} {c} SE");
                Assert.That(Doubles(root, "df", c)[f], Is.EqualTo(Col(r, $"df{j}")), $"{method} feature {f} {c} df");
            }
            // As mzLib's own test: a very small p-value magnifies a tiny difference in t, so compare log p.
            Close(Math.Log(Doubles(root, "p_value", "age")[f]), Math.Log(Col(r, "p1")), $"{method} feature {f} log p(age)");
            Close(Math.Log(Doubles(root, "p_value", "wp")[f]), Math.Log(Col(r, "p2")), $"{method} feature {f} log p(wp)");
            Close(Doubles(root, "residual_variance", "age")[f], sigma2, $"{method} feature {f} sigma2");
            // nlme cannot reach tau2 = 0 (it works on log tau), so compare on the scale of sigma2.
            Assert.That(Math.Abs(Doubles(root, "group_variance", "age")[f] - Col(r, "tau2")), Is.LessThan(Lme * sigma2),
                $"{method} feature {f} tau2");
            Close(Doubles(root, "log_likelihood", "age")[f], Col(r, "loglik"), $"{method} feature {f} logLik");
        }
        Assert.That(compared, Is.GreaterThan(rows.Count / 2), "most of mzLib's reference features are fittable");
        Assert.That(Caveats(root), Has.Some.StartsWith(estimator + ": "));
    }

    [Test]
    public void FitMixed_IsMzLibsFit_CellForCell_AtMzLibsDefaultEstimator()
    {
        var (responses, design, groups) = WriteLmmReference();
        JsonElement root = FitMixed("age\nwp\n", "--responses", responses, "--design", design, "--groups", groups);
        string[] labels = ReadReference("lmm_groups.tsv").Rows.Select(r => r[0]).ToArray();
        MixedModelFit direct = MixedModel.Fit(ReferenceMatrix("lmm_responses.tsv"), ReferenceMatrix("lmm_design.tsv"), labels,
            new[] { "Intercept", "age", "wp" }, maxThreads: 1);

        Assert.That(root.GetProperty("estimator").GetString(), Is.EqualTo(Statistics.EstimatorName(direct.Estimator)),
            "absent, the estimator is mzLib's default, as mzLib reports it");
        for (int j = 1; j <= 2; j++)
        {
            string c = j == 1 ? "age" : "wp";
            int n = direct.FeatureCount;
            double[] Expected(Func<int, double> value) => Enumerable.Range(0, n).Select(value).ToArray();
            Assert.That(Doubles(root, "estimate", c), Is.EqualTo(Expected(f => direct.Coefficient(f, j))));
            Assert.That(Doubles(root, "standard_error", c), Is.EqualTo(Expected(f => direct.StandardError(f, j))));
            Assert.That(Doubles(root, "df", c), Is.EqualTo(Expected(f => direct.DegreesOfFreedom(f, j))));
            Assert.That(Doubles(root, "t", c), Is.EqualTo(Expected(f => direct.TValue(f, j))));
            Assert.That(Doubles(root, "p_value", c), Is.EqualTo(Expected(f => direct.PValue(f, j))));
            Assert.That(Doubles(root, "residual_variance", c), Is.EqualTo(Expected(f => direct.ResidualVariance[f])));
            Assert.That(Doubles(root, "group_variance", c), Is.EqualTo(Expected(f => direct.GroupVariance[f])));
            Assert.That(Doubles(root, "log_likelihood", c), Is.EqualTo(Expected(f => direct.LogLikelihood[f])));
            Assert.That(Doubles(root, "observed", c), Is.EqualTo(Expected(f => direct.Observed[f])));
            Assert.That(Doubles(root, "groups", c), Is.EqualTo(Expected(f => direct.Groups[f])));
        }
        Assert.That(root.GetProperty("group_count").GetInt32(), Is.EqualTo(7));
    }

    // ---- statuses --------------------------------------------------------------------------------

    [Test]
    public void FitMixed_ReportsEveryOutcome_AndKeysStatusCountsToThem()
    {
        // Eight samples in four groups of two; design (intercept, x), x varying within groups.
        string responses = WriteTemp("mixed.tsv",
            "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n" +
            "too_few\t1.0\t\t\t\t\t\t\t2.0\n" +
            "rank\t1.0\t\t1.5\t\t2.2\t\t\t\n" +
            "one_each\t1.0\t\t\t1.9\t0.7\t\t\t2.2\n" +
            "replicates\t1.0\t1.000001\t5.0\t5.000003\t-2.0\t-2.000002\t3.0\t3.000001\n" +
            "fits\t1.1\t2.3\t0.4\t1.9\t2.8\t3.1\t0.2\t1.0\n");
        string design = WriteTemp("mixed_design.tsv",
            "sample\tintercept\tx\ns1\t1\t0\ns2\t1\t1\ns3\t1\t0\ns4\t1\t1\ns5\t1\t0\ns6\t1\t1\ns7\t1\t0\ns8\t1\t1\n");
        string groups = WriteTemp("mixed_groups.tsv", "sample\tdonor\ns1\ta\ns2\ta\ns3\tb\ns4\tb\ns5\tc\ns6\tc\ns7\td\ns8\td\n");
        JsonElement root = FitMixed("x\n", "--responses", responses, "--design", design, "--groups", groups);

        Assert.That(Strings(root, "status"),
            Is.EqualTo(new[] { "too_few_observations", "rank_deficient", "too_few_groups", "variance_ratio_at_limit", "fitted" }));
        Assert.That(Doubles(root, "groups"), Is.EqualTo(new double[] { 2, 3, 4, 4, 4 }),
            "one_each is seen in four groups, but in none twice");
        Assert.That(Strings(root, "estimate").Take(4), Is.All.Null, "only a fitted feature carries statistics");
        var counts = root.GetProperty("status_counts").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetInt32());
        Assert.That(counts.Keys, Is.EqualTo(new[]
            {
                "fitted", "too_few_observations", "rank_deficient", "too_few_groups", "variance_ratio_at_limit", "not_converged",
            }),
            "the outcomes MixedModel.Fit's body can return, in Statistics.MixedFitStatuses order (stop-gap for mzLib#1487)");
        Assert.That(counts.Values.Sum(), Is.EqualTo(5));
        Assert.That(root.GetProperty("grouping").GetString(), Is.EqualTo("donor"));
        Assert.That(Caveats(root), Has.Some.StartsWith(
            "4 of 5 features were not fitted (1 too_few_observations, 1 rank_deficient, 1 too_few_groups, 1 variance_ratio_at_limit)"));
    }

    /// <summary>
    /// mzLib's behaviour at d10017ce, pinned so a fix upstream is noticed: a feature the fixed effects
    /// fit exactly (here 0, 1, 0, 1 against x) leaves a residual sum of squares at rounding level, not
    /// exactly 0, so MixedModel's not_converged guard does not catch it. It is reported fitted, with a
    /// residual variance near 1e-16 and a t in the millions. The spec and the wire say so in a caveat.
    /// </summary>
    [Test]
    public void FitMixed_AnExactFit_IsReportedFitted_WithRoundingLevelVariance()
    {
        string responses = WriteTemp("exact.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nexact\t0\t1\t0\t1\t0\t1\t0\t1\n");
        string design = WriteTemp("exact_design.tsv",
            "sample\tintercept\tx\ns1\t1\t0\ns2\t1\t1\ns3\t1\t0\ns4\t1\t1\ns5\t1\t0\ns6\t1\t1\ns7\t1\t0\ns8\t1\t1\n");
        string groups = WriteTemp("exact_groups.tsv", "sample\tdonor\ns1\ta\ns2\ta\ns3\tb\ns4\tb\ns5\tc\ns6\tc\ns7\td\ns8\td\n");
        JsonElement root = FitMixed("x\n", "--responses", responses, "--design", design, "--groups", groups);
        Assert.That(Strings(root, "status")[0], Is.EqualTo("fitted"));
        Assert.That(Doubles(root, "residual_variance")[0], Is.LessThan(1e-12));
        Assert.That(Math.Abs(Doubles(root, "t")[0]), Is.GreaterThan(1e6));
        Assert.That(Caveats(root), Has.Some.StartsWith("A feature the fixed effects fit exactly"));
    }

    [Test]
    public void FitMixed_AGroupLevelCoefficientWithTooFewGroups_HasNoDf()
    {
        // Two groups and a covariate constant within each: containment df = G - 1 - 1 = 0.
        string responses = WriteTemp("between.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\nA\t1.0\t1.4\t0.9\t2.2\t2.7\t2.1\n");
        string design = WriteTemp("between_design.tsv", "sample\tintercept\tw\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\ns4\t1\t1\ns5\t1\t1\ns6\t1\t1\n");
        string groups = WriteTemp("between_groups.tsv", "sample\tbatch\ns1\tA\ns2\tA\ns3\tA\ns4\tB\ns5\tB\ns6\tB\n");
        JsonElement root = FitMixed("intercept\nw\n", "--responses", responses, "--design", design, "--groups", groups);
        Assert.That(Strings(root, "status"), Is.All.EqualTo("fitted"));
        Assert.That(Doubles(root, "df", "intercept")[0], Is.EqualTo(4), "within-group df: 6 - 2 groups - 0");
        Assert.That(Strings(root, "df")[1], Is.Null);
        Assert.That(Strings(root, "p_value")[1], Is.Null);
        Assert.That(Strings(root, "t")[1], Is.Not.Null, "t is reported; only its test has no df");
        Assert.That(Caveats(root), Has.Some.StartsWith("1 fitted rows have no positive containment df"));
    }

    [Test]
    public void MixedFitStatuses_AndEstimatorNames_CoverMzLibsValues()
    {
        Assert.That(Statistics.MixedFitStatuses, Does.Contain(FeatureFitStatus.Fitted));
        Assert.That(Statistics.MixedFitStatuses, Is.Unique);
        Assert.That(Statistics.MixedFitStatuses, Does.Not.Contain(FeatureFitStatus.Separated), "a logistic outcome");
        // Fails on the repin that adds a VarianceEstimator, before its C# name can reach the wire.
        Assert.That(Enum.GetValues<VarianceEstimator>().Select(Statistics.EstimatorName), Is.EqualTo(new[] { "reml", "ml" }));
    }

    // ---- refusals --------------------------------------------------------------------------------

    [Test]
    public void FitMixed_RefusesBadInput_AsUsage()
    {
        string responses = WriteTemp("r.tsv", "feature\ts1\ts2\ts3\ts4\nA\t1\t2\t3\t4\n");
        string design = WriteTemp("d.tsv", "sample\tintercept\ns1\t1\ns2\t1\ns3\t1\ns4\t1\n");
        string groups = WriteTemp("g.tsv", "sample\tday\ns1\tmon\ns2\tmon\ns3\ttue\ns4\ttue\n");
        string Refusal(string groupsPath, params string[] options) =>
            Usage(Statistics.FitMixed, "intercept\n", new[] { "stats", "fit-mixed", "--responses", responses, "--design", design,
                "--groups", groupsPath }.Concat(options).ToArray());

        Assert.That(Refusal(groups, "--estimator", "REML"), Does.Contain("must be 'reml' or 'ml'"));
        Assert.That(Refusal(groups, "--estimator"), Does.Contain("without a value"));
        Assert.That(Refusal(WriteTemp("g3.tsv", "sample\tday\tsite\ns1\tmon\tx\n")), Does.Contain("needs exactly two"));
        Assert.That(Refusal(WriteTemp("gmissing.tsv", "sample\tday\ns1\tmon\ns2\tmon\ns3\ttue\n")),
            Does.Contain("Not in the groups file: s4."));
        Assert.That(Refusal(WriteTemp("gextra.tsv", "sample\tday\ns1\tmon\ns2\tmon\ns3\ttue\ns4\ttue\ns5\twed\n")),
            Does.Contain("In the groups file but not in the responses: s5."));
        Assert.That(Refusal(WriteTemp("gdup.tsv", "sample\tday\ns1\tmon\ns1\tmon\ns3\ttue\ns4\ttue\n")),
            Does.Contain("more than once"));
        // An empty label is mzLib's refusal, told by the sample name the caller wrote.
        string blank = Refusal(WriteTemp("gblank.tsv", "sample\tday\ns1\tmon\ns2\tmon\ns3\t\ns4\ttue\n"));
        Assert.That(blank, Does.StartWith("Sample 2 has no group label.").And.EndWith("That is sample 's3'."));
        Assert.That(Usage(Statistics.FitMixed, "intercept\n", "stats", "fit-mixed", "--responses", responses, "--design", design),
            Does.Contain("--groups"));
    }

    // ---- the fixtures the Python docs run on -----------------------------------------------------

    /// <summary>
    /// mzLib's lme reference input as the verb's three TSVs (feature ids f1.. and sample columns added,
    /// cells verbatim) and nlme's output for it, verbatim; and the MALAT1 dilution series' acquisition
    /// day: each run's name begins with its date (2026_09_18_...), which is its group.
    /// </summary>
    [Test]
    public void MixedFixtures_AreMzLibsLmeReference_AndTheDilutionSeriesDays()
    {
        var (responses, design, groups) = LmmReferenceTables(shuffle: false);
        AssertFixture("lmm_reference_responses.tsv", responses);
        AssertFixture("lmm_reference_design.tsv", design);
        AssertFixture("lmm_reference_groups.tsv", groups);
        AssertFixture("lmm_reference_fit.tsv", File.ReadAllText(Reference("lmm_fit.tsv")));

        var runs = File.ReadAllLines(Path.Combine(StatsFixtures, "malat_dilution_design.tsv")).Where(l => l.Length > 0)
            .Skip(1).Select(l => l.Split('\t')[0]).ToList();
        Assert.That(runs, Has.All.Match(@"^\d{4}_\d{2}_\d{2}_"));
        var days = new StringBuilder("run\tday\n");
        foreach (string run in runs)
            days.Append(run).Append('\t').Append(run[..10]).Append('\n');
        AssertFixture("malat_dilution_day.tsv", days.ToString());
    }
}
