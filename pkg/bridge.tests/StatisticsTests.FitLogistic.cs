using System.Globalization;
using System.Text;
using System.Text.Json;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats fit-logistic</c>. Its numbers must be mzLib's LogisticRegression, which mzLib holds to
/// R 4.6.1's glm(family = binomial) with its own reference fixture (make_glm_lmm_fixtures.R). These
/// tests hold the WIRE to that same fixture, so a parsing slip, a sample re-ordering or a column
/// mix-up in the bridge shows as a disagreement with R, not as a silent change.
/// </summary>
public partial class StatisticsTests
{
    private static JsonElement FitLogistic(string stdin, params string[] options) =>
        Run(Statistics.FitLogistic, stdin, new[] { "stats", "fit-logistic" }.Concat(options).ToArray());

    /// <summary>
    /// mzLib's glm reference input as the verb's two TSVs: feature ids <c>f1</c>.. and the sample
    /// names its header already has (V1..V40) added, cells verbatim. With <paramref name="reverseDesign"/>
    /// the design rows are written in reverse, so a test also proves samples are matched by name.
    /// </summary>
    private static (string Responses, string Design) GlmReferenceTables(bool reverseDesign)
    {
        var (sampleHeader, responseRows) = ReadReference("glm_responses.tsv");
        var (designHeader, designRows) = ReadReference("glm_design.tsv");
        var responses = new StringBuilder("feature\t" + string.Join('\t', sampleHeader) + "\n");
        for (int f = 0; f < responseRows.Count; f++)
            responses.Append('f').Append(f + 1).Append('\t').Append(string.Join('\t', responseRows[f])).Append('\n');
        var design = new StringBuilder("sample\t" + string.Join('\t', designHeader) + "\n");
        IEnumerable<int> order = Enumerable.Range(0, designRows.Count);
        foreach (int s in reverseDesign ? order.Reverse() : order)
            design.Append(sampleHeader[s]).Append('\t').Append(string.Join('\t', designRows[s])).Append('\n');
        return (responses.ToString(), design.ToString());
    }

    private (string Responses, string Design) WriteGlmReference()
    {
        var (responses, design) = GlmReferenceTables(reverseDesign: true);
        return (WriteTemp("glm_responses_with_ids.tsv", responses), WriteTemp("glm_design_with_samples.tsv", design));
    }

    private static double[,] ReferenceMatrix(string name)
    {
        var (header, rows) = ReadReference(name);
        var m = new double[rows.Count, header.Length];
        for (int i = 0; i < rows.Count; i++)
            for (int j = 0; j < header.Length; j++)
                m[i, j] = Parse(rows[i][j]);
        return m;
    }

    // ---- agreement with R's glm ------------------------------------------------------------------

    [Test]
    public void FitLogistic_ReproducesRsGlm_ThroughTheWire()
    {
        var (responses, design) = WriteGlmReference();
        JsonElement root = FitLogistic("Intercept\nx\ngrp\n", "--responses", responses, "--design", design);

        var (h, rows) = ReadReference("glm_fit.tsv");
        double Col(string[] r, string c) => Parse(r[Array.IndexOf(h, c)]);
        string[] coefficients = { "Intercept", "x", "grp" };
        List<string?> status = Strings(root, "status");
        int fitted = 0, separated = 0;
        for (int f = 0; f < rows.Count; f++)
        {
            string reference = rows[f][Array.IndexOf(h, "status")];
            if (reference == "separated")
            {
                Assert.That(status[f], Is.EqualTo("separated"), $"feature {f}");
                separated++;
                continue;
            }
            if (reference == "not_fitted")
            {
                Assert.That(status[f], Is.Not.EqualTo("fitted"), $"feature {f}");
                continue;
            }
            Assert.That(status[f], Is.EqualTo("fitted"), $"feature {f}");
            fitted++;
            for (int j = 0; j < 3; j++)
            {
                AssertNear(Doubles(root, "estimate", coefficients[j])[f], Col(rows[f], $"b{j}"), $"feature {f} {coefficients[j]} estimate");
                // summary.glm's standard errors come from the weights of the step before the estimate;
                // mzLib's from the estimate itself. mzLib's own reference test allows 1e-6 for them.
                double se = Doubles(root, "standard_error", coefficients[j])[f];
                Assert.That(Math.Abs(se - Col(rows[f], $"se{j}")) / Col(rows[f], $"se{j}"), Is.LessThan(1e-6), $"feature {f} SE {j}");
            }
            double p = Doubles(root, "p_value", "x")[f];
            Assert.That(Math.Abs(p - Col(rows[f], "p1")) / Col(rows[f], "p1"), Is.LessThan(1e-6), $"feature {f} p(x)");
            AssertNear(Doubles(root, "deviance", "x")[f], Col(rows[f], "deviance"), $"feature {f} deviance");
        }
        Assert.That(fitted, Is.GreaterThan(rows.Count / 2), "most of mzLib's reference features are fittable");
        Assert.That(separated, Is.GreaterThan(0), "the reference holds separated features");

        // The deviance, iterations, observed and events of a feature are the same in every block.
        Assert.That(Doubles(root, "deviance", "Intercept"), Is.EqualTo(Doubles(root, "deviance", "grp")));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(3 * rows.Count));
        Assert.That(Strings(root, "coefficient").Take(rows.Count), Is.All.EqualTo("Intercept"), "blocks in stdin order");
    }

    [Test]
    public void FitLogistic_IsMzLibsFit_CellForCell_AtMzLibsDefaults()
    {
        // Absent options are mzLib's defaults: the wire must equal the overload called without them.
        var (responses, design) = WriteGlmReference();
        JsonElement root = FitLogistic("x\ngrp\n", "--responses", responses, "--design", design);
        LogisticModelFit direct = LogisticRegression.Fit(ReferenceMatrix("glm_responses.tsv"), ReferenceMatrix("glm_design.tsv"),
            new[] { "Intercept", "x", "grp" }, 1);

        Assert.That(root.GetProperty("max_iterations").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("tolerance").ValueKind, Is.EqualTo(JsonValueKind.Null));
        AssertSameAsDirect(root, direct, "x", 1);
        AssertSameAsDirect(root, direct, "grp", 2);

        // Given, both reach mzLib, and the echo says so.
        JsonElement tight = FitLogistic("x\n", "--responses", responses, "--design", design,
            "--max-iterations", "7", "--tolerance", "1e-6", "--threads", "-1");
        LogisticModelFit given = LogisticRegression.Fit(ReferenceMatrix("glm_responses.tsv"), ReferenceMatrix("glm_design.tsv"),
            new[] { "Intercept", "x", "grp" }, 1, maxIterations: 7, tolerance: 1e-6);
        Assert.That(tight.GetProperty("max_iterations").GetInt32(), Is.EqualTo(7));
        Assert.That(tight.GetProperty("tolerance").GetDouble(), Is.EqualTo(1e-6));
        Assert.That(tight.GetProperty("threads").GetInt32(), Is.EqualTo(-1));
        AssertSameAsDirect(tight, given, "x", 1);
    }

    private static void AssertSameAsDirect(JsonElement root, LogisticModelFit direct, string coefficient, int j)
    {
        int n = direct.FeatureCount;
        double[] Expected(Func<int, double> value) => Enumerable.Range(0, n).Select(value).ToArray();
        Assert.That(Doubles(root, "estimate", coefficient), Is.EqualTo(Expected(f => direct.Coefficient(f, j))));
        Assert.That(Doubles(root, "standard_error", coefficient), Is.EqualTo(Expected(f => direct.StandardError(f, j))));
        Assert.That(Doubles(root, "z", coefficient), Is.EqualTo(Expected(f => direct.ZValue(f, j))));
        Assert.That(Doubles(root, "p_value", coefficient), Is.EqualTo(Expected(f => direct.PValue(f, j))));
        Assert.That(Doubles(root, "deviance", coefficient), Is.EqualTo(Expected(f => direct.Deviance[f])));
        Assert.That(Doubles(root, "iterations", coefficient), Is.EqualTo(Expected(f => direct.Iterations[f])));
        Assert.That(Doubles(root, "observed", coefficient), Is.EqualTo(Expected(f => direct.Observed[f])));
        Assert.That(Doubles(root, "events", coefficient), Is.EqualTo(Expected(f => direct.Events[f])));
        var statuses = root.GetProperty("columns").GetProperty("status").EnumerateArray().Select(e => e.GetString()).ToList();
        var coefficients = root.GetProperty("columns").GetProperty("coefficient").EnumerateArray().Select(e => e.GetString()).ToList();
        Assert.That(Enumerable.Range(0, statuses.Count).Where(i => coefficients[i] == coefficient).Select(i => statuses[i]),
            Is.EqualTo(direct.Status.Select(Statistics.StatusName)));
    }

    // ---- statuses --------------------------------------------------------------------------------

    /// <summary>One feature per outcome LogisticRegression.Fit can give, over six samples and the design
    /// (intercept, x): too few observed, x constant where observed, every outcome 1, and a fittable one.</summary>
    private (string Responses, string Design) WriteLogisticOutcomes()
    {
        string responses = WriteTemp("outcomes.tsv",
            "feature\ts1\ts2\ts3\ts4\ts5\ts6\n" +
            "too_few\t1\t\t\t\t\t0\n" +
            "rank\t1\t0\t1\t\t\t\n" +
            "all_ones\t1\t1\t1\t1\t1\t1\n" +
            "fits\t0\t1\t0\t1\t1\t0\n");
        string design = WriteTemp("outcomes_design.tsv",
            "sample\tintercept\tx\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\ns4\t1\t1\ns5\t1\t2\ns6\t1\t3\n");
        return (responses, design);
    }

    [Test]
    public void FitLogistic_ReportsEveryOutcome_AndKeysStatusCountsToThem()
    {
        var (responses, design) = WriteLogisticOutcomes();
        JsonElement root = FitLogistic("x\n", "--responses", responses, "--design", design);
        Assert.That(Strings(root, "status"), Is.EqualTo(new[] { "too_few_observations", "rank_deficient", "separated", "fitted" }));
        Assert.That(Doubles(root, "iterations"), Is.EqualTo(new double[] { 0, 0, 0, Doubles(root, "iterations")[3] }),
            "no fit is attempted when a feature cannot be fitted");
        Assert.That(Doubles(root, "observed"), Is.EqualTo(new double[] { 2, 3, 6, 6 }));
        Assert.That(Doubles(root, "events"), Is.EqualTo(new double[] { 1, 2, 6, 3 }));
        Assert.That(Strings(root, "estimate").Take(3), Is.All.Null, "only a fitted feature carries statistics");

        var counts = root.GetProperty("status_counts").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetInt32());
        Assert.That(counts.Keys, Is.EqualTo(new[] { "fitted", "too_few_observations", "rank_deficient", "separated", "not_converged" }),
            "the outcomes LogisticRegression.Fit's body can return, in Statistics.LogisticFitStatuses order (stop-gap for mzLib#1487)");
        Assert.That(counts.Values.Sum(), Is.EqualTo(4));
        Assert.That(Caveats(root), Has.Some.StartsWith("3 of 4 features were not fitted (1 too_few_observations, 1 rank_deficient, 1 separated)"));

        // One iteration is not enough for the fittable feature: mzLib reports it not_converged.
        JsonElement once = FitLogistic("x\n", "--responses", responses, "--design", design, "--max-iterations", "1");
        Assert.That(Strings(once, "status")[3], Is.EqualTo("not_converged"));
        Assert.That(Doubles(once, "iterations")[3], Is.EqualTo(1));
        Assert.That(once.GetProperty("status_counts").GetProperty("not_converged").GetInt32(), Is.EqualTo(1));
    }

    [Test]
    public void FitLogistic_ReportsAFeatureTheDesignSeparates()
    {
        // x splits the outcomes perfectly: 0 below 1.5, 1 above.
        string responses = WriteTemp("sep.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\nsplit\t0\t0\t0\t1\t1\t1\nmixed\t0\t1\t0\t1\t0\t1\n");
        string design = WriteTemp("sep_design.tsv", "sample\tintercept\tx\ns1\t1\t0\ns2\t1\t1\ns3\t1\t1.4\ns4\t1\t1.6\ns5\t1\t2\ns6\t1\t3\n");
        JsonElement root = FitLogistic("x\n", "--responses", responses, "--design", design);
        Assert.That(Strings(root, "status"), Is.EqualTo(new[] { "separated", "fitted" }));
        Assert.That(Doubles(root, "iterations")[0], Is.GreaterThan(0), "a separated feature that IRLS ran on reports its iterations");
        Assert.That(Strings(root, "estimate")[0], Is.Null);
        Assert.That(Caveats(root), Has.Some.StartsWith("separated: the design splits").And.Some.Contains("1E-10"));
    }

    [Test]
    public void LogisticFitStatuses_AreFeatureFitStatusValues_WithFitted()
    {
        Assert.That(Statistics.LogisticFitStatuses, Does.Contain(FeatureFitStatus.Fitted));
        Assert.That(Statistics.LogisticFitStatuses, Is.Unique);
        Assert.That(Statistics.LogisticFitStatuses, Does.Not.Contain(FeatureFitStatus.TooFewGroups), "a mixed-model outcome");
    }

    // ---- refusals --------------------------------------------------------------------------------

    [Test]
    public void FitLogistic_RefusesBadInput_AsUsage()
    {
        var (responses, design) = WriteLogisticOutcomes();
        string Refusal(string stdin, params string[] options) =>
            Usage(Statistics.FitLogistic, stdin, new[] { "stats", "fit-logistic", "--responses", responses, "--design", design }
                .Concat(options).ToArray());

        Assert.That(Refusal(""), Does.Contain("No coefficient to test"));
        Assert.That(Refusal("y\n"), Does.Contain("No coefficient named 'y'"));
        Assert.That(Refusal("x\nx\n"), Does.Contain("more than once"));
        Assert.That(Refusal("x\n", "--max-iterations", "0"),
            Does.StartWith("mzLib's LogisticRegression.Fit refuses --max-iterations 0:"));
        Assert.That(Refusal("x\n", "--tolerance", "0"), Does.StartWith("mzLib's LogisticRegression.Fit refuses --tolerance 0:"));
        Assert.That(Refusal("x\n", "--tolerance", "NaN"), Does.StartWith("mzLib's LogisticRegression.Fit refuses --tolerance NaN:"));
        Assert.That(Refusal("x\n", "--max-iterations"), Does.Contain("without a value"));
        Assert.That(Refusal("x\n", "--tolerance", "tight"), Does.Contain("must be a number"));
        Assert.That(Refusal("x\n", "--threads", "0"), Does.Contain("--threads"));

        // Not a binary outcome: mzLib's refusal, told by the names the caller wrote.
        string proportions = WriteTemp("proportions.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\nA\t1\t0\t1\t0\t1\t0\nB\t1\t0\t0.5\t1\t0\t1\n");
        string message = Usage(Statistics.FitLogistic, "x\n", "stats", "fit-logistic", "--responses", proportions, "--design", design);
        Assert.That(message, Does.Contain("Response at feature 1, sample 2 is 0.5; expected 0, 1 or missing."));
        Assert.That(message, Does.EndWith("That is feature 'B', sample 's3'."));

        string collinear = WriteTemp("collinear.tsv", "sample\tintercept\tx\tx2\ns1\t1\t0\t0\ns2\t1\t0\t0\ns3\t1\t0\t0\ns4\t1\t1\t2\ns5\t1\t2\t4\ns6\t1\t3\t6\n");
        Assert.That(Usage(Statistics.FitLogistic, "x\n", "stats", "fit-logistic", "--responses", responses, "--design", collinear),
            Does.Contain("not of full column rank"));
    }

    // ---- the fixtures the Python docs run on -----------------------------------------------------

    /// <summary>
    /// mzLib's glm reference input as the verb's two TSVs (feature ids f1.. and a sample column added,
    /// cells verbatim) and R's glm output for it, verbatim; and the MALAT1 dilution series' detection
    /// table: 1 where <c>malat_dilution_log2_vs_fluc.tsv</c> has a value (FlashLFQ measured the oligo
    /// in that run), 0 where it is blank.
    /// </summary>
    [Test]
    public void LogisticFixtures_AreMzLibsGlmReference_AndTheDilutionSeriesDetection()
    {
        var (responses, design) = GlmReferenceTables(reverseDesign: false);
        AssertFixture("glm_reference_responses.tsv", responses);
        AssertFixture("glm_reference_design.tsv", design);
        AssertFixture("glm_reference_fit.tsv", File.ReadAllText(Reference("glm_fit.tsv")));

        var detected = new StringBuilder();
        foreach (var (line, i) in File.ReadAllLines(Path.Combine(StatsFixtures, "malat_dilution_log2_vs_fluc.tsv"))
                     .Where(l => l.Length > 0).Select((l, i) => (l, i)))
        {
            string[] cells = line.Split('\t');
            detected.Append(i == 0
                ? line
                : cells[0] + "\t" + string.Join('\t', cells.Skip(1).Select(c => c.Length == 0 ? "0" : "1")));
            detected.Append('\n');
        }
        AssertFixture("malat_dilution_detected.tsv", detected.ToString());
    }
}
