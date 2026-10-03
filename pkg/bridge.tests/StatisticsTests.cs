using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.Json;
using MzLibBridge;
using MathNet.Numerics.Statistics;
using Readers;

namespace MzLibBridge.Tests;

/// <summary>
/// The stats verbs: their numbers must be mzLib's StatisticalModels, which are themselves held to
/// limma 3.68.5 and metafor 5.2.1 by mzLib's own reference fixtures. These tests hold the WIRE to
/// those same fixtures, so a JSON round trip, a sample re-ordering or a column mix-up in the bridge
/// would show as a disagreement with limma, not as a silent change.
/// </summary>
[TestFixture]
[NonParallelizable] // swaps Console.In
public class StatisticsTests
{
    private const double Tolerance = 1e-8;

    private string _temp = "";

    [SetUp]
    public void CreateTemp()
    {
        _temp = Path.Combine(Path.GetTempPath(), "mzlib-bridge-stats-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_temp);
    }

    [TearDown]
    public void DeleteTemp()
    {
        if (Directory.Exists(_temp))
            Directory.Delete(_temp, recursive: true);
    }

    // ---- locations -------------------------------------------------------------------------------

    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    private static string MzLibTest(params string[] parts)
    {
        string path = Path.Combine(new[] { RepoRoot(), "code", "mzLib", "mzLib", "Test" }.Concat(parts).ToArray());
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present in the worktree: {path}");
        return path;
    }

    private static string Reference(string name) => MzLibTest("StatisticalModels", "ReferenceData", name);

    /// <summary>The derived real-data fixtures the Python guide and docstrings run on.</summary>
    internal static string StatsFixtures => Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "stats");

    // ---- running a verb --------------------------------------------------------------------------

    private static JsonElement Run(Func<Program.Arguments, object> verb, string stdin, params string[] args)
    {
        TextReader original = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            object result = verb(new Program.Arguments(args));
            return JsonDocument.Parse(JsonSerializer.Serialize(result, Program.JsonOptions)).RootElement;
        }
        finally
        {
            Console.SetIn(original);
        }
    }

    private static string Usage(Func<Program.Arguments, object> verb, string stdin, params string[] args)
    {
        TextReader original = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            return Assert.Throws<Program.UsageException>(() => verb(new Program.Arguments(args)))!.Message;
        }
        finally
        {
            Console.SetIn(original);
        }
    }

    private static double[] Doubles(JsonElement root, string column, string? coefficient = null)
    {
        JsonElement columns = root.GetProperty("columns");
        var values = columns.GetProperty(column).EnumerateArray().ToList();
        IEnumerable<int> rows = Enumerable.Range(0, values.Count);
        if (coefficient != null)
        {
            var coefficients = columns.GetProperty("coefficient").EnumerateArray().Select(e => e.GetString()).ToList();
            rows = rows.Where(i => coefficients[i] == coefficient);
        }
        return rows.Select(i => values[i].ValueKind == JsonValueKind.Null ? double.NaN : values[i].GetDouble()).ToArray();
    }

    private static double Parse(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static (string[] Header, List<string[]> Rows) ReadReference(string name)
    {
        string[] lines = File.ReadAllLines(Reference(name)).Where(l => l.Length > 0).ToArray();
        return (lines[0].Split('\t'), lines.Skip(1).Select(l => l.Split('\t')).ToList());
    }

    private static double[] ReferenceColumn(string name, string column)
    {
        var (header, rows) = ReadReference(name);
        int j = Array.IndexOf(header, column);
        Assert.That(j, Is.GreaterThanOrEqualTo(0), $"{name} has no column {column}");
        return rows.Select(r => Parse(r[j])).ToArray();
    }

    private static void AssertClose(IReadOnlyList<double> actual, IReadOnlyList<double> expected, string what)
    {
        Assert.That(actual.Count, Is.EqualTo(expected.Count), what);
        for (int i = 0; i < expected.Count; i++)
        {
            if (double.IsNaN(expected[i]))
            {
                Assert.That(actual[i], Is.NaN, $"{what}[{i}]");
                continue;
            }
            double scale = Math.Max(Math.Abs(expected[i]), 1e-300);
            Assert.That(Math.Abs(actual[i] - expected[i]) / scale, Is.LessThan(Tolerance),
                $"{what}[{i}]: {actual[i]:R} vs reference {expected[i]:R}");
        }
    }

    /// <summary>
    /// mzLib's limma reference tables as the verb's two TSVs. limma_responses.tsv has no feature id
    /// column and limma_design.tsv no sample column, so both are added; the cells are copied verbatim.
    /// The design rows are written REVERSED, so the test also proves the bridge matches samples by
    /// name rather than by position.
    /// </summary>
    private (string Responses, string Design) WriteLimmaReference()
    {
        var (sampleHeader, responseRows) = ReadReference("limma_responses.tsv");
        var (designHeader, designRows) = ReadReference("limma_design.tsv");

        var responses = new StringBuilder("feature\t" + string.Join('\t', sampleHeader) + "\n");
        for (int f = 0; f < responseRows.Count; f++)
            responses.Append("g").Append(f + 1).Append('\t').Append(string.Join('\t', responseRows[f])).Append('\n');

        var design = new StringBuilder("sample\t" + string.Join('\t', designHeader) + "\n");
        for (int s = designRows.Count - 1; s >= 0; s--)
            design.Append(sampleHeader[s]).Append('\t').Append(string.Join('\t', designRows[s])).Append('\n');

        string responsesPath = Path.Combine(_temp, "limma_responses_with_ids.tsv");
        string designPath = Path.Combine(_temp, "limma_design_with_samples.tsv");
        File.WriteAllText(responsesPath, responses.ToString());
        File.WriteAllText(designPath, design.ToString());
        return (responsesPath, designPath);
    }

    // ---- stats fit: agreement with limma ---------------------------------------------------------

    [TestCase(false, "notrend")]
    [TestCase(true, "trend")]
    public void Fit_ReproducesLimmaEBayesLegacy_ThroughTheWire(bool trend, string tag)
    {
        var (responses, design) = WriteLimmaReference();
        var args = new List<string> { "stats", "fit", "--responses", responses, "--design", design };
        if (trend) args.Add("--trend");

        JsonElement root = Run(Statistics.Fit, "age_decades\n", args.ToArray());

        Assert.That(root.GetProperty("residual_df_differ").GetBoolean(), Is.True,
            "the fixture omits missing values, so default limma would use its newer estimator");
        Assert.That(root.GetProperty("prior").GetProperty("trended").GetBoolean(), Is.EqualTo(trend));
        AssertClose(new[] { root.GetProperty("prior").GetProperty("df").GetDouble() },
            ReferenceColumn($"limma_ebayes_{tag}_prior.tsv", "df_prior"), "df.prior");

        string file = $"limma_ebayes_{tag}.tsv";
        AssertClose(Doubles(root, "prior_variance"), ReferenceColumn(file, "s2_prior"), "s2.prior");
        AssertClose(Doubles(root, "posterior_variance"), ReferenceColumn(file, "s2_post"), "s2.post");
        AssertClose(Doubles(root, "df_total"), ReferenceColumn(file, "df_total"), "df.total");
        AssertClose(Doubles(root, "t"), ReferenceColumn(file, "t_age"), "t");
        AssertClose(Doubles(root, "p_value"), ReferenceColumn(file, "p_age"), "p.value");
        AssertClose(Doubles(root, "bh_adjusted"), ReferenceColumn(file, "bh_age"), "BH");

        AssertClose(Doubles(root, "estimate"), ReferenceColumn("limma_fit.tsv", "coef_age"), "coef");
        AssertClose(Doubles(root, "sigma"), ReferenceColumn("limma_fit.tsv", "sigma"), "sigma");
        AssertClose(Doubles(root, "df_residual"), ReferenceColumn("limma_fit.tsv", "df_residual"), "df.residual");
        AssertClose(Doubles(root, "average_response"), ReferenceColumn("limma_fit.tsv", "amean"), "Amean");

        if (trend)
            Assert.That(root.GetProperty("prior").GetProperty("scale").ValueKind, Is.EqualTo(JsonValueKind.Null),
                "a trended prior varies per feature, so it is the prior_variance column, not one scale");
        else
            AssertClose(new[] { root.GetProperty("prior").GetProperty("scale").GetDouble() },
                new[] { ReferenceColumn(file, "s2_prior")[0] }, "scale");

        List<string> caveats = root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();
        Assert.That(caveats.Any(c => c.StartsWith("residual_df_differ:", StringComparison.Ordinal)), Is.True);
        Assert.That(caveats.Any(c => c.Contains("eBayes(legacy = TRUE)", StringComparison.Ordinal)), Is.True);
    }

    [Test]
    public void Fit_SeveralCoefficients_GiveOneBlockEach_AndTheSamePrior()
    {
        var (responses, design) = WriteLimmaReference();
        JsonElement both = Run(Statistics.Fit, "age_decades\nsex\n", "stats", "fit", "--responses", responses, "--design", design);
        JsonElement age = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design);
        JsonElement sex = Run(Statistics.Fit, "sex\n", "stats", "fit", "--responses", responses, "--design", design);

        Assert.That(both.GetProperty("row_count").GetInt32(), Is.EqualTo(800));
        Assert.That(Doubles(both, "t", "age_decades"), Is.EqualTo(Doubles(age, "t")));
        Assert.That(Doubles(both, "t", "sex"), Is.EqualTo(Doubles(sex, "t")));
        Assert.That(both.GetProperty("prior").GetProperty("df").GetDouble(),
            Is.EqualTo(sex.GetProperty("prior").GetProperty("df").GetDouble()),
            "the prior depends on the fit, not on the coefficient");
        Assert.That(both.GetProperty("tested_coefficients").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "age_decades", "sex" }));
        Assert.That(both.GetProperty("coefficient_names").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "Intercept", "age_decades", "sex" }));
    }

    [Test]
    public void Fit_IsIdenticalAtAnyThreadCount()
    {
        var (responses, design) = WriteLimmaReference();
        string one = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design, "--threads", "1").GetProperty("columns").GetRawText();
        string four = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design, "--threads", "4").GetProperty("columns").GetRawText();
        Assert.That(four, Is.EqualTo(one));
    }

    // ---- stats fit: statuses, zeros and the envelope ---------------------------------------------

    private (string, string) WriteSmall(string responses, string design)
    {
        string r = Path.Combine(_temp, "responses.tsv");
        string d = Path.Combine(_temp, "design.tsv");
        File.WriteAllText(r, responses);
        File.WriteAllText(d, design);
        return (r, d);
    }

    private const string TwoGroupDesign =
        "sample\tintercept\ttreated\n" +
        "a1\t1\t0\na2\t1\t0\na3\t1\t0\nb1\t1\t1\nb2\t1\t1\nb3\t1\t1\n";

    [Test]
    public void Fit_ReportsUnfittedFeatures_AsNullStatistics_WithACountedCaveat()
    {
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\n" +
            "up\t10\t10.2\t9.9\t12\t12.1\t11.8\n" +
            "flat\t20\t20.1\t19.8\t20.2\t19.9\t20\n" +
            "noise\t15\t15.5\t14.6\t15.2\t14.9\t15.4\n" +
            "sparse\t5\t\tNA\tNaN\t\t\n" +          // one observation: too few for two coefficients
            "onlya\t7\t7.1\t6.9\t\t\t\n",           // observed in one group only: 'treated' is not estimable
            TwoGroupDesign);

        JsonElement root = Run(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d);

        var status = root.GetProperty("columns").GetProperty("status").EnumerateArray().Select(e => e.GetString()).ToList();
        Assert.That(status, Is.EqualTo(new[] { "fitted", "fitted", "fitted", "too_few_observations", "rank_deficient" }));
        JsonElement counts = root.GetProperty("status_counts");
        Assert.That(counts.GetProperty("fitted").GetInt32(), Is.EqualTo(3));
        Assert.That(counts.GetProperty("too_few_observations").GetInt32(), Is.EqualTo(1));
        Assert.That(counts.GetProperty("rank_deficient").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("missing_count").GetInt32(), Is.EqualTo(8));
        Assert.That(root.GetProperty("columns").GetProperty("p_value")[3].ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("columns").GetProperty("df_residual")[4].ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("columns").GetProperty("observed")[3].GetInt32(), Is.EqualTo(1));
        Assert.That(Doubles(root, "estimate")[0], Is.EqualTo(1.9333333333333336).Within(1e-12));
        Assert.That(root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString())
            .Any(c => c!.StartsWith("2 of 5 features were not fitted (1 too_few_observations, 1 rank_deficient)", StringComparison.Ordinal)), Is.True);
    }

    [Test]
    public void Fit_CountsZeroCells_AsObservations_AndSaysSo()
    {
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\n" +
            "f1\t0\t10.2\t9.9\t12\t12.1\t11.8\n" +
            "f2\t20\t20.1\t19.8\t20.2\t19.9\t20\n",
            TwoGroupDesign);
        JsonElement root = Run(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d);
        Assert.That(root.GetProperty("zero_count").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("columns").GetProperty("observed")[0].GetInt32(), Is.EqualTo(6));
        Assert.That(root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString())
            .Any(c => c!.StartsWith("1 response cells are exactly 0", StringComparison.Ordinal)), Is.True);
    }

    [Test]
    public void Fit_EqualResidualDf_DoesNotRaiseTheLegacyCaveat()
    {
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\n" +
            "f1\t10\t10.2\t9.9\t12\t12.1\t11.8\n" +
            "f2\t20\t20.1\t19.8\t20.2\t19.9\t20\n" +
            "f3\t15\t15.5\t14.6\t15.2\t14.9\t15.4\n",
            TwoGroupDesign);
        JsonElement root = Run(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d);
        Assert.That(root.GetProperty("residual_df_differ").GetBoolean(), Is.False);
        Assert.That(root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString())
            .Any(c => c!.StartsWith("residual_df_differ", StringComparison.Ordinal)), Is.False);
        Assert.That(root.GetProperty("sample_names").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "a1", "a2", "a3", "b1", "b2", "b3" }));
        Assert.That(root.GetProperty("responses_file").GetString(), Is.EqualTo(Path.GetFullPath(r)));
        Assert.That(root.GetProperty("threads").GetInt32(), Is.EqualTo(1), "the bridge default is 1 (BULK.md), not mzLib's -1");
    }

    [Test]
    public void Fit_InfinitePriorDf_CrossesAsNull_WithTheFlag()
    {
        // Every feature has the same residual variance, so the variances are less dispersed than
        // sampling alone predicts and the prior df is infinite.
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\n" +
            "f1\t10\t11\t12\t20\t21\t22\n" +
            "f2\t30\t31\t32\t40\t41\t42\n" +
            "f3\t50\t51\t52\t50\t51\t52\n",
            TwoGroupDesign);
        JsonElement root = Run(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d);
        JsonElement prior = root.GetProperty("prior");
        Assert.That(prior.GetProperty("df_infinite").GetBoolean(), Is.True);
        Assert.That(prior.GetProperty("df").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString())
            .Any(c => c!.StartsWith("The prior df is infinite", StringComparison.Ordinal)), Is.True);
    }

    // ---- stats fit: usage errors -----------------------------------------------------------------

    [Test]
    public void Fit_RefusesBadInput_AsUsage_NamingTheProblem()
    {
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\n" +
            "f1\t10\t10.2\t9.9\t12\t12.1\t11.8\n" +
            "f2\t20\t20.1\t19.8\t20.2\t19.9\t20\n",
            TwoGroupDesign);

        Assert.That(Usage(Statistics.Fit, "", "stats", "fit", "--responses", r, "--design", d),
            Does.Contain("No coefficient to test").And.Contain("intercept, treated"));
        Assert.That(Usage(Statistics.Fit, "dose\n", "stats", "fit", "--responses", r, "--design", d),
            Does.Contain("No coefficient named 'dose'"));
        Assert.That(Usage(Statistics.Fit, "treated\ntreated\n", "stats", "fit", "--responses", r, "--design", d),
            Does.Contain("more than once"));
        Assert.That(Usage(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d, "--spline-basis", "3"),
            Does.Contain("--trend"));
        Assert.That(Usage(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d, "--trend", "--spline-basis", "0"),
            Does.Contain("1 or more"));
        Assert.That(Usage(Statistics.Fit, "treated\n", "stats", "fit", "--responses", Path.Combine(_temp, "missing.tsv"), "--design", d),
            Does.Contain("--responses file was not found"));
        Assert.That(Usage(Statistics.Fit, "treated\n", "stats", "fit", "--design", d),
            Does.Contain("--responses"));
    }

    [Test]
    public void Fit_RefusesMismatchedOrMalformedTables()
    {
        string Fit(string responses, string design)
        {
            var (r, d) = WriteSmall(responses, design);
            return Usage(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d);
        }

        const string goodResponses = "id\ta1\ta2\ta3\tb1\tb2\tb3\nf1\t10\t10\t9\t12\t12\t11\nf2\t20\t20\t19\t20\t19\t20\n";

        Assert.That(Fit(goodResponses.Replace("\tb3\n", "\tb4\n"), TwoGroupDesign),
            Does.Contain("Not in the design: b4").And.Contain("In the design but not in the responses: b3"));
        Assert.That(Fit(goodResponses.Replace("f2\t20\t20", "f2\t20\toops"), TwoGroupDesign),
            Does.Contain("'oops' is not a number"));
        Assert.That(Fit(goodResponses.Replace("f2\t", "f1\t"), TwoGroupDesign),
            Does.Contain("feature id 'f1' appears more than once"));
        Assert.That(Fit(goodResponses, TwoGroupDesign.Replace("b3\t1\t1", "b3\t1\tNA")),
            Does.Contain("must be a finite number"));
        Assert.That(Fit(goodResponses, TwoGroupDesign + "a1\t1\t0\n"),
            Does.Contain("Sample 'a1' appears more than once"));
        Assert.That(Fit("id\ta1\ta2\nf1\t1\n", "sample\ti\nA\t1\n"),
            Does.Contain("has 2 fields; the header has 3"));
        Assert.That(Fit("id\n", TwoGroupDesign), Does.Contain("needs a feature id column and at least one sample column"));
        Assert.That(Fit("id\ta1\n", TwoGroupDesign), Does.Contain("no feature rows"));
        // mzLib's own refusal of a rank-deficient design is reclassified as usage, without .NET's suffix.
        string rank = Fit(goodResponses,
            "sample\tintercept\ttreated\tcopy\na1\t1\t0\t0\na2\t1\t0\t0\na3\t1\t0\t0\nb1\t1\t1\t1\nb2\t1\t1\t1\nb3\t1\t1\t1\n");
        Assert.That(rank, Does.Contain("not of full column rank").And.Not.Contain("(Parameter"));
    }

    [Test]
    public void Fit_TooSparseForAPrior_IsUsage()
    {
        var (r, d) = WriteSmall(
            "id\ta1\ta2\ta3\tb1\tb2\tb3\nf1\t10\t10.2\t9.9\t12\t12.1\t11.8\nf2\t5\t\t\t\t\t\n",
            TwoGroupDesign);
        Assert.That(Usage(Statistics.Fit, "treated\n", "stats", "fit", "--responses", r, "--design", d),
            Does.Contain("at least 2 fitted features"));
    }

    // ---- stats adjust ----------------------------------------------------------------------------

    [Test]
    public void Adjust_ReproducesLimmasBenjaminiHochberg()
    {
        double[] p = ReferenceColumn("limma_ebayes_notrend.tsv", "p_age");
        string stdin = string.Join('\n', p.Select(v => v.ToString("R", CultureInfo.InvariantCulture))) + "\n";
        JsonElement root = Run(Statistics.Adjust, stdin, "stats", "adjust");
        AssertClose(Doubles(root, "bh_adjusted"), ReferenceColumn("limma_ebayes_notrend.tsv", "bh_age"), "BH");
        Assert.That(root.GetProperty("tested_count").GetInt32(), Is.EqualTo(p.Length));
    }

    [Test]
    public void Adjust_KeepsLinePositions_AndLeavesUntestedOut()
    {
        JsonElement root = Run(Statistics.Adjust, "0.01\n\nNA\n0.04\n0.03\n", "stats", "adjust");
        Assert.That(root.GetProperty("line_count").GetInt32(), Is.EqualTo(5));
        Assert.That(root.GetProperty("tested_count").GetInt32(), Is.EqualTo(3));
        double[] adjusted = Doubles(root, "bh_adjusted");
        Assert.That(adjusted[1], Is.NaN);
        Assert.That(adjusted[2], Is.NaN);
        Assert.That(adjusted[0], Is.EqualTo(0.03).Within(1e-15));
        Assert.That(adjusted[3], Is.EqualTo(0.04).Within(1e-15));
        Assert.That(adjusted[4], Is.EqualTo(0.04).Within(1e-15));
    }

    [Test]
    public void Adjust_RefusesBadInput_AsUsage()
    {
        Assert.That(Usage(Statistics.Adjust, "", "stats", "adjust"), Does.Contain("No p-values"));
        Assert.That(Usage(Statistics.Adjust, "0.1\nabc\n", "stats", "adjust"), Does.Contain("Line 2").And.Contain("not a number"));
        string outOfRange = Usage(Statistics.Adjust, "0.1\n1.5\n", "stats", "adjust");
        Assert.That(outOfRange, Does.Contain("outside [0, 1]").And.Not.Contain("(Parameter"));
    }

    // ---- stats meta ------------------------------------------------------------------------------

    [Test]
    public void Meta_ReproducesMetaforDerSimonianLaird()
    {
        var (inHeader, inRows) = ReadReference("dl_inputs.tsv");
        var (outHeader, outRows) = ReadReference("dl_results.tsv");
        int c = Array.IndexOf(inHeader, "case"), yi = Array.IndexOf(inHeader, "yi"), sei = Array.IndexOf(inHeader, "sei");
        string stdin = string.Join('\n', inRows.Select(r => $"{r[c]}\t{r[yi]}\t{r[sei]}")) + "\n";

        JsonElement root = Run(Statistics.Meta, stdin, "stats", "meta");
        var features = root.GetProperty("columns").GetProperty("feature").EnumerateArray().Select(e => e.GetString()!).ToList();
        Assert.That(features, Is.EqualTo(outRows.Select(r => r[Array.IndexOf(outHeader, "case")])),
            "features come back in first-appearance order");

        double Out(string[] row, string col) => Parse(row[Array.IndexOf(outHeader, col)]);
        foreach (var (column, reference) in new[] {
            ("estimate", "estimate"), ("standard_error", "se"), ("q", "q"), ("p_value", "pval"),
            ("confidence_low", "ci_lb"), ("confidence_high", "ci_ub") })
            AssertClose(Doubles(root, column), outRows.Select(r => Out(r, reference)).ToArray(), column);
        double[] tau2 = Doubles(root, "tau2"), i2 = Doubles(root, "i_squared");
        for (int i = 0; i < outRows.Count; i++)
        {
            Assert.That(tau2[i], Is.EqualTo(Out(outRows[i], "tau2")).Within(1e-12));
            Assert.That(i2[i], Is.EqualTo(Out(outRows[i], "i2")).Within(1e-12));
        }
        Assert.That(root.GetProperty("study_count").GetInt32(), Is.EqualTo(inRows.Count));
    }

    [Test]
    public void Meta_OneStudy_HasNoHeterogeneity_AndInterleavedFeaturesGroup()
    {
        JsonElement root = Run(Statistics.Meta, "B\t0.5\t0.1\nA\t1\t0.2\nB\t0.7\t0.1\n", "stats", "meta", "--confidence", "0.9");
        Assert.That(root.GetProperty("columns").GetProperty("feature").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "B", "A" }));
        Assert.That(root.GetProperty("columns").GetProperty("studies").EnumerateArray().Select(e => e.GetInt32()),
            Is.EqualTo(new[] { 2, 1 }));
        Assert.That(root.GetProperty("columns").GetProperty("q")[1].ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("columns").GetProperty("i_squared")[1].ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Doubles(root, "tau2")[1], Is.EqualTo(0));
        Assert.That(root.GetProperty("confidence").GetDouble(), Is.EqualTo(0.9));
        Assert.That(Doubles(root, "confidence_low")[1], Is.EqualTo(1 - 1.6448536269514722 * 0.2).Within(1e-12));
    }

    [Test]
    public void Meta_RefusesBadInput_AsUsage()
    {
        Assert.That(Usage(Statistics.Meta, "", "stats", "meta"), Does.Contain("No studies"));
        Assert.That(Usage(Statistics.Meta, "A\t1\n", "stats", "meta"), Does.Contain("expected 3"));
        Assert.That(Usage(Statistics.Meta, "\t1\t0.1\n", "stats", "meta"), Does.Contain("empty feature"));
        Assert.That(Usage(Statistics.Meta, "A\tx\t0.1\n", "stats", "meta"), Does.Contain("estimate must be a finite number"));
        Assert.That(Usage(Statistics.Meta, "A\t1\t0\n", "stats", "meta"), Does.Contain("Feature 'A'").And.Contain("finite and positive"));
        Assert.That(Usage(Statistics.Meta, "A\t1\t0.1\n", "stats", "meta", "--confidence", "1"), Does.Contain("between 0 and 1"));
        Assert.That(Usage(Statistics.Meta, "A\t1\t0.1\n", "stats", "meta", "--confidence"), Does.Contain("without a value"));
    }

    // ---- the reference fixtures the Python docs run on --------------------------------------------

    private static bool Regenerate => Environment.GetEnvironmentVariable("PYMZLIB_REGENERATE_STATS_FIXTURES") == "1";

    /// <summary>Writes <paramref name="expected"/> to the stats fixture folder when regenerating, then
    /// holds the committed file to it byte for byte (line endings aside).</summary>
    private static void AssertFixture(string name, string expected)
    {
        string path = Path.Combine(StatsFixtures, name);
        if (Regenerate)
        {
            Directory.CreateDirectory(StatsFixtures);
            File.WriteAllText(path, expected);
        }
        Assert.That(File.Exists(path), Is.True, $"{path} is missing; regenerate with PYMZLIB_REGENERATE_STATS_FIXTURES=1");
        Assert.That(File.ReadAllText(path).Replace("\r\n", "\n"), Is.EqualTo(expected.Replace("\r\n", "\n")), name);
    }

    /// <summary>
    /// mzLib's limma reference input as the verb's two TSVs (feature ids <c>g1</c>.. and a sample
    /// column added, cells verbatim), and limma's own eBayes output and metafor's DerSimonian-Laird
    /// inputs and results copied verbatim. The Python guide fits the first and compares with the
    /// second, so the claim "matches limma" is something a reader runs, not something they are told.
    /// </summary>
    [Test]
    public void ReferenceFixtures_AreMzLibsLimmaAndMetaforData()
    {
        var (sampleHeader, responseRows) = ReadReference("limma_responses.tsv");
        var (designHeader, designRows) = ReadReference("limma_design.tsv");

        var responses = new StringBuilder("feature\t" + string.Join('\t', sampleHeader) + "\n");
        for (int f = 0; f < responseRows.Count; f++)
            responses.Append("g").Append(f + 1).Append('\t').Append(string.Join('\t', responseRows[f])).Append('\n');
        var design = new StringBuilder("sample\t" + string.Join('\t', designHeader) + "\n");
        for (int s = 0; s < designRows.Count; s++)
            design.Append(sampleHeader[s]).Append('\t').Append(string.Join('\t', designRows[s])).Append('\n');

        AssertFixture("limma_reference_responses.tsv", responses.ToString());
        AssertFixture("limma_reference_design.tsv", design.ToString());
        foreach (string verbatim in new[] { "limma_ebayes_notrend.tsv", "dl_inputs.tsv", "dl_results.tsv", "PROVENANCE.txt" })
            AssertFixture(verbatim switch
            {
                "PROVENANCE.txt" => "limma_metafor_PROVENANCE.txt",
                _ when verbatim.StartsWith("dl_", StringComparison.Ordinal) => "metafor_" + verbatim,
                _ => verbatim,
            }, File.ReadAllText(Reference(verbatim)));
    }

    // ---- the real-data fixture -------------------------------------------------------------------

    /// <summary>
    /// The dilution series the Python guide analyses, derived from mzLib's own RNA test table
    /// (MetaMorpheus_RNA_AllQuantifiedOligos.tsv, mzLib #1388) by a rule written here and nowhere
    /// else:
    /// <list type="number">
    /// <item>keep the 27 runs of the 500/250/125 ng MALAT1 series at a fixed 500 ng FLuc and 50 ng
    /// 20mer-2 on the 147 min gradient, and drop runs labelled "dontuse";</item>
    /// <item>take log2(intensity), leaving a cell blank where FlashLFQ reported 0 (not detected);</item>
    /// <item>in each run, subtract the median log2 intensity of the oligos that map to FLuc alone. FLuc
    /// is loaded at 500 ng in every run, so this removes run-to-run loading (the raw run medians span
    /// about 7 log2 units) and leaves MALAT1 measured against a constant spike.</item>
    /// </list>
    /// The feature id is <c>&lt;protein groups&gt;:&lt;full sequence&gt;</c>. Set
    /// PYMZLIB_REGENERATE_STATS_FIXTURES=1 to rewrite the files from mzLib's table; otherwise this
    /// test fails if the committed files are not exactly what the rule produces.
    /// </summary>
    [Test]
    public void MalatDilutionFixture_IsLog2OfMzLibsOligoTable()
    {
        string table = MzLibTest("FileReadingTests", "ExternalFileTypes", "MetaMorpheus_RNA_AllQuantifiedOligos.tsv");
        var rows = new QuantifiedOligoFile(table).ToList();

        static int? MalatNanograms(string label)
        {
            if (label.Contains("dontuse", StringComparison.Ordinal) || !label.Contains("_20mer-2_50ng_147mingrad", StringComparison.Ordinal))
                return null;
            if (label.Contains("_1MALAT_500ng_1FLuc_500ng_", StringComparison.Ordinal)) return 500;
            if (label.Contains("_1MALAT_250ng_2FLuc_500ng_", StringComparison.Ordinal)) return 250;
            if (label.Contains("_1MALAT_125ng_4FLuc_500ng_", StringComparison.Ordinal)) return 125;
            return null;
        }

        List<string> samples = rows[0].Samples.Keys.Where(l => MalatNanograms(l) != null).ToList();
        Assert.That(samples, Has.Count.EqualTo(27));

        static double? Log2(double? intensity) => intensity is > 0 ? Math.Log2(intensity.Value) : null;
        Dictionary<string, double> flucMedian = samples.ToDictionary(sample => sample, sample => rows
            .Where(row => row.ProteinGroups == "FLuc")
            .Select(row => Log2(row.Samples[sample].Intensity))
            .Where(value => value.HasValue)
            .Select(value => value!.Value)
            .Median());

        var responses = new StringBuilder("oligo\t" + string.Join('\t', samples) + "\n");
        foreach (var row in rows)
        {
            responses.Append(row.ProteinGroups).Append(':').Append(row.Sequence);
            foreach (string sample in samples)
            {
                double? log2 = Log2(row.Samples[sample].Intensity);
                responses.Append('\t');
                if (log2.HasValue)
                    responses.Append((log2.Value - flucMedian[sample]).ToString("R", CultureInfo.InvariantCulture));
            }
            responses.Append('\n');
        }

        var design = new StringBuilder("run\tintercept\tmalat_250ng\tmalat_125ng\n");
        foreach (string sample in samples)
        {
            int ng = MalatNanograms(sample)!.Value;
            design.Append(sample).Append("\t1\t").Append(ng == 250 ? 1 : 0).Append('\t').Append(ng == 125 ? 1 : 0).Append('\n');
        }

        string responsesPath = Path.Combine(StatsFixtures, "malat_dilution_log2_vs_fluc.tsv");
        string designPath = Path.Combine(StatsFixtures, "malat_dilution_design.tsv");
        if (Environment.GetEnvironmentVariable("PYMZLIB_REGENERATE_STATS_FIXTURES") == "1")
        {
            Directory.CreateDirectory(StatsFixtures);
            File.WriteAllText(responsesPath, responses.ToString());
            File.WriteAllText(designPath, design.ToString());
        }

        Assert.That(File.ReadAllText(responsesPath).Replace("\r\n", "\n"), Is.EqualTo(responses.ToString()));
        Assert.That(File.ReadAllText(designPath).Replace("\r\n", "\n"), Is.EqualTo(design.ToString()));
    }
}
