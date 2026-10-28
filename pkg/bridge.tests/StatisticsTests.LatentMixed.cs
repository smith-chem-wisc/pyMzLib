using System.Globalization;
using System.Text;
using System.Text.Json;
using MathNet.Numerics.Distributions;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats fit-latent-mixed</c> and <c>stats fit-link</c>. Their numbers must be mzLib's
/// LogitNormalMixedModel, which mzLib holds to lme4's lmer and glmer and to glmmTMB with its own
/// reference fixture (make_lnmm_fixtures.R, PROVENANCE_lnmm.txt). These tests hold the WIRE to that
/// same fixture at mzLib's own tolerances, and to mzLib's direct calls cell for cell.
/// </summary>
public partial class StatisticsTests
{
    /// <summary>
    /// The σ² mzLib's reference test gives in the counted cases, where R fits no run scatter
    /// (ReferenceComparisonLogitNormalTests.TinyResidualVariance at cc67b7a2).
    /// </summary>
    private const string TinyResidualVariance = "2e-7";

    private static JsonElement FitLatentMixed(string stdin, params string[] options) =>
        Run(Statistics.FitLatentMixed, stdin, new[] { "stats", "fit-latent-mixed" }.Concat(options).ToArray());

    private static JsonElement FitLinkVerb(params string[] options) =>
        Run(Statistics.FitLink, "", new[] { "stats", "fit-link" }.Concat(options).ToArray());

    private static double Expit(double x) => 1 / (1 + Math.Exp(-x));

    // ---- mzLib's reference input as the verb's tables --------------------------------------------

    /// <summary>The sample names: the V1..V23 header R wrote on every lnmm matrix.</summary>
    private static string[] LnmmSamples => ReadReference("lnmm_scatter_successes.tsv").Header;

    /// <summary>lnmm_design.tsv with the sample names added, cells verbatim; reversed with <paramref name="shuffle"/>.</summary>
    private static string LnmmDesignTable(bool shuffle)
    {
        var (header, rows) = ReadReference("lnmm_design.tsv");
        string[] samples = LnmmSamples;
        var text = new StringBuilder("sample\t" + string.Join('\t', header) + "\n");
        foreach (int s in shuffle ? Enumerable.Range(0, rows.Count).Reverse() : Enumerable.Range(0, rows.Count))
            text.Append(samples[s]).Append('\t').Append(string.Join('\t', rows[s])).Append('\n');
        return text.ToString();
    }

    /// <summary>lnmm_groups.tsv with the sample names added; rotated with <paramref name="shuffle"/>.</summary>
    private static string LnmmGroupsTable(bool shuffle)
    {
        var (header, rows) = ReadReference("lnmm_groups.tsv");
        string[] samples = LnmmSamples;
        int n = rows.Count;
        var text = new StringBuilder("sample\t" + header[0] + "\n");
        foreach (int s in Enumerable.Range(0, n).Select(i => shuffle ? (i + 5) % n : i))
            text.Append(samples[s]).Append('\t').Append(rows[s][0]).Append('\n');
        return text.ToString();
    }

    /// <summary>One fixture's successes and trials joined as k/n per cell ("NaN" = blank), feature ids f1...</summary>
    private static string LnmmCountsTable(string label)
    {
        var (header, successes) = ReadReference($"lnmm_{label}_successes.tsv");
        var (_, trials) = ReadReference($"lnmm_{label}_trials.tsv");
        var text = new StringBuilder("feature\t" + string.Join('\t', header) + "\n");
        for (int f = 0; f < successes.Count; f++)
        {
            text.Append('f').Append(f + 1);
            for (int s = 0; s < header.Length; s++)
                text.Append('\t').Append(successes[f][s] == "NaN" ? "" : successes[f][s] + "/" + trials[f][s]);
            text.Append('\n');
        }
        return text.ToString();
    }

    /// <summary>lnmm_measured_responses.tsv with feature ids f1.. added, cells verbatim ("NaN" kept: it is missing).</summary>
    private static string LnmmResponsesTable()
    {
        var (header, rows) = ReadReference("lnmm_measured_responses.tsv");
        var text = new StringBuilder("feature\t" + string.Join('\t', header) + "\n");
        for (int f = 0; f < rows.Count; f++)
            text.Append('f').Append(f + 1).Append('\t').Append(string.Join('\t', rows[f])).Append('\n');
        return text.ToString();
    }

    /// <summary>lnmm_counted_depth.tsv with the sample names added, cells verbatim.</summary>
    private static string LnmmDepthTable()
    {
        var (header, rows) = ReadReference("lnmm_counted_depth.tsv");
        string[] samples = LnmmSamples;
        var text = new StringBuilder("sample\t" + header[0] + "\n");
        for (int s = 0; s < rows.Count; s++)
            text.Append(samples[s]).Append('\t').Append(rows[s][0]).Append('\n');
        return text.ToString();
    }

    /// <summary>mzLib's TinyResidualVariance for every feature of a counted fixture.</summary>
    private static string LnmmResidualVarianceTable(int features)
    {
        var text = new StringBuilder("feature\tresidual_variance\n");
        for (int f = 0; f < features; f++)
            text.Append('f').Append(f + 1).Append('\t').Append(TinyResidualVariance).Append('\n');
        return text.ToString();
    }

    /// <summary>
    /// mzLib's measured-and-counted simulation (LogitNormalMixedModelTests.Simulated at cc67b7a2, one
    /// feature per seed): 16 animals of 4 runs, a group-level trait, latent = -0.5 + 1.0 trait + u + e,
    /// u ~ N(0, 0.6^2), e ~ N(0, 0.5^2); every other run measured, the rest counted, 5 trials each through
    /// the link (0.4, 0.6). One feature per seed, named seed1...; values written round-trip.
    /// </summary>
    private static (string Responses, string Counts, string Design, string Groups) SimulatedTables(params int[] seeds)
    {
        const int G = 16, k = 4, n = G * k;
        string Sample(int s) => $"s{s + 1}";
        var header = "feature\t" + string.Join('\t', Enumerable.Range(0, n).Select(Sample)) + "\n";
        var responses = new StringBuilder(header);
        var counts = new StringBuilder(header);
        foreach (int seed in seeds)
        {
            var rng = new Random(seed);
            var u = Enumerable.Range(0, G).Select(_ => Normal.Sample(rng, 0, 0.6)).ToArray();
            responses.Append("seed").Append(seed);
            counts.Append("seed").Append(seed);
            for (int s = 0; s < n; s++)
            {
                double latent = -0.5 + 1.0 * ((s / k) % 2) + u[s / k] + Normal.Sample(rng, 0, 0.5);
                if (s % 2 == 0)
                {
                    responses.Append('\t').Append(R(latent));
                    counts.Append('\t');
                    continue;
                }
                double g = Expit(0.4 + 0.6 * latent);
                int successes = Enumerable.Range(0, 5).Count(_ => rng.NextDouble() < g);
                responses.Append('\t');
                counts.Append('\t').Append(successes).Append("/5");
            }
            responses.Append('\n');
            counts.Append('\n');
        }
        var design = new StringBuilder("sample\tIntercept\ttrait\n");
        var groups = new StringBuilder("sample\tanimal\n");
        for (int s = 0; s < n; s++)
        {
            design.Append(Sample(s)).Append("\t1\t").Append((s / k) % 2).Append('\n');
            groups.Append(Sample(s)).Append("\ta").Append(s / k).Append('\n');
        }
        return (responses.ToString(), counts.ToString(), design.ToString(), groups.ToString());
    }

    /// <summary>
    /// mzLib's FitLink test cells (LogitNormalMixedModelTests.LinkCells at cc67b7a2, as
    /// FitLinkRecoversASimulatedLink calls it): the truth (-0.5, 0.5, 0.3), binomial, 4000 cells, seed 5.
    /// Values written round-trip.
    /// </summary>
    private static string LinkCellsTable(int count = 4000, int seed = 5, bool withDepth = true)
    {
        var rng = new Random(seed);
        var text = new StringBuilder("latent_mean\tlatent_sd\tsuccesses\ttrials" + (withDepth ? "\tdepth" : "") + "\n");
        for (int i = 0; i < count; i++)
        {
            double mean = Normal.Sample(rng, 0, 1.5), sd = 0.5, depth = Normal.Sample(rng, 0, 1);
            double latent = mean + Normal.Sample(rng, 0, sd);
            double g = Expit(-0.5 + 0.5 * latent + 0.3 * depth);
            int trials = 1 + rng.Next(6);
            int successes = Enumerable.Range(0, trials).Count(_ => rng.NextDouble() < g);
            text.Append(R(mean)).Append('\t').Append(R(sd)).Append('\t').Append(successes).Append('\t').Append(trials);
            if (withDepth) text.Append('\t').Append(R(depth));
            text.Append('\n');
        }
        return text.ToString();
    }

    private static List<LinkCell> ParseLinkCells(string table)
    {
        string[] lines = table.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        string[] header = lines[0].Split('\t');
        int depth = Array.IndexOf(header, "depth");
        return lines.Skip(1).Select(l => l.Split('\t')).Select(c => depth < 0
            ? new LinkCell(Parse(c[0]), Parse(c[1]), int.Parse(c[2]), int.Parse(c[3]))
            : new LinkCell(Parse(c[0]), Parse(c[1]), int.Parse(c[2]), int.Parse(c[3]), Parse(c[depth]))).ToList();
    }

    /// <summary>The verb's tables as the Python fixtures hold them: mzLib's own, rendered.</summary>
    [Test]
    public void LatentMixedFixtures_AreMzLibsLnmmReference()
    {
        AssertFixture("lnmm_reference_design.tsv", LnmmDesignTable(shuffle: false));
        AssertFixture("lnmm_reference_groups.tsv", LnmmGroupsTable(shuffle: false));
        AssertFixture("lnmm_scatter_counts.tsv", LnmmCountsTable("scatter"));
        AssertFixture("lnmm_scatter_fit.tsv", File.ReadAllText(Reference("lnmm_scatter_fit.tsv")));
        AssertFixture("lnmm_link_counts.tsv", LnmmCountsTable("counted_link"));
        AssertFixture("lnmm_link_depth.tsv", LnmmDepthTable());
        AssertFixture("lnmm_link_residual_variance.tsv", LnmmResidualVarianceTable(16));
        AssertFixture("lnmm_counted_fit.tsv", File.ReadAllText(Reference("lnmm_counted_fit.tsv")));
        AssertFixture("lnmm_PROVENANCE.txt", File.ReadAllText(Reference("PROVENANCE_lnmm.txt")));
        var (responses, counts, design, groups) = SimulatedTables(1, 2, 3);
        AssertSimulatedFixture("lnmm_simulated_responses.tsv", responses);
        AssertSimulatedFixture("lnmm_simulated_counts.tsv", counts);
        AssertSimulatedFixture("lnmm_simulated_design.tsv", design);
        AssertSimulatedFixture("lnmm_simulated_groups.tsv", groups);
        AssertSimulatedFixture("link_cells.tsv", LinkCellsTable());
    }

    /// <summary>
    /// <see cref="AssertFixture"/> for a table rendered from a simulation: every field equal, except that a
    /// non-integer number may differ by 1e-12 relative. MathNet's normal draws go through the platform's
    /// Math.Log and Math.Sin, whose last digits differ between Windows and Linux, so the committed file
    /// (rendered on Windows) cannot be byte for byte what Linux renders. Counts and trial numbers are exact.
    /// </summary>
    private static void AssertSimulatedFixture(string name, string expected)
    {
        if (Regenerate)
        {
            AssertFixture(name, expected);
            return;
        }
        string path = Path.Combine(StatsFixtures, name);
        Assert.That(File.Exists(path), Is.True, $"{path} is missing; regenerate with PYMZLIB_REGENERATE_STATS_FIXTURES=1");
        string[] committed = File.ReadAllText(path).Replace("\r\n", "\n").TrimEnd('\n').Split('\n');
        string[] rendered = expected.Replace("\r\n", "\n").TrimEnd('\n').Split('\n');
        Assert.That(committed, Has.Length.EqualTo(rendered.Length), $"{name}: line count");
        for (int line = 0; line < rendered.Length; line++)
        {
            string[] a = committed[line].Split('\t'), b = rendered[line].Split('\t');
            Assert.That(a, Has.Length.EqualTo(b.Length), $"{name} line {line + 1}: field count");
            for (int i = 0; i < b.Length; i++)
            {
                if (a[i] == b[i]) continue;
                bool numbers = double.TryParse(a[i], NumberStyles.Float, CultureInfo.InvariantCulture, out double x)
                    & double.TryParse(b[i], NumberStyles.Float, CultureInfo.InvariantCulture, out double y);
                Assert.That(numbers && b[i].Contains('.') && Math.Abs(x - y) <= 1e-12 * Math.Max(Math.Abs(y), 1e-300), Is.True,
                    $"{name} line {line + 1} field {i + 1}: committed '{a[i]}', rendered '{b[i]}'");
            }
        }
    }

    // ---- agreement with R, through the wire -------------------------------------------------------

    private (string Design, string Groups) WriteLnmmLayout() =>
        (WriteTemp("lnmm_design.tsv", LnmmDesignTable(shuffle: true)), WriteTemp("lnmm_groups.tsv", LnmmGroupsTable(shuffle: true)));

    /// <summary>One R reference row's columns by name.</summary>
    private static Func<string, double> LnmmRow(string[] header, string[] row) =>
        column => Parse(row[Array.IndexOf(header, column)]);

    /// <summary>
    /// With run scatter (σ² estimated) and every sample counted: glmmTMB, to the Laplace error, at mzLib's
    /// own tolerances (ReferenceComparisonLogitNormalTests.WithRunScatterItMatchesGlmmTmbToTheLaplaceError).
    /// </summary>
    [TestCase("REML", "reml")]
    [TestCase("ML", "ml")]
    public void FitLatentMixed_WithRunScatter_MatchesGlmmTmb_ThroughTheWire(string method, string estimator)
    {
        var (design, groups) = WriteLnmmLayout();
        string counts = WriteTemp("scatter_counts.tsv", LnmmCountsTable("scatter"));
        JsonElement root = FitLatentMixed("Intercept\ntrait\nx\n", "--counts", counts, "--design", design, "--groups", groups,
            "--link", "identity", "--estimator", estimator, "--threads", "-1");
        Assert.That(root.GetProperty("estimator").GetString(), Is.EqualTo(estimator));
        Assert.That(root.GetProperty("responses_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("measured_cells").GetInt32(), Is.EqualTo(0));

        var (h, all) = ReadReference("lnmm_scatter_fit.tsv");
        List<string?> status = Strings(root, "status");
        int compared = 0;
        foreach (string[] r in all.Where(r => r[Array.IndexOf(h, "method")] == method && r[Array.IndexOf(h, "status")] == "fitted"))
        {
            var col = LnmmRow(h, r);
            int f = (int)col("feature");
            Assert.That(status[f], Is.EqualTo("fitted"), $"feature {f}");
            compared++;
            Assert.That(Math.Abs(Doubles(root, "group_variance", "trait")[f] - col("tau2")), Is.LessThan(0.03 * col("tau2") + 0.003), $"{method} feature {f} tau2");
            Assert.That(Math.Abs(Doubles(root, "residual_variance", "trait")[f] - col("sigma2")), Is.LessThan(0.1 * col("sigma2") + 0.005), $"{method} feature {f} sigma2");
            Assert.That(Math.Abs(Doubles(root, "log_likelihood", "trait")[f] - col("loglik")), Is.LessThan(0.25), $"{method} feature {f} logLik");
            string[] coefficients = { "Intercept", "trait", "x" };
            for (int j = 0; j < 3; j++)
            {
                double se = col($"se{j}");
                Assert.That(Math.Abs(Doubles(root, "standard_error", coefficients[j])[f] - se) / se, Is.LessThan(0.05), $"{method} feature {f} SE {j}");
                if (method == "ML")
                    Assert.That(Math.Abs(Doubles(root, "estimate", coefficients[j])[f] - col($"b{j}")), Is.LessThan(0.03 * se), $"ML feature {f} b{j}");
            }
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(18), "glmmTMB fitted most of mzLib's reference features");
    }

    /// <summary>
    /// Every sample counted, σ² given (tiny): glmer's 25-point adaptive quadrature, through the identity link and
    /// one with an intercept, a slope of 0.6 and a depth term, at mzLib's own tolerances
    /// (ReferenceComparisonLogitNormalTests.WithEverySampleCountedItMatchesGlmerQuadrature).
    /// </summary>
    [TestCase("identity", "ML", "ml")]
    [TestCase("link", "ML", "ml")]
    [TestCase("link", "REML", "reml")]
    public void FitLatentMixed_EveryCounted_MatchesGlmerQuadrature_ThroughTheWire(string linkName, string method, string estimator)
    {
        var (design, groups) = WriteLnmmLayout();
        var (lh, links) = ReadReference("lnmm_counted_links.tsv");
        string[] l = links.Single(r => r[0] == linkName);
        string link = $"{l[Array.IndexOf(lh, "a")]},{l[Array.IndexOf(lh, "b")]},{l[Array.IndexOf(lh, "c")]}";
        string counts = WriteTemp("counted.tsv", LnmmCountsTable("counted_" + linkName));
        string depth = WriteTemp("depth.tsv", LnmmDepthTable());
        string sigma = WriteTemp("sigma.tsv", LnmmResidualVarianceTable(16));
        JsonElement root = FitLatentMixed("Intercept\ntrait\nx\n", "--counts", counts, "--design", design, "--groups", groups,
            "--link", link, "--depth", depth, "--residual-variance", sigma, "--estimator", estimator, "--threads", "-1");
        JsonElement echoed = root.GetProperty("link");
        Assert.That(echoed.GetProperty("slope").GetDouble(), Is.EqualTo(Parse(l[Array.IndexOf(lh, "b")])));
        Assert.That(echoed.GetProperty("depth_slope").GetDouble(), Is.EqualTo(Parse(l[Array.IndexOf(lh, "c")])));
        Assert.That(echoed.GetProperty("overdispersion").GetDouble(), Is.EqualTo(0), "not given: DetectionLink's own default");
        Assert.That(Strings(root, "residual_variance_given"), Is.All.EqualTo("True"));

        bool ml = method == "ML";
        double beta = ml ? 1e-5 : 2e-3, seTol = ml ? 5e-5 : 2e-3, variance = ml ? 5e-5 : 2e-3, logLik = ml ? 1e-5 : 5e-5;
        var (h, all) = ReadReference("lnmm_counted_fit.tsv");
        var rows = all.Where(r => r[Array.IndexOf(h, "link")] == linkName && r[Array.IndexOf(h, "method")] == method).ToList();
        string[] coefficients = { "Intercept", "trait", "x" };
        int compared = 0;
        foreach (string[] r in rows.Where(r => r[Array.IndexOf(h, "status")] == "fitted"))
        {
            var col = LnmmRow(h, r);
            int f = (int)col("feature");
            Assert.That(Strings(root, "status")[f], Is.EqualTo("fitted"), $"feature {f}");
            compared++;
            for (int j = 0; j < 3; j++)
            {
                double se = col($"se{j}");
                Assert.That(Math.Abs(Doubles(root, "estimate", coefficients[j])[f] - col($"b{j}")), Is.LessThan(beta * se), $"{linkName} {method} feature {f} b{j}");
                Assert.That(Math.Abs(Doubles(root, "standard_error", coefficients[j])[f] - se) / se, Is.LessThan(seTol), $"{linkName} {method} feature {f} se{j}");
            }
            Assert.That(Math.Abs(Doubles(root, "group_variance", "trait")[f] - col("tau2")), Is.LessThan(variance * Math.Max(1, col("tau2"))), $"{linkName} {method} feature {f} tau2");
            Assert.That(Math.Abs(Doubles(root, "log_likelihood", "trait")[f] - col("loglik")), Is.LessThan(logLik), $"{linkName} {method} feature {f} logLik");
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(rows.Count * 9 / 10), "glmer's quadrature is exact for most fixture features");
    }

    /// <summary>
    /// Nothing counted: lmer (REML) at mzLib's tolerance, and stats fit-mixed's model. Under reml the two
    /// verbs agree to 1e-5; under ml estimates and variances agree but the standard errors do not, which the
    /// caveats say.
    /// </summary>
    [Test]
    public void FitLatentMixed_WithNothingCounted_IsLmer_AndFitMixedsModel()
    {
        var (design, groups) = WriteLnmmLayout();
        string responses = WriteTemp("measured.tsv", LnmmResponsesTable());
        string[] coefficients = { "Intercept", "trait", "x" };
        JsonElement reml = FitLatentMixed("Intercept\ntrait\nx\n", "--responses", responses, "--design", design, "--groups", groups,
            "--link", "identity", "--threads", "-1");
        Assert.That(reml.GetProperty("counted_cells").GetInt32(), Is.EqualTo(0));
        Assert.That(reml.GetProperty("counts_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Caveats(reml), Has.Some.StartsWith("No cell was counted"));

        var (h, all) = ReadReference("lnmm_measured_fit.tsv");
        foreach (string[] r in all.Where(r => r[Array.IndexOf(h, "method")] == "REML" && r[Array.IndexOf(h, "status")] == "fitted"))
        {
            var col = LnmmRow(h, r);
            int f = (int)col("feature");
            for (int j = 0; j < 3; j++)
            {
                double se = col($"se{j}");
                Assert.That(Math.Abs(Doubles(reml, "estimate", coefficients[j])[f] - col($"b{j}")), Is.LessThan(1e-5 * se), $"feature {f} b{j}");
                Assert.That(Math.Abs(Doubles(reml, "standard_error", coefficients[j])[f] - se) / se, Is.LessThan(1e-4), $"feature {f} se{j}");
            }
            Assert.That(Math.Abs(Doubles(reml, "log_likelihood", "trait")[f] - col("loglik")), Is.LessThan(1e-5), $"feature {f} logLik");
        }

        JsonElement mixedReml = FitMixed("Intercept\ntrait\nx\n", "--responses", responses, "--design", design, "--groups", groups);
        JsonElement latentMl = FitLatentMixed("trait\n", "--responses", responses, "--design", design, "--groups", groups,
            "--link", "identity", "--estimator", "ml", "--threads", "-1");
        JsonElement mixedMl = FitMixed("trait\n", "--responses", responses, "--design", design, "--groups", groups, "--estimator", "ml");
        double Worst(JsonElement a, JsonElement b, string column, string? coefficient = null) =>
            Doubles(a, column, coefficient).Zip(Doubles(b, column, coefficient)).Max(p => Math.Abs(p.First - p.Second) / Math.Abs(p.Second));
        Assert.That(Worst(reml, mixedReml, "estimate", "trait"), Is.LessThan(1e-4));
        Assert.That(Worst(reml, mixedReml, "standard_error", "trait"), Is.LessThan(1e-4));
        Assert.That(Worst(latentMl, mixedMl, "estimate"), Is.LessThan(1e-4));
        Assert.That(Worst(latentMl, mixedMl, "residual_variance"), Is.LessThan(1e-4));
        Assert.That(Worst(latentMl, mixedMl, "standard_error"), Is.GreaterThan(1e-2),
            "under ml the standard errors differ: at the ML variances here, RSS / (observed - coefficients) in fit-mixed");
        Assert.That(Caveats(latentMl), Has.Some.Contains("standard errors do not"));
    }

    // ---- the wire is mzLib's fit, cell for cell -------------------------------------------------

    /// <summary>
    /// Measured and counted samples together (mzLib's own simulation), with the counts table's features and
    /// samples in a different order from the responses' and the design rows reversed: every column is
    /// mzLib's direct call, and mzLib's own claims hold through the wire.
    /// </summary>
    [Test]
    public void FitLatentMixed_MeasuredAndCounted_IsMzLibsFit_CellForCell()
    {
        var (responsesText, countsText, designText, groupsText) = SimulatedTables(1, 2, 3);
        string responses = WriteTemp("sim_responses.tsv", responsesText);
        string[] countLines = countsText.TrimEnd('\n').Split('\n');
        // Counts: features reversed, sample columns rotated.
        string[][] cells = countLines.Select(l => l.Split('\t')).ToArray();
        int n = cells[0].Length - 1;
        int[] order = Enumerable.Range(0, n).Select(i => 1 + (i + 17) % n).ToArray();
        var shuffled = new StringBuilder();
        foreach (string[] row in new[] { cells[0] }.Concat(cells.Skip(1).Reverse()))
            shuffled.Append(row[0]).Append('\t').Append(string.Join('\t', order.Select(i => row[i]))).Append('\n');
        string counts = WriteTemp("sim_counts.tsv", shuffled.ToString());
        string[] designLines = designText.TrimEnd('\n').Split('\n');
        string design = WriteTemp("sim_design.tsv", designLines[0] + "\n" + string.Join('\n', designLines.Skip(1).Reverse()) + "\n");
        string groups = WriteTemp("sim_groups.tsv", groupsText);

        JsonElement root = FitLatentMixed("trait\nIntercept\n", "--responses", responses, "--counts", counts, "--design", design,
            "--groups", groups, "--link", "0.4,0.6");
        Assert.That(root.GetProperty("measured_cells").GetInt32(), Is.EqualTo(96));
        Assert.That(root.GetProperty("counted_cells").GetInt32(), Is.EqualTo(96));
        Assert.That(root.GetProperty("missing_count").GetInt32(), Is.EqualTo(0));
        Assert.That(root.GetProperty("estimator").GetString(), Is.EqualTo("reml"), "absent: mzLib's default, as mzLib reports it");
        JsonElement link = root.GetProperty("link");
        Assert.That(new[] { "intercept", "slope", "depth_slope", "overdispersion" }.Select(k => link.GetProperty(k).GetDouble()),
            Is.EqualTo(new[] { 0.4, 0.6, 0, 0 }));
        Assert.That(Strings(root, "feature").Take(3), Is.EqualTo(new[] { "seed1", "seed2", "seed3" }), "the responses' order");

        // mzLib's direct call on the same simulation.
        double[,] y = ParseMatrix(responsesText, out _);
        var k = new BinomialCount?[3, 64];
        string[][] original = countsText.TrimEnd('\n').Split('\n').Skip(1).Select(l => l.Split('\t')).ToArray();
        for (int f = 0; f < 3; f++)
            for (int s = 0; s < 64; s++)
                if (original[f][s + 1].Length > 0)
                {
                    string[] kn = original[f][s + 1].Split('/');
                    k[f, s] = new BinomialCount(int.Parse(kn[0]), int.Parse(kn[1]));
                }
        var x = new double[64, 2];
        var labels = new string[64];
        for (int s = 0; s < 64; s++) { x[s, 0] = 1; x[s, 1] = (s / 4) % 2; labels[s] = $"a{s / 4}"; }
        LogitNormalMixedModelFit direct = LogitNormalMixedModel.Fit(y, k, x, labels, new DetectionLink(0.4, 0.6),
            coefficientNames: new[] { "Intercept", "trait" }, maxThreads: 1);
        foreach (var (c, j) in new[] { ("trait", 1), ("Intercept", 0) })
        {
            double[] Expected(Func<int, double> value) => Enumerable.Range(0, 3).Select(value).ToArray();
            Assert.That(Doubles(root, "estimate", c), Is.EqualTo(Expected(f => direct.Coefficient(f, j))));
            Assert.That(Doubles(root, "standard_error", c), Is.EqualTo(Expected(f => direct.StandardError(f, j))));
            Assert.That(Doubles(root, "df", c), Is.EqualTo(Expected(f => direct.DegreesOfFreedom(f, j))));
            Assert.That(Doubles(root, "t", c), Is.EqualTo(Expected(f => direct.TValue(f, j))));
            Assert.That(Doubles(root, "p_value", c), Is.EqualTo(Expected(f => direct.PValue(f, j))));
            Assert.That(Doubles(root, "residual_variance", c), Is.EqualTo(Expected(f => direct.ResidualVariance[f])));
            Assert.That(Doubles(root, "group_variance", c), Is.EqualTo(Expected(f => direct.GroupVariance[f])));
            Assert.That(Doubles(root, "log_likelihood", c), Is.EqualTo(Expected(f => direct.LogLikelihood[f])));
            Assert.That(Doubles(root, "measured", c), Is.EqualTo(Expected(f => direct.Measured[f])));
            Assert.That(Doubles(root, "counted", c), Is.EqualTo(Expected(f => direct.Counted[f])));
            Assert.That(Doubles(root, "groups", c), Is.EqualTo(Expected(f => direct.Groups[f])));
        }

        // mzLib's own claims (MeasuredAndCountedSamplesRecoverTheTraitEffect), through the wire: the trait is
        // recovered within 3 standard errors, tested on G - 2 = 14 df, and the counted samples add information.
        double[] trait = Doubles(root, "estimate", "trait"), se = Doubles(root, "standard_error", "trait");
        for (int f = 0; f < 3; f++)
            Assert.That(Math.Abs(trait[f] - 1.0), Is.LessThan(3 * se[f]), $"seed {f + 1}");
        Assert.That(Doubles(root, "df", "trait"), Is.All.EqualTo(14));
        JsonElement measuredOnly = FitMixed("trait\n", "--responses", responses, "--design", design, "--groups", groups);
        Assert.That(se[0], Is.LessThan(Doubles(measuredOnly, "standard_error")[0]));
    }

    private static double[,] ParseMatrix(string table, out List<string> ids)
    {
        string[][] rows = table.TrimEnd('\n').Split('\n').Skip(1).Select(l => l.Split('\t')).ToArray();
        ids = rows.Select(r => r[0]).ToList();
        var m = new double[rows.Length, rows[0].Length - 1];
        for (int i = 0; i < rows.Length; i++)
            for (int j = 1; j < rows[i].Length; j++)
                m[i, j - 1] = rows[i][j].Length == 0 ? double.NaN : Parse(rows[i][j]);
        return m;
    }

    // ---- statuses --------------------------------------------------------------------------------

    private const string EightSampleDesign =
        "sample\tintercept\tx\ns1\t1\t0\ns2\t1\t1\ns3\t1\t0\ns4\t1\t1\ns5\t1\t0\ns6\t1\t1\ns7\t1\t0\ns8\t1\t1\n";

    private const string EightSampleGroups = "sample\tdonor\ns1\ta\ns2\ta\ns3\tb\ns4\tb\ns5\tc\ns6\tc\ns7\td\ns8\td\n";

    [Test]
    public void FitLatentMixed_ReportsEveryOutcome_AndKeysStatusCountsToThem()
    {
        // Eight samples in four groups of two; design (intercept, x), x varying within groups.
        string responses = WriteTemp("latent.tsv",
            "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n" +
            "too_few\t1.0\t\t\t\t\t\t\t\n" +
            "rank\t1.0\t\t1.5\t\t2.2\t\t\t\n" +
            "one_each\t1.0\t\t\t1.9\t0.7\t\t\t2.2\n" +
            "one_sided\t\t\t\t\t\t\t\t\n" +
            "fits\t1.1\t2.3\t0.4\t1.9\t2.8\t3.1\t0.2\t1.0\n");
        string counts = WriteTemp("latent_counts.tsv",
            "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n" +
            "too_few\t\t\t\t\t\t\t\t2/4\n" +
            "rank\t\t\t\t\t\t\t\t\n" +
            "one_each\t\t\t\t\t\t\t\t\n" +
            "one_sided\t0/3\t0/3\t0/3\t0/3\t3/3\t3/3\t3/3\t3/3\n" +
            "fits\t\t\t\t\t\t\t\t\n");
        string design = WriteTemp("latent_design.tsv", EightSampleDesign);
        string groups = WriteTemp("latent_groups.tsv", EightSampleGroups);
        JsonElement root = FitLatentMixed("x\n", "--responses", responses, "--counts", counts, "--design", design, "--groups", groups,
            "--link", "identity");

        Assert.That(Strings(root, "status"),
            Is.EqualTo(new[] { "too_few_observations", "rank_deficient", "too_few_groups", "separated", "fitted" }));
        Assert.That(Ints(root, "measured"), Is.EqualTo(new[] { 1, 3, 4, 0, 8 }));
        Assert.That(Ints(root, "counted"), Is.EqualTo(new[] { 1, 0, 0, 8, 0 }));
        Assert.That(Strings(root, "estimate").Take(4), Is.All.Null, "only a fitted feature carries statistics");
        Assert.That(Strings(root, "residual_variance_given"), Is.All.EqualTo("False"));
        var statusCounts = root.GetProperty("status_counts").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetInt32());
        Assert.That(statusCounts.Keys, Is.EqualTo(new[]
            {
                "fitted", "too_few_observations", "rank_deficient", "separated", "not_converged", "too_few_groups",
            }),
            "the outcomes LogitNormalMixedModel.Fit's body can return, in FeatureFitStatus order (stop-gap for mzLib#1487)");
        Assert.That(statusCounts.Values.Sum(), Is.EqualTo(5));
        Assert.That(root.GetProperty("missing_count").GetInt32(), Is.EqualTo(40 - 16 - 9));
        Assert.That(Caveats(root), Has.Some.StartsWith(
            "4 of 5 features were not fitted (1 too_few_observations, 1 rank_deficient, 1 separated, 1 too_few_groups)"));
    }

    [Test]
    public void LatentMixedFitStatuses_CoverMzLibsOutcomes()
    {
        Assert.That(Statistics.LatentMixedFitStatuses, Is.Unique);
        Assert.That(Statistics.LatentMixedFitStatuses, Does.Contain(FeatureFitStatus.Separated));
        Assert.That(Statistics.LatentMixedFitStatuses, Does.Not.Contain(FeatureFitStatus.VarianceRatioAtLimit),
            "LogitNormalMixedModel never returns it (cc67b7a2)");
        Assert.That(Statistics.LatentMixedFitStatuses.Select(Statistics.StatusName), Is.All.Not.Null);
    }

    [Test]
    public void FitLatentMixed_AGivenResidualVariance_IsUsedForThatFeatureOnly()
    {
        string responses = WriteTemp("given.tsv",
            "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\n" +
            "A\t1.1\t2.3\t0.4\t1.9\t2.8\t3.1\t0.2\t1.0\n" +
            "B\t0.1\t1.3\t1.4\t0.9\t1.8\t2.1\t1.2\t2.0\n");
        string design = WriteTemp("given_design.tsv", EightSampleDesign);
        string groups = WriteTemp("given_groups.tsv", EightSampleGroups);
        string sigma = WriteTemp("given_sigma.tsv", "feature\tsigma2\nB\t0.25\n");
        JsonElement root = FitLatentMixed("x\n", "--responses", responses, "--design", design, "--groups", groups,
            "--link", "identity", "--residual-variance", sigma);
        Assert.That(Strings(root, "residual_variance_given"), Is.EqualTo(new[] { "False", "True" }));
        Assert.That(Doubles(root, "residual_variance")[1], Is.EqualTo(0.25).Within(1e-15));
        Assert.That(Doubles(root, "residual_variance")[0], Is.Not.EqualTo(0.25), "A, not listed, is estimated");
        Assert.That(Caveats(root), Has.Some.StartsWith("1 of 2 features were fitted at the residual_variance you gave"));
    }

    // ---- refusals --------------------------------------------------------------------------------

    [Test]
    public void FitLatentMixed_RefusesBadInput_AsUsage()
    {
        string responses = WriteTemp("r.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nA\t1\t2\t3\t\t5\t6\t7\t8\n");
        string counts = WriteTemp("c.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nA\t\t\t\t2/5\t\t\t\t\n");
        string design = WriteTemp("d.tsv", EightSampleDesign);
        string groups = WriteTemp("g.tsv", EightSampleGroups);
        string Refusal(params string[] options) =>
            Usage(Statistics.FitLatentMixed, "x\n", new[] { "stats", "fit-latent-mixed", "--design", design, "--groups", groups }
                .Concat(options).ToArray());
        string[] both = { "--responses", responses, "--counts", counts };

        Assert.That(Refusal("--link", "identity"), Does.Contain("Give --responses"));
        Assert.That(Refusal(both), Does.Contain("--link"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "1" }).ToArray()), Does.Contain("2 to 4 comma-separated numbers"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "0,one" }).ToArray()), Does.Contain("slope must be a number"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--estimator", "REML" }).ToArray()), Does.Contain("must be 'reml' or 'ml'"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--depth" }).ToArray()), Does.Contain("without a value"));
        Assert.That(Refusal("--counts", WriteTemp("bad.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nA\t3 of 5\t\t\t\t\t\t\t\n"), "--link", "identity"),
            Does.Contain("'3 of 5' is not a count"));
        Assert.That(Refusal("--responses", responses, "--counts", WriteTemp("other.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts9\nB\t\t\t\t\t\t\t\t\n"), "--link", "identity"),
            Does.Contain("Features not in the counts: A.").And.Contain("Features in the counts but not in the responses: B.")
                .And.Contain("Samples not in the counts: s8.").And.Contain("Samples in the counts but not in the responses: s9."));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--depth", WriteTemp("dep.tsv", "sample\tdepth\ns1\t1\n") }).ToArray()),
            Does.Contain("Not in the depth file: s2"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--residual-variance", WriteTemp("rv.tsv", "feature\tv\nZ\t1\n") }).ToArray()),
            Does.Contain("features the fit does not have: Z."));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--residual-variance", WriteTemp("rv2.tsv", "feature\tv\nA\tbig\n") }).ToArray()),
            Does.Contain("'big' is not a number"));

        // mzLib's refusals, told by the names the caller wrote.
        string measuredAndCounted = WriteTemp("mc.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nA\t\t2/5\t\t\t\t\t\t\n");
        Assert.That(Refusal("--responses", responses, "--counts", measuredAndCounted, "--link", "identity"),
            Does.StartWith("Feature 0, sample 1 is both measured and counted.").And.EndWith("That is feature 'A', sample 's2'."));
        string tooMany = WriteTemp("tm.tsv", "feature\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\nA\t\t\t\t6/5\t\t\t\t\n");
        Assert.That(Refusal("--responses", responses, "--counts", tooMany, "--link", "identity"),
            Does.Contain("6 of 5 is not a count").And.EndWith("That is feature 'A', sample 's4'."));
        Assert.That(Refusal(both.Concat(new[] { "--link", "0,0" }).ToArray()), Does.Contain("non-zero slope"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "0,1,0,1" }).ToArray()), Does.Contain("overdispersion in [0, 1)"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "0,1,0.3" }).ToArray()), Does.Contain("supply one depth per sample"));
        Assert.That(Refusal(both.Concat(new[] { "--link", "identity", "--residual-variance", WriteTemp("rv3.tsv", "feature\tv\nA\t-1\n") }).ToArray()),
            Does.Contain("positive and finite"));
    }

    // ---- fit-link --------------------------------------------------------------------------------

    /// <summary>
    /// mzLib's own FitLink test cells through the wire: the link is mzLib's direct call, exactly, and recovers
    /// the simulated truth at mzLib's tolerances (FitLinkRecoversASimulatedLink).
    /// </summary>
    [Test]
    public void FitLink_IsMzLibsFitLink_AndRecoversTheSimulatedLink()
    {
        string table = LinkCellsTable();
        string cells = WriteTemp("cells.tsv", table);
        JsonElement root = FitLinkVerb("--cells", cells, "--with-depth", "true");
        DetectionLink direct = LogitNormalMixedModel.FitLink(ParseLinkCells(table), withDepth: true)!;
        JsonElement link = root.GetProperty("link");
        Assert.That(link.GetProperty("intercept").GetDouble(), Is.EqualTo(direct.Intercept));
        Assert.That(link.GetProperty("slope").GetDouble(), Is.EqualTo(direct.Slope));
        Assert.That(link.GetProperty("depth_slope").GetDouble(), Is.EqualTo(direct.DepthSlope));
        Assert.That(link.GetProperty("overdispersion").GetDouble(), Is.EqualTo(0));
        Assert.That(direct.Intercept, Is.EqualTo(-0.5).Within(0.15));
        Assert.That(direct.Slope, Is.EqualTo(0.5).Within(0.1));
        Assert.That(direct.DepthSlope, Is.EqualTo(0.3).Within(0.1));
        Assert.That(root.GetProperty("cell_count").GetInt32(), Is.EqualTo(4000));
        Assert.That(root.GetProperty("depth_given").GetBoolean(), Is.True);
        Assert.That(root.GetProperty("with_depth").GetBoolean(), Is.True);
        Assert.That(root.GetProperty("with_overdispersion").GetBoolean(), Is.False, "absent: mzLib's declared default");
        Assert.That(root.GetProperty("detection_only").GetBoolean(), Is.True, "absent: mzLib's declared default");
        Assert.That(Caveats(root), Has.Some.StartsWith("detection_only:"));
        Assert.That(Caveats(root), Has.Some.StartsWith("overdispersion was not fitted"));
    }

    [Test]
    public void FitLink_PassesEachOption_AndEchoesMzLibsDeclaredDefaults()
    {
        Assert.That(Statistics.FitLinkDefault<bool>("withDepth"), Is.False, "FitLink's declared default at cc67b7a2");
        Assert.That(Statistics.FitLinkDefault<bool>("withOverdispersion"), Is.False);
        Assert.That(Statistics.FitLinkDefault<bool>("detectionOnly"), Is.True);

        string table = LinkCellsTable(count: 600, seed: 11);
        string cells = WriteTemp("cells_small.tsv", table);
        List<LinkCell> parsed = ParseLinkCells(table);
        JsonElement full = FitLinkVerb("--cells", cells, "--with-overdispersion", "true", "--detection-only", "FALSE");
        DetectionLink direct = LogitNormalMixedModel.FitLink(parsed, withDepth: false, withOverdispersion: true, detectionOnly: false)!;
        JsonElement link = full.GetProperty("link");
        Assert.That(link.GetProperty("slope").GetDouble(), Is.EqualTo(direct.Slope));
        Assert.That(link.GetProperty("overdispersion").GetDouble(), Is.EqualTo(direct.Overdispersion));
        Assert.That(link.GetProperty("depth_slope").GetDouble(), Is.EqualTo(0));
        Assert.That(full.GetProperty("detection_only").GetBoolean(), Is.False);
        Assert.That(Caveats(full), Has.Some.StartsWith("Fitted on the full counts"));
        Assert.That(Caveats(full), Has.Some.StartsWith("depth_slope was not fitted"));

        JsonElement plain = FitLinkVerb("--cells", cells);
        DetectionLink defaults = LogitNormalMixedModel.FitLink(parsed)!;
        Assert.That(plain.GetProperty("link").GetProperty("intercept").GetDouble(), Is.EqualTo(defaults.Intercept));
        Assert.That(plain.GetProperty("with_depth").GetBoolean(), Is.False);
    }

    [Test]
    public void FitLink_WithoutADepthColumn_UsesLinkCellsOwnDepth_AndSaysSo()
    {
        string table = LinkCellsTable(count: 300, seed: 3, withDepth: false);
        string cells = WriteTemp("cells_nodepth.tsv", table);
        JsonElement root = FitLinkVerb("--cells", cells, "--with-depth", "true");
        DetectionLink direct = LogitNormalMixedModel.FitLink(ParseLinkCells(table), withDepth: true)!;
        Assert.That(root.GetProperty("link").GetProperty("intercept").GetDouble(), Is.EqualTo(direct.Intercept));
        Assert.That(root.GetProperty("depth_given").GetBoolean(), Is.False);
        Assert.That(Caveats(root), Has.Some.StartsWith("The depth slope was fitted, but the cells file has no depth column"));
    }

    [Test]
    public void FitLink_RefusesBadInput_AsUsage()
    {
        string Refusal(string content, params string[] options) =>
            Usage(Statistics.FitLink, "", new[] { "stats", "fit-link", "--cells", WriteTemp("bad_cells.tsv", content) }.Concat(options).ToArray());
        const string Header = "latent_mean\tlatent_sd\tsuccesses\ttrials\n";
        Assert.That(Refusal("latent_mean\tlatent_sd\tsuccesses\ttrials\tdeapth\n0\t1\t1\t2\t0\n"),
            Does.Contain("Not a column this file takes: deapth."));
        Assert.That(Refusal("latent_mean\tsuccesses\ttrials\n0\t1\t2\n"), Does.Contain("Missing: latent_sd."));
        Assert.That(Refusal(Header), Does.Contain("no cells"));
        Assert.That(Refusal(Header + "zero\t1\t1\t2\n"), Does.Contain("latent_mean 'zero' is not a number"));
        Assert.That(Refusal(Header + "0\t1\t1.5\t2\n"), Does.Contain("successes '1.5' is not a whole number"));
        Assert.That(Refusal(Header + "0\t1\t1\t2\n", "--with-depth", "yes"), Does.Contain("must be 'true' or 'false'"));
        Assert.That(Refusal(Header + "0\t1\t1\t2\n", "--detection-only"), Does.Contain("without a value"));
        // mzLib's refusal: more successes than trials.
        Assert.That(Refusal(Header + "0\t1\t3\t2\n"), Does.Contain("0 ≤ successes ≤ trials"));
        Assert.That(Usage(Statistics.FitLink, "", "stats", "fit-link"), Does.Contain("--cells"));
    }
}
