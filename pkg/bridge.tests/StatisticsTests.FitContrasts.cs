using System.Globalization;
using System.Text;
using System.Text.Json;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>stats fit --contrasts</c>, <c>--confidence</c> and <c>--estimator</c>. The contrasts are held through the
/// wire to mzLib's frozen limma 3.68.5 reference (lmFit, contrasts.fit, eBayes(legacy = TRUE),
/// topTable(confint = 0.95), and the exact per-feature reference for missing values), and every default to
/// what mzLib itself reports when the bridge sends nothing.
/// </summary>
public partial class StatisticsTests
{
    // ---- inputs ----------------------------------------------------------------------------------

    /// <summary>
    /// mzLib's contrast reference as the verb's three TSVs. The responses carry feature ids and sample columns
    /// 1..12 already; the design gets those sample names, its rows written REVERSED; and the weights file's
    /// coefficient columns are written in REVERSE order, so the test also proves the bridge matches a
    /// contrasts file's columns to the design by name, not by position.
    /// </summary>
    private (string Responses, string Design, string Contrasts) WriteContrastReference(string data)
    {
        string responses = Reference($"limma_contrast_responses_{data}.tsv");
        string[] samples = File.ReadLines(responses).First().Split('\t').Skip(1).ToArray();

        var (designHeader, designRows) = ReadReference("limma_contrast_design.tsv");
        var design = new StringBuilder("sample\t" + string.Join('\t', designHeader) + "\n");
        for (int s = designRows.Count - 1; s >= 0; s--)
            design.Append(samples[s]).Append('\t').Append(string.Join('\t', designRows[s])).Append('\n');

        var (weightHeader, weightRows) = ReadReference("limma_contrast_weights.tsv");
        int[] order = new[] { 0 }.Concat(Enumerable.Range(1, weightHeader.Length - 1).Reverse()).ToArray();
        var contrasts = new StringBuilder(string.Join('\t', order.Select(j => weightHeader[j])) + "\n");
        foreach (string[] row in weightRows)
            contrasts.Append(string.Join('\t', order.Select(j => row[j]))).Append('\n');

        string designPath = Path.Combine(_temp, "limma_contrast_design_with_samples.tsv");
        string contrastsPath = Path.Combine(_temp, "limma_contrast_weights_reordered.tsv");
        File.WriteAllText(designPath, design.ToString());
        File.WriteAllText(contrastsPath, contrasts.ToString());
        return (responses, designPath, contrastsPath);
    }

    /// <summary>The MALAT1 dilution series with a contrasts file on its committed design's columns.</summary>
    private (string Responses, string Design, string Contrasts) WriteMalatContrasts(string body)
    {
        string path = Path.Combine(_temp, "malat_contrasts.tsv");
        File.WriteAllText(path, body);
        return (Path.Combine(StatsFixtures, "malat_dilution_log2_vs_fluc.tsv"),
            Path.Combine(StatsFixtures, "malat_dilution_design.tsv"), path);
    }

    private const string MalatContrasts =
        "contrast\tintercept\tmalat_250ng\tmalat_125ng\n" +
        "half\t0\t1\t0\n" +
        "quarter_vs_half\t0\t-1\t1\n";

    /// <summary>One column of the contrasts table, optionally one contrast's block, null as NaN.</summary>
    private static double[] ContrastDoubles(JsonElement root, string column, string? contrast = null)
    {
        JsonElement table = root.GetProperty("contrasts").GetProperty("columns");
        var values = table.GetProperty(column).EnumerateArray().ToList();
        var names = table.GetProperty("contrast").EnumerateArray().Select(e => e.GetString()).ToList();
        return Enumerable.Range(0, values.Count)
            .Where(i => contrast is null || names[i] == contrast)
            .Select(i => values[i].ValueKind == JsonValueKind.Null ? double.NaN : values[i].GetDouble())
            .ToArray();
    }

    /// <summary>mzLib called directly on the same inputs, for the numbers the wire must carry.</summary>
    private static (LinearModelFit Fit, List<ContrastWeights> Weights) DirectFit(string responsesPath, string designPath,
        string contrastsPath)
    {
        List<string[]> table = File.ReadAllLines(responsesPath).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        List<string[]> design = File.ReadAllLines(designPath).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        List<string[]> weights = File.ReadAllLines(contrastsPath).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        string[] samples = table[0].Skip(1).ToArray();
        string[] coefficients = design[0].Skip(1).ToArray();
        var bySample = design.Skip(1).ToDictionary(r => r[0], r => r.Skip(1).Select(Parse).ToArray());

        var y = new double[table.Count - 1, samples.Length];
        for (int f = 0; f < table.Count - 1; f++)
            for (int s = 0; s < samples.Length; s++)
            {
                string cell = table[f + 1][s + 1];
                y[f, s] = cell.Length == 0 || cell == "NaN" ? double.NaN : Parse(cell);
            }
        var x = new double[samples.Length, coefficients.Length];
        for (int s = 0; s < samples.Length; s++)
            for (int j = 0; j < coefficients.Length; j++)
                x[s, j] = bySample[samples[s]][j];

        int[] column = coefficients.Select(c => Array.IndexOf(weights[0], c)).ToArray();
        var contrasts = weights.Skip(1)
            .Select(r => new ContrastWeights(r[0], column.Select(k => Parse(r[k])).ToArray())).ToList();
        return (LinearModel.Fit(y, x, coefficients, 1), contrasts);
    }

    // ---- agreement with limma ----------------------------------------------------------------------

    [TestCase("complete", false)]
    [TestCase("complete", true)]
    [TestCase("na", false)]
    [TestCase("na", true)]
    public void FitContrasts_ReproduceLimmasContrastsFit_ThroughTheWire(string data, bool trend)
    {
        var (responses, design, contrasts) = WriteContrastReference(data);
        var args = new List<string> { "stats", "fit", "--responses", responses, "--design", design, "--contrasts", contrasts };
        if (trend) args.Add("--trend");

        JsonElement root = Run(Statistics.Fit, "", args.ToArray());

        JsonElement block = root.GetProperty("contrasts");
        Assert.That(block.GetProperty("tested_contrasts").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { "C_vs_B", "age", "BC_vs_A" }), "file order");
        Assert.That(block.GetProperty("weights").GetProperty("C_vs_B").EnumerateArray().Select(e => e.GetDouble()),
            Is.EqualTo(new[] { 0.0, -1, 1, 0 }), "echoed in the design's coefficient order, whatever the file's column order");
        Assert.That(block.GetProperty("row_count").GetInt32(), Is.EqualTo(3 * 300));
        Assert.That(block.GetProperty("column_names").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(Statistics.ContrastColumns));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(0), "no coefficient on stdin, so no coefficient rows");
        Assert.That(root.GetProperty("tested_coefficients").GetArrayLength(), Is.EqualTo(0));
        Assert.That(root.GetProperty("contrasts_file").GetString(), Is.EqualTo(Path.GetFullPath(contrasts)));

        string file = $"limma_contrast_{data}_{(trend ? "trend" : "none")}.tsv";
        var (header, rows) = ReadReference(file);
        double[] Expected(string contrast, string column) =>
            rows.Where(r => r[0] == contrast).Select(r => Parse(r[Array.IndexOf(header, column)])).ToArray();
        foreach (string contrast in new[] { "C_vs_B", "age", "BC_vs_A" })
        {
            string where = $"{file} {contrast} ";
            Assert.That(rows.Where(r => r[0] == contrast).Select(r => r[1]),
                Is.EqualTo(ReadReference($"limma_contrast_responses_{data}.tsv").Rows.Select(r => r[0])), where + "feature order");
            AssertClose(ContrastDoubles(root, "estimate", contrast), Expected(contrast, "estimate"), where + "estimate");
            AssertClose(ContrastDoubles(root, "standard_error", contrast), Expected(contrast, "se"), where + "se");
            AssertClose(ContrastDoubles(root, "t", contrast), Expected(contrast, "t"), where + "t");
            AssertClose(ContrastDoubles(root, "df_total", contrast), Expected(contrast, "df_total"), where + "df.total");
            AssertClose(ContrastDoubles(root, "p_value", contrast), Expected(contrast, "p"), where + "p");
            AssertClose(ContrastDoubles(root, "bh_adjusted", contrast), Expected(contrast, "adj_p"), where + "adj.p");
            AssertClose(ContrastDoubles(root, "confidence_low", contrast), Expected(contrast, "ci_low"), where + "CI.L");
            AssertClose(ContrastDoubles(root, "confidence_high", contrast), Expected(contrast, "ci_high"), where + "CI.R");
            AssertClose(ContrastDoubles(root, "posterior_variance", contrast), Expected(contrast, "s2_post"), where + "s2.post");
            AssertClose(ContrastDoubles(root, "prior_variance", contrast), Expected(contrast, "s2_prior"), where + "s2.prior");
        }
        AssertClose(new[] { root.GetProperty("prior").GetProperty("df").GetDouble() },
            new[] { Expected("C_vs_B", "df_prior")[0] }, "df.prior");
        Assert.That(root.GetProperty("prior").GetProperty("estimator").GetString(), Is.EqualTo("moments_legacy"));
    }

    // ---- defaults are mzLib's -----------------------------------------------------------------------

    [Test]
    public void FitContrasts_WithoutConfidence_TakeTheLevelMzLibReports()
    {
        var (responses, design, contrasts) = WriteMalatContrasts(MalatContrasts);
        JsonElement root = Run(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts);

        var (fit, weights) = DirectFit(responses, design, contrasts);
        IReadOnlyList<ModeratedContrast> direct = EmpiricalBayes.ModerateContrasts(fit, weights, trend: false);
        Assert.That(root.GetProperty("contrasts").GetProperty("confidence_level").GetDouble(),
            Is.EqualTo(direct[0].ConfidenceLevel), "the bridge sends no level, so the level is the one mzLib reports");
        for (int c = 0; c < direct.Count; c++)
        {
            Assert.That(ContrastDoubles(root, "confidence_low", direct[c].Name), Is.EqualTo(direct[c].ConfidenceLow), direct[c].Name);
            Assert.That(ContrastDoubles(root, "confidence_high", direct[c].Name), Is.EqualTo(direct[c].ConfidenceHigh), direct[c].Name);
            Assert.That(ContrastDoubles(root, "t", direct[c].Name), Is.EqualTo(direct[c].T), direct[c].Name);
        }
        List<string> caveats = root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();
        Assert.That(caveats.Any(x => x.StartsWith("Contrasts are tested with the same variance prior", StringComparison.Ordinal)
            && x.Contains($"two-sided {direct[0].ConfidenceLevel.ToString("R", CultureInfo.InvariantCulture)} interval", StringComparison.Ordinal)));
        Assert.That(caveats.Any(x => x.Contains("limma's contrasts.fit approximates it", StringComparison.Ordinal)));
    }

    [Test]
    public void FitContrasts_AtAGivenConfidence_AreMzLibsIntervalsAtThatLevel()
    {
        var (responses, design, contrasts) = WriteMalatContrasts(MalatContrasts);
        JsonElement root = Run(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts, "--confidence", "0.8");

        var (fit, weights) = DirectFit(responses, design, contrasts);
        IReadOnlyList<ModeratedContrast> direct = EmpiricalBayes.ModerateContrasts(fit, weights, trend: false, confidenceLevel: 0.8);
        Assert.That(root.GetProperty("contrasts").GetProperty("confidence_level").GetDouble(), Is.EqualTo(0.8));
        foreach (ModeratedContrast contrast in direct)
        {
            Assert.That(ContrastDoubles(root, "confidence_low", contrast.Name), Is.EqualTo(contrast.ConfidenceLow), contrast.Name);
            Assert.That(ContrastDoubles(root, "confidence_high", contrast.Name), Is.EqualTo(contrast.ConfidenceHigh), contrast.Name);
        }

        JsonElement wider = Run(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts);
        double[] narrow = ContrastDoubles(root, "confidence_high", "half").Zip(ContrastDoubles(root, "confidence_low", "half"), (h, l) => h - l).ToArray();
        double[] wide = ContrastDoubles(wider, "confidence_high", "half").Zip(ContrastDoubles(wider, "confidence_low", "half"), (h, l) => h - l).ToArray();
        Assert.That(Enumerable.Range(0, narrow.Length).Where(i => double.IsFinite(narrow[i])).All(i => narrow[i] < wide[i]),
            "an 80% interval is narrower than the default one");
    }

    [Test]
    public void Fit_WithoutEstimator_UsesMzLibsDefault_AndEchoesIt()
    {
        var (responses, design) = WriteLimmaReference();
        JsonElement absent = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design);
        JsonElement legacy = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design,
            "--estimator", "moments_legacy");

        var (y, x, names) = LimmaReferenceArrays();
        VariancePriorEstimator mzLibDefault = EmpiricalBayes.Moderate(LinearModel.Fit(y, x, names, 1), "age_decades", trend: false).Prior.Estimator;
        Assert.That(absent.GetProperty("prior").GetProperty("estimator").GetString(),
            Is.EqualTo(Statistics.PriorEstimatorName(mzLibDefault)));
        Assert.That(Doubles(absent, "t"), Is.EqualTo(Doubles(legacy, "t")));
        Assert.That(absent.GetProperty("contrasts_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(absent.GetProperty("contrasts").ValueKind, Is.EqualTo(JsonValueKind.Null), "no --contrasts, no contrasts table");
        Assert.That(absent.GetProperty("column_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(Statistics.FitColumns));
    }

    [Test]
    public void Fit_WithMarginalLikelihood_IsMzLibsEstimator_ForCoefficientsAndContrasts()
    {
        var (responses, design) = WriteLimmaReference();
        string contrasts = Path.Combine(_temp, "limma_reference_contrasts.tsv");
        File.WriteAllText(contrasts, "contrast\tIntercept\tage_decades\tsex\nage\t0\t1\t0\nage_minus_sex\t0\t1\t-1\n");
        JsonElement root = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts, "--estimator", "marginal_likelihood");

        var (y, x, names) = LimmaReferenceArrays();
        LinearModelFit fit = LinearModel.Fit(y, x, names, 1);
        ModeratedTest direct = EmpiricalBayes.Moderate(fit, "age_decades", trend: false, estimator: VariancePriorEstimator.MarginalLikelihood);
        ModeratedTest legacy = EmpiricalBayes.Moderate(fit, "age_decades", trend: false);
        var weights = new List<ContrastWeights>
        {
            new("age", new[] { 0.0, 1, 0 }), new("age_minus_sex", new[] { 0.0, 1, -1 }),
        };
        IReadOnlyList<ModeratedContrast> directContrasts =
            EmpiricalBayes.ModerateContrasts(fit, weights, trend: false, estimator: VariancePriorEstimator.MarginalLikelihood);

        Assert.That(root.GetProperty("prior").GetProperty("estimator").GetString(), Is.EqualTo("marginal_likelihood"));
        Assert.That(root.GetProperty("prior").GetProperty("df").GetDouble(), Is.EqualTo(direct.Prior.Df));
        Assert.That(root.GetProperty("prior").GetProperty("df").GetDouble(), Is.Not.EqualTo(legacy.Prior.Df), "a different estimator");
        Assert.That(Doubles(root, "t"), Is.EqualTo(direct.T.ToArray()));
        foreach (ModeratedContrast contrast in directContrasts)
        {
            Assert.That(ContrastDoubles(root, "t", contrast.Name), Is.EqualTo(contrast.T), contrast.Name);
            Assert.That(ContrastDoubles(root, "confidence_low", contrast.Name), Is.EqualTo(contrast.ConfidenceLow), contrast.Name);
        }
        Assert.That(ContrastDoubles(root, "t", "age"), Is.EqualTo(Doubles(root, "t")),
            "a one-coefficient contrast gives exactly that coefficient's row, under either estimator");

        List<string> caveats = root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();
        Assert.That(caveats[0], Does.StartWith("marginal_likelihood:"));
        Assert.That(caveats.Single(c => c.StartsWith("residual_df_differ:", StringComparison.Ordinal)),
            Does.Contain("marginal_likelihood handles that exactly"));
    }

    [Test]
    public void Fit_ResidualDfDiffer_UnderTheDefault_PointsToTheEstimatorThatHandlesIt()
    {
        var (responses, design) = WriteLimmaReference();
        JsonElement root = Run(Statistics.Fit, "age_decades\n", "stats", "fit", "--responses", responses, "--design", design);
        List<string> caveats = root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();
        Assert.That(caveats[0], Does.StartWith("moments_legacy:"));
        Assert.That(caveats.Single(c => c.StartsWith("residual_df_differ:", StringComparison.Ordinal)),
            Does.Contain("--estimator marginal_likelihood"));
    }

    [Test]
    public void EveryVariancePriorEstimator_HasAWireName()
    {
        string[] names = Enum.GetValues<VariancePriorEstimator>().Select(Statistics.PriorEstimatorName).ToArray();
        Assert.That(names, Is.EqualTo(new[] { "moments_legacy", "marginal_likelihood" }));
    }

    // ---- the coefficient table does not move ---------------------------------------------------------

    [Test]
    public void FitContrasts_LeaveTheCoefficientTable_AsItIsWithoutThem()
    {
        var (responses, design, contrasts) = WriteMalatContrasts(MalatContrasts);
        JsonElement plain = Run(Statistics.Fit, "malat_250ng\nmalat_125ng\n", "stats", "fit",
            "--responses", responses, "--design", design);
        JsonElement both = Run(Statistics.Fit, "malat_250ng\nmalat_125ng\n", "stats", "fit",
            "--responses", responses, "--design", design, "--contrasts", contrasts);

        Assert.That(both.GetProperty("columns").GetRawText(), Is.EqualTo(plain.GetProperty("columns").GetRawText()));
        Assert.That(both.GetProperty("prior").GetRawText(), Is.EqualTo(plain.GetProperty("prior").GetRawText()),
            "the prior depends on the fit, not on what is tested");
        Assert.That(ContrastDoubles(both, "t", "half"), Is.EqualTo(Doubles(both, "t", "malat_250ng")),
            "a one-coefficient contrast is exactly that coefficient's test");
        Assert.That(ContrastDoubles(both, "bh_adjusted", "half"), Is.EqualTo(Doubles(both, "bh_adjusted", "malat_250ng")));
        double[] difference = Doubles(both, "estimate", "malat_125ng").Zip(Doubles(both, "estimate", "malat_250ng"), (q, h) => q - h).ToArray();
        double[] contrast = ContrastDoubles(both, "estimate", "quarter_vs_half");
        for (int f = 0; f < contrast.Length; f++)
        {
            if (double.IsNaN(difference[f]))
                Assert.That(contrast[f], Is.NaN, $"feature {f}");
            else
                Assert.That(contrast[f], Is.EqualTo(difference[f]).Within(1e-12), $"feature {f}");
        }
    }

    // ---- refusals ------------------------------------------------------------------------------------

    [TestCase("contrast\tintercept\tmalat_250ng\nhalf\t0\t1\n", "Not in the contrasts file: malat_125ng.")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\theight\nhalf\t0\t1\t0\t0\n", "In the contrasts file but not in the design: height.")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_250ng\nhalf\t0\t1\t0\n", "appears more than once in the contrasts header")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\n", "has a header but no contrast rows")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\n\t0\t1\t0\n", "has an empty contrast id")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\nhalf\t0\t1\n", "has 3 fields; the header has 4")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\nhalf\t0\tone\t0\n", "Contrast 'half', coefficient 'malat_250ng' must be a finite number")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\nhalf\t0\tNA\t0\n", "must be a finite number")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\nnothing\t0\t0\t0\n", "Every weight of contrast 'nothing' is zero.")]
    [TestCase("contrast\tintercept\tmalat_250ng\tmalat_125ng\nhalf\t0\t1\t0\nhalf\t0\t0\t1\n", "Two contrasts are named 'half'.")]
    public void FitContrasts_RefuseAMalformedContrastsFile(string body, string message)
    {
        var (responses, design, contrasts) = WriteMalatContrasts(body);
        string text = Usage(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design, "--contrasts", contrasts);
        Assert.That(text, Does.Contain(message));
        Assert.That(text, Does.Not.Contain("(Parameter"));
    }

    [TestCase(new[] { "--confidence", "0.9" }, "give --contrasts as well")]
    [TestCase(new[] { "--contrasts" }, "Option --contrasts was given without a value.")]
    [TestCase(new[] { "--contrasts", "no_such_file.tsv" }, "The --contrasts file was not found")]
    [TestCase(new[] { "--estimator", "reml" }, "Option --estimator must be one of moments_legacy, marginal_likelihood")]
    [TestCase(new[] { "--estimator" }, "Option --estimator was given without a value.")]
    public void Fit_RefusesMalformedContrastAndEstimatorOptions(string[] extra, string message)
    {
        var (responses, design, _) = WriteMalatContrasts(MalatContrasts);
        string text = Usage(Statistics.Fit, "malat_250ng\n",
            new[] { "stats", "fit", "--responses", responses, "--design", design }.Concat(extra).ToArray());
        Assert.That(text, Does.Contain(message));
    }

    [TestCase("1.5", "The confidence level must be strictly between 0 and 1.")]
    [TestCase("0", "The confidence level must be strictly between 0 and 1.")]
    [TestCase("NaN", "The confidence level must be strictly between 0 and 1.")]
    [TestCase("ninety", "Option --confidence must be a number; got 'ninety'.")]
    public void FitContrasts_RefuseAConfidenceMzLibRefuses(string level, string message)
    {
        var (responses, design, contrasts) = WriteMalatContrasts(MalatContrasts);
        string text = Usage(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts, "--confidence", level);
        Assert.That(text, Is.EqualTo(message), "mzLib's own rule and words, without .NET's parameter tail");
    }

    [Test]
    public void FitContrasts_ConfidenceWithoutAValue_IsRefused()
    {
        var (responses, design, contrasts) = WriteMalatContrasts(MalatContrasts);
        string text = Usage(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design,
            "--contrasts", contrasts, "--confidence");
        Assert.That(text, Is.EqualTo("Option --confidence was given without a value."));
    }

    [Test]
    public void Fit_NeedsACoefficientOrAContrastsFile()
    {
        var (responses, design, _) = WriteMalatContrasts(MalatContrasts);
        string text = Usage(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", design);
        Assert.That(text, Does.StartWith("No coefficient to test was given on stdin."));
    }

    // ---- helpers ---------------------------------------------------------------------------------------

    /// <summary>mzLib's limma reference input as arrays, in its own sample order.</summary>
    private static (double[,] Y, double[,] X, string[] Names) LimmaReferenceArrays()
    {
        var (_, responseRows) = ReadReference("limma_responses.tsv");
        var (designHeader, designRows) = ReadReference("limma_design.tsv");
        var y = new double[responseRows.Count, responseRows[0].Length];
        for (int f = 0; f < responseRows.Count; f++)
            for (int s = 0; s < responseRows[f].Length; s++)
                y[f, s] = Parse(responseRows[f][s]);
        var x = new double[designRows.Count, designHeader.Length];
        for (int s = 0; s < designRows.Count; s++)
            for (int j = 0; j < designHeader.Length; j++)
                x[s, j] = Parse(designRows[s][j]);
        return (y, x, designHeader);
    }
}
