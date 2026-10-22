using System.Globalization;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats fit-mixed</c>: one random-intercept linear mixed model per feature of a feature-by-sample
/// table, by mzLib's <see cref="MixedModel"/> (mzLib #1367, 1.0.594). The bridge reads the tables as
/// <c>stats fit</c> does, plus one group label per sample, and makes one mzLib call; every number is
/// mzLib's.
/// </summary>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats fit-mixed</c> table, in order.</summary>
    internal static readonly string[] MixedColumns =
    {
        "feature", "coefficient", "status", "estimate", "standard_error", "df", "t", "p_value",
        "residual_variance", "group_variance", "log_likelihood", "observed", "groups",
    };

    /// <summary>
    /// The outcomes <see cref="MixedModel.Fit"/> can give a feature, and so the keys of
    /// <c>status_counts</c>. Read from the method body at mzLib d10017ce (FitOne): too few observed
    /// samples, a rank-deficient observed design, fewer than two observed groups or none observed
    /// twice, the variance ratio at the top of its search, a failed final evaluation, or fitted. mzLib
    /// does not yet say which of <see cref="FeatureFitStatus"/>'s values each model returns
    /// (mzLib#1487 asks it to); until it does, this list is the bridge's reading of the code, held by
    /// a test.
    /// </summary>
    internal static readonly FeatureFitStatus[] MixedFitStatuses =
    {
        FeatureFitStatus.Fitted, FeatureFitStatus.TooFewObservations, FeatureFitStatus.RankDeficient,
        FeatureFitStatus.TooFewGroups, FeatureFitStatus.VarianceRatioAtLimit, FeatureFitStatus.NotConverged,
    };

    /// <summary>
    /// <c>stats fit-mixed --responses FILE --design FILE --groups FILE [--estimator reml|ml]
    /// [--threads N]</c>, with the coefficients to report on stdin, one name per line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>responses</b> and <b>design</b> are <c>stats fit</c>'s. <b>groups</b> is a TSV with a header
    /// row and one row per sample: the sample name, then exactly one column holding the sample's
    /// group label (a dataset, a donor, a day), whose header names the grouping. It must list exactly
    /// the responses' samples, in any order. Labels are compared as text, exactly.
    /// </para>
    /// <para>
    /// <b>estimator</b> is sent only when given. Absent, mzLib's default applies
    /// (<see cref="MixedModel.Fit"/>'s <c>estimator</c> parameter), and the envelope echoes the
    /// estimator mzLib reports it used (<see cref="MixedModelFit.Estimator"/>).
    /// </para>
    /// </remarks>
    public static object FitMixed(Program.Arguments arguments)
    {
        string responsesPath = RequireFile(arguments, "responses");
        string designPath = RequireFile(arguments, "design");
        string groupsPath = RequireFile(arguments, "groups");
        RequireValueIfProvided(arguments, "estimator");
        VarianceEstimator? estimator = arguments.Optional("estimator") switch
        {
            null => null,
            "reml" => VarianceEstimator.Reml,
            "ml" => VarianceEstimator.MaximumLikelihood,
            var other => throw new Program.UsageException(
                $"Option --estimator must be 'reml' or 'ml' (mzLib's VarianceEstimator Reml or MaximumLikelihood); got '{other}'."),
        };
        int threads = Proteins.ThreadsFrom(arguments);

        List<string> tested = ReadTestedCoefficients();
        ResponseTable responses = ReadResponses(responsesPath);
        DesignTable design = ReadDesign(designPath, responses.SampleNames);
        GroupTable groups = ReadGroups(groupsPath, responses.SampleNames);
        CheckTestedCoefficients(tested, design);

        MixedModelFit fit;
        try
        {
            // The estimator is passed only when the caller named one, so an absent one is mzLib's default.
            fit = estimator is { } chosen
                ? MixedModel.Fit(responses.Values, design.Values, groups.Labels, design.CoefficientNames, chosen, threads)
                : MixedModel.Fit(responses.Values, design.Values, groups.Labels, design.CoefficientNames, maxThreads: threads);
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            // A design that is not of full rank, or a sample with an empty group label: the caller's
            // input, refused by mzLib, so a usage error like every malformed argument.
            throw new Program.UsageException(MzLibRefusal(invalid, responses));
        }

        var table = new Sdrf.ColumnTable(MixedColumns);
        int untestable = 0;
        foreach (string coefficient in tested)
        {
            int j = fit.IndexOf(coefficient);
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                double df = fit.DegreesOfFreedom(f, j);
                if (fit.Status[f] == FeatureFitStatus.Fitted && !double.IsFinite(df))
                    untestable++;
                table.Add(
                    responses.FeatureIds[f],
                    coefficient,
                    StatusName(fit.Status[f]),
                    Finite(fit.Coefficient(f, j)),
                    Finite(fit.StandardError(f, j)),
                    Finite(df),
                    Finite(fit.TValue(f, j)),
                    Finite(fit.PValue(f, j)),
                    Finite(fit.ResidualVariance[f]),
                    Finite(fit.GroupVariance[f]),
                    Finite(fit.LogLikelihood[f]),
                    fit.Observed[f],
                    fit.Groups[f]);
            }
        }

        var statusCounts = MixedFitStatuses.ToDictionary(StatusName, s => fit.Status.Count(x => x == s));
        List<string> groupNames = groups.Labels.Distinct(StringComparer.Ordinal).ToList();

        return new
        {
            responses_file = Path.GetFullPath(responsesPath),
            design_file = Path.GetFullPath(designPath),
            groups_file = Path.GetFullPath(groupsPath),
            grouping = groups.Name,
            group_names = groupNames,
            group_count = groupNames.Count,
            feature_count = fit.FeatureCount,
            sample_count = responses.SampleNames.Count,
            sample_names = responses.SampleNames,
            coefficient_names = design.CoefficientNames,
            tested_coefficients = tested,
            estimator = EstimatorName(fit.Estimator),
            threads,
            missing_count = responses.MissingCount,
            zero_count = responses.ZeroCount,
            status_counts = statusCounts,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = MixedCaveats(fit, untestable, responses.ZeroCount),
        };
    }

    /// <summary>
    /// mzLib's <see cref="VarianceEstimator"/> as the wire writes it, the same words --estimator
    /// takes. Every value is named here, so a value mzLib adds fails a test on the repin instead of
    /// reaching the wire under its C# name.
    /// </summary>
    internal static string EstimatorName(VarianceEstimator estimator) => estimator switch
    {
        VarianceEstimator.Reml => "reml",
        VarianceEstimator.MaximumLikelihood => "ml",
        _ => throw new InvalidOperationException(
            $"VarianceEstimator.{estimator} has no wire name. Name it in Statistics.EstimatorName."),
    };

    /// <summary>One group label per sample, in the responses header's order, and the grouping's name.</summary>
    private sealed record GroupTable(string Name, List<string> Labels);

    /// <summary>
    /// A TSV of one row per sample: the sample name, then exactly one column, the group label, whose
    /// header names the grouping. Rows are matched to the responses' samples by name, which they must
    /// be exactly. Cells are trimmed. An empty label is passed on for mzLib to refuse.
    /// </summary>
    private static GroupTable ReadGroups(string path, List<string> sampleNames)
    {
        List<string[]> rows = ReadTsv(path, "groups");
        string[] header = rows[0];
        if (header.Length != 2)
            throw new Program.UsageException(
                $"The groups file '{path}' has {header.Length} column{(header.Length == 1 ? "" : "s")}; " +
                "it needs exactly two: the sample name, then its group label.");
        string name = header[1].Trim();
        if (name.Length == 0)
            throw new Program.UsageException($"The groups file '{path}' has an empty header for its group column.");

        var bySample = new Dictionary<string, string>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the groups file has {row.Length} fields; the header has {header.Length}.");
            string sample = row[0].Trim();
            if (bySample.ContainsKey(sample))
                throw new Program.UsageException($"Sample '{sample}' appears more than once in the groups file.");
            bySample[sample] = row[1].Trim();
        }

        List<string> absent = sampleNames.Where(s => !bySample.ContainsKey(s)).ToList();
        List<string> extra = bySample.Keys.Where(s => !sampleNames.Contains(s, StringComparer.Ordinal)).ToList();
        if (absent.Count > 0 || extra.Count > 0)
            throw new Program.UsageException(
                "The groups file's samples must be exactly the responses' sample columns." +
                (absent.Count > 0 ? " Not in the groups file: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? " In the groups file but not in the responses: " + string.Join(", ", extra) + "." : ""));

        return new GroupTable(name, sampleNames.Select(s => bySample[s]).ToList());
    }

    private static List<string> MixedCaveats(MixedModelFit fit, int untestable, int zeroCount)
    {
        bool reml = fit.Estimator == VarianceEstimator.Reml;
        var caveats = new List<string>
        {
            "A random intercept per group and nothing more: y = X beta + u_group + e, with u ~ N(0, group_variance) and e ~ N(0, residual_variance), as nlme::lme(y ~ ..., random = ~ 1 | group). No random slopes, nested or crossed groupings, or residual correlation.",
            "Fixed effects are tested by Wald t on nlme's containment degrees of freedom: observed samples minus observed groups minus the within-group coefficients, for a coefficient that varies within some group and for the intercept; observed groups minus the group-level coefficients, minus 1 more when the design has an intercept, for one constant within every group. Satterthwaite and Kenward-Roger (lmerTest, pbkrtest) are not implemented, so p-values differ from lmerTest's. mzLib matches nlme to 1e-5 relative, nlme's own precision, with identical df.",
            "group_variance 0 is the boundary estimate: the groups' baselines are no more different than the residual variance explains. nlme works on the log scale, cannot reach 0, and reports a tiny positive value there instead.",
            reml
                ? "reml: the variances are restricted maximum-likelihood estimates, and log_likelihood is the restricted log-likelihood. Compare reml log-likelihoods only between models with the same fixed effects; use --estimator ml to compare designs."
                : "ml: residual_variance is the maximum-likelihood estimate RSS / observed, but standard errors use RSS / (observed - coefficients), as nlme's summary reports ML fits. log_likelihood is the full log-likelihood, comparable between designs fitted to the same samples.",
            "p_value is not adjusted across features. Send one coefficient's p_value column to stats adjust for Benjamini-Hochberg.",
            "Missing cells are omitted per feature, never imputed. A 0 is an observation: log-transform first, or blank out unmeasured values.",
            "A feature the fixed effects fit exactly (a constant one, say) is reported fitted, with a residual_variance at rounding level (about 1e-16) and a meaningless t and p_value: mzLib's not_converged guard catches only a residual sum of squares of exactly 0. Check residual_variance before trusting a very small p_value.",
        };
        int notFitted = fit.Status.Count(s => s != FeatureFitStatus.Fitted);
        if (notFitted > 0)
        {
            string counts = string.Join(", ", MixedFitStatuses.Where(s => s != FeatureFitStatus.Fitted)
                .Select(s => (Status: s, Count: fit.Status.Count(x => x == s)))
                .Where(c => c.Count > 0)
                .Select(c => $"{c.Count} {StatusName(c.Status)}"));
            caveats.Add($"{notFitted} of {fit.FeatureCount} features were not fitted ({counts}); their statistics are null.");
        }
        if (untestable > 0)
            caveats.Add(
                $"{untestable} fitted rows have no positive containment df: a group-level coefficient with too few observed groups to test it. Their df and p_value are null; estimate, standard_error and t are reported.");
        int boundary = Enumerable.Range(0, fit.FeatureCount)
            .Count(f => fit.Status[f] == FeatureFitStatus.Fitted && fit.GroupVariance[f] == 0);
        if (boundary > 0)
            caveats.Add(
                $"{boundary} of the fitted features have group_variance 0 (the boundary): for them the model is ordinary least squares, tested on the containment df.");
        if (zeroCount > 0)
            caveats.Add(
                $"{zeroCount} response cells are exactly 0 and were fitted as observations. If 0 means 'not measured' in your table, blank those cells.");
        return caveats;
    }
}
