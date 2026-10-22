using System.Globalization;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats fit-logistic</c>: one logistic regression per feature of a 0/1 feature-by-sample table,
/// by mzLib's <see cref="LogisticRegression"/> (mzLib #1367, 1.0.594). The bridge reads the tables
/// as <c>stats fit</c> does and makes one mzLib call; every number is mzLib's.
/// </summary>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats fit-logistic</c> table, in order.</summary>
    internal static readonly string[] LogisticColumns =
    {
        "feature", "coefficient", "status", "estimate", "standard_error", "z", "p_value", "deviance",
        "iterations", "observed", "events",
    };

    /// <summary>
    /// The outcomes <see cref="LogisticRegression.Fit"/> can give a feature, and so the keys of
    /// <c>status_counts</c>. Read from the method body at mzLib d10017ce (FitOne): too few observed
    /// samples, a rank-deficient observed design, separation (including an all-0 or all-1 feature),
    /// no convergence within the iteration limit, or fitted. mzLib does not yet say which of
    /// <see cref="FeatureFitStatus"/>'s values each model returns (mzLib#1487 asks it to); until it
    /// does, this list is the bridge's reading of the code, held by a test.
    /// </summary>
    internal static readonly FeatureFitStatus[] LogisticFitStatuses =
    {
        FeatureFitStatus.Fitted, FeatureFitStatus.TooFewObservations, FeatureFitStatus.RankDeficient,
        FeatureFitStatus.Separated, FeatureFitStatus.NotConverged,
    };

    /// <summary>
    /// <c>stats fit-logistic --responses FILE --design FILE [--max-iterations N] [--tolerance X]
    /// [--threads N]</c>, with the coefficients to report on stdin, one name per line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>responses</b> is <c>stats fit</c>'s table, but every cell is an outcome: 0, 1, or blank /
    /// NA / NaN for missing (omitted from that feature's fit). mzLib refuses any other value, since a
    /// proportion is not a binary outcome. <b>design</b> is <c>stats fit</c>'s design.
    /// </para>
    /// <para>
    /// <b>max-iterations</b> and <b>tolerance</b> are sent only when given. Absent, mzLib's own
    /// defaults apply (<see cref="LogisticRegression.Fit"/>'s optional parameters), and the envelope
    /// echoes null for them: the bridge does not restate mzLib's numbers.
    /// </para>
    /// </remarks>
    public static object FitLogistic(Program.Arguments arguments)
    {
        string responsesPath = RequireFile(arguments, "responses");
        string designPath = RequireFile(arguments, "design");
        RequireValueIfProvided(arguments, "max-iterations");
        RequireValueIfProvided(arguments, "tolerance");
        int? maxIterations = arguments.WasProvided("max-iterations") ? arguments.OptionalInt("max-iterations", 0) : null;
        double? tolerance = arguments.WasProvided("tolerance") ? arguments.OptionalDouble("tolerance", 0) : null;
        int threads = Proteins.ThreadsFrom(arguments);

        List<string> tested = ReadTestedCoefficients();
        ResponseTable responses = ReadResponses(responsesPath);
        DesignTable design = ReadDesign(designPath, responses.SampleNames);
        CheckTestedCoefficients(tested, design);

        LogisticModelFit fit;
        try
        {
            // Each optional mzLib parameter is passed only when the caller gave it, so an absent one
            // is mzLib's default and never a number the bridge chose.
            fit = (maxIterations, tolerance) switch
            {
                (null, null) => LogisticRegression.Fit(responses.Values, design.Values, design.CoefficientNames, threads),
                ({ } limit, null) => LogisticRegression.Fit(responses.Values, design.Values, design.CoefficientNames, threads,
                    maxIterations: limit),
                (null, { } rule) => LogisticRegression.Fit(responses.Values, design.Values, design.CoefficientNames, threads,
                    tolerance: rule),
                ({ } limit, { } rule) => LogisticRegression.Fit(responses.Values, design.Values, design.CoefficientNames, threads,
                    limit, rule),
            };
        }
        catch (ArgumentOutOfRangeException outOfRange)
        {
            // mzLib's own range rules on the two IRLS controls (maxIterations >= 1, tolerance > 0).
            string option = outOfRange.ParamName == "maxIterations" ? "max-iterations" : "tolerance";
            string given = option == "max-iterations"
                ? maxIterations?.ToString(CultureInfo.InvariantCulture) ?? ""
                : tolerance?.ToString("R", CultureInfo.InvariantCulture) ?? "";
            throw new Program.UsageException(
                $"mzLib's LogisticRegression.Fit refuses --{option} {given}: {MzLibRefusal(outOfRange, responses)}");
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            // A response that is not 0, 1 or missing, or a design that is not of full rank: the
            // caller's input, refused by mzLib, so a usage error like every malformed argument.
            throw new Program.UsageException(MzLibRefusal(invalid, responses));
        }

        var table = new Sdrf.ColumnTable(LogisticColumns);
        foreach (string coefficient in tested)
        {
            int j = fit.IndexOf(coefficient);
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                table.Add(
                    responses.FeatureIds[f],
                    coefficient,
                    StatusName(fit.Status[f]),
                    Finite(fit.Coefficient(f, j)),
                    Finite(fit.StandardError(f, j)),
                    Finite(fit.ZValue(f, j)),
                    Finite(fit.PValue(f, j)),
                    Finite(fit.Deviance[f]),
                    fit.Iterations[f],
                    fit.Observed[f],
                    fit.Events[f]);
            }
        }

        var statusCounts = LogisticFitStatuses.ToDictionary(StatusName, s => fit.Status.Count(x => x == s));

        return new
        {
            responses_file = Path.GetFullPath(responsesPath),
            design_file = Path.GetFullPath(designPath),
            feature_count = fit.FeatureCount,
            sample_count = responses.SampleNames.Count,
            sample_names = responses.SampleNames,
            coefficient_names = design.CoefficientNames,
            tested_coefficients = tested,
            max_iterations = maxIterations,
            tolerance,
            threads,
            missing_count = responses.MissingCount,
            status_counts = statusCounts,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = LogisticCaveats(fit),
        };
    }

    private static List<string> LogisticCaveats(LogisticModelFit fit)
    {
        string threshold = LogisticRegression.SeparationThreshold.ToString("R", CultureInfo.InvariantCulture);
        var caveats = new List<string>
        {
            "Maximum likelihood by iteratively reweighted least squares, started and stopped by R glm's rules (relative change in deviance below the tolerance); standard errors from the inverse Fisher information at the estimate, and a two-sided Wald z test against the standard normal. mzLib matches R's glm(family = binomial, epsilon = 1e-14) to 1e-8 relative in coefficients and deviance, and to 1e-6 relative in standard errors and p-values, because summary.glm takes its weights from the step before the estimate. Wald tests are the least reliable of the likelihood tests near separation and in small samples; no likelihood-ratio test is offered.",
            $"separated: the design splits the feature's 0s from its 1s, so the maximum-likelihood estimate does not exist (a coefficient runs to infinity). mzLib reports the feature separated, with null statistics, when a fitted probability comes within {threshold} of 0 or 1 (LogisticRegression.SeparationThreshold), and also when every observed outcome is the same. R's glm instead prints large, meaningless coefficients and a warning. A separated feature can be the strongest effect in the table (detected in every run of one condition and in none of another): read its events and observed counts. Firth's penalized likelihood (logistf), which gives finite estimates there, is not implemented.",
            "p_value is not adjusted across features. Send one coefficient's p_value column to stats adjust for Benjamini-Hochberg.",
            "Missing cells are omitted per feature, never imputed. A 0 is an outcome (for a detection table: measured as absent), not a missing value.",
        };
        int notFitted = fit.Status.Count(s => s != FeatureFitStatus.Fitted);
        if (notFitted > 0)
        {
            string counts = string.Join(", ", LogisticFitStatuses.Where(s => s != FeatureFitStatus.Fitted)
                .Select(s => (Status: s, Count: fit.Status.Count(x => x == s)))
                .Where(c => c.Count > 0)
                .Select(c => $"{c.Count} {StatusName(c.Status)}"));
            caveats.Add($"{notFitted} of {fit.FeatureCount} features were not fitted ({counts}); their estimate, standard_error, z, p_value and deviance are null.");
        }
        return caveats;
    }
}
