using System.Globalization;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats fit</c>'s contrasts and prior estimator: <c>--contrasts FILE [--confidence LEVEL]</c>, tested by
/// mzLib's <see cref="EmpiricalBayes.ModerateContrasts"/> (mzLib #1452, 1.0.595), and
/// <c>--estimator</c>, mzLib's <see cref="VariancePriorEstimator"/> (#1368, 1.0.594). The bridge reads
/// the contrasts table and shapes the result; every number is mzLib's.
/// </summary>
/// <remarks>
/// <para>
/// <b>The contrasts file is what <c>diff design --out</c> writes</b> (<c>contrasts_i.tsv</c>): a header of
/// the contrast id column, then the design's coefficient columns, then one row per contrast, its id and one
/// weight per coefficient. Its columns are matched to the design's by name, so their order does not matter;
/// a column the design lacks, or a design column the file lacks, is refused, because a weight cannot be
/// guessed. The weights are passed to mzLib in the design's coefficient order.
/// </para>
/// <para>
/// <b>Defaults are mzLib's.</b> Without <c>--estimator</c> or <c>--confidence</c> the bridge calls mzLib's
/// overload without that argument, so mzLib's default applies, and the envelope echoes what mzLib reports it
/// used (<see cref="VariancePrior.Estimator"/>, <see cref="ModeratedContrast.ConfidenceLevel"/>). A range
/// rule mzLib owns (the confidence level strictly between 0 and 1, a contrast not all zero, names unique)
/// is left to mzLib, whose refusal is reported as the caller's mistake.
/// </para>
/// <para>
/// <b>A separate table, not more rows.</b> The contrasts come back in their own <c>contrasts</c> object,
/// with its own columns, so the coefficient table is exactly what it was without <c>--contrasts</c>. A
/// contrast row carries an interval a coefficient row does not, and is keyed by a contrast id rather than a
/// coefficient name; one table holding both would leave half of every row's keys null.
/// </para>
/// </remarks>
internal static partial class Statistics
{
    /// <summary>Row columns of <c>stats fit</c>'s <c>contrasts</c> table, in order.</summary>
    internal static readonly string[] ContrastColumns =
    {
        "feature", "contrast", "status", "estimate", "standard_error", "confidence_low", "confidence_high",
        "t", "df_total", "p_value", "bh_adjusted", "posterior_variance", "prior_variance", "sigma",
        "df_residual", "observed", "average_response",
    };

    /// <summary>
    /// mzLib's <see cref="VariancePriorEstimator"/> as the wire writes it, the same words <c>--estimator</c>
    /// takes. Every value is named here, so a value mzLib adds fails a test on the repin instead of reaching
    /// the wire under its C# name.
    /// </summary>
    internal static string PriorEstimatorName(VariancePriorEstimator estimator) => estimator switch
    {
        VariancePriorEstimator.MomentsLegacy => "moments_legacy",
        VariancePriorEstimator.MarginalLikelihood => "marginal_likelihood",
        _ => throw new InvalidOperationException(
            $"VariancePriorEstimator.{estimator} has no wire name. Name it in Statistics.PriorEstimatorName."),
    };

    /// <summary>The estimator <c>--estimator</c> names, or null when it was not given (mzLib's default then).</summary>
    private static VariancePriorEstimator? PriorEstimatorFrom(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "estimator");
        string? text = arguments.Optional("estimator");
        if (text is null)
            return null;
        VariancePriorEstimator[] all = Enum.GetValues<VariancePriorEstimator>();
        foreach (VariancePriorEstimator estimator in all)
        {
            if (PriorEstimatorName(estimator) == text)
                return estimator;
        }
        throw new Program.UsageException(
            $"Option --estimator must be one of {string.Join(", ", all.Select(PriorEstimatorName))} " +
            $"(mzLib's VariancePriorEstimator {string.Join(" or ", all)}); got '{text}'.");
    }

    /// <summary>
    /// <c>--confidence</c>, or null when it was not given (mzLib's default then). Only a number is checked
    /// here; the range is mzLib's rule, refused by <see cref="EmpiricalBayes.ModerateContrasts"/>.
    /// </summary>
    private static double? ConfidenceFrom(Program.Arguments arguments, bool contrastsGiven)
    {
        if (!arguments.WasProvided("confidence"))
            return null;
        if (!contrastsGiven)
            throw new Program.UsageException(
                "Option --confidence sets the contrasts' interval level; give --contrasts as well, or omit --confidence.");
        return arguments.OptionalDouble("confidence", double.NaN);
    }

    /// <summary>
    /// The contrasts file: a header (the contrast id column, any name, then exactly the design's coefficient
    /// columns, in any order) and one row per contrast, its id and a finite weight per column. Returned in
    /// file order, each contrast's weights in the design's coefficient order.
    /// </summary>
    private static List<ContrastWeights> ReadContrasts(string path, DesignTable design)
    {
        List<string[]> rows = ReadTsv(path, "contrasts");
        string[] header = rows[0];
        List<string> columns = header.Skip(1).Select(h => h.Trim()).ToList();
        RequireUnique(columns, "coefficient", "contrasts header");

        List<string> absent = design.CoefficientNames.Where(c => !columns.Contains(c, StringComparer.Ordinal)).ToList();
        List<string> extra = columns.Where(c => !design.CoefficientNames.Contains(c, StringComparer.Ordinal)).ToList();
        if (absent.Count > 0 || extra.Count > 0)
            throw new Program.UsageException(
                "The contrasts file's columns after the contrast id must be exactly the design's coefficients (" +
                string.Join(", ", design.CoefficientNames) + "), one weight each." +
                (absent.Count > 0 ? " Not in the contrasts file: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? " In the contrasts file but not in the design: " + string.Join(", ", extra) + "." : ""));
        if (rows.Count == 1)
            throw new Program.UsageException($"The contrasts file '{path}' has a header but no contrast rows.");

        int[] position = design.CoefficientNames.Select(c => columns.IndexOf(c) + 1).ToArray();
        var contrasts = new List<ContrastWeights>(rows.Count - 1);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the contrasts file has {row.Length} fields; the header has {header.Length}.");
            string name = row[0].Trim();
            if (name.Length == 0)
                throw new Program.UsageException($"Line {r + 1} of the contrasts file has an empty contrast id.");
            var weights = new double[position.Length];
            for (int j = 0; j < position.Length; j++)
                weights[j] = ParseRequired(row[position[j]], $"Contrast '{name}', coefficient '{design.CoefficientNames[j]}'");
            contrasts.Add(new ContrastWeights(name, weights));
        }
        return contrasts;
    }

    /// <summary>One coefficient's moderated test, with the estimator only when the caller named one.</summary>
    private static ModeratedTest ModerateAsAsked(LinearModelFit fit, string coefficient, bool trend, int? splineBasis,
        VariancePriorEstimator? estimator) =>
        estimator is { } chosen
            ? EmpiricalBayes.Moderate(fit, coefficient, trend, splineBasis, chosen)
            : EmpiricalBayes.Moderate(fit, coefficient, trend, splineBasis);

    /// <summary>
    /// The contrasts' moderated tests, passing the estimator and the confidence level only when the caller
    /// gave them, so an absent one is mzLib's default.
    /// </summary>
    private static IReadOnlyList<ModeratedContrast> ModerateContrastsAsAsked(LinearModelFit fit,
        IReadOnlyList<ContrastWeights> contrasts, bool trend, int? splineBasis, VariancePriorEstimator? estimator,
        double? confidence) =>
        (estimator, confidence) switch
        {
            ({ } e, { } level) => EmpiricalBayes.ModerateContrasts(fit, contrasts, trend, splineBasis, e, level),
            ({ } e, null) => EmpiricalBayes.ModerateContrasts(fit, contrasts, trend, splineBasis, e),
            (null, { } level) => EmpiricalBayes.ModerateContrasts(fit, contrasts, trend, splineBasis, confidenceLevel: level),
            (null, null) => EmpiricalBayes.ModerateContrasts(fit, contrasts, trend, splineBasis),
        };

    /// <summary>The <c>contrasts</c> object: one row per (contrast, feature), contrasts in file order.</summary>
    private static object ContrastResult(IReadOnlyList<ModeratedContrast> contrasts, LinearModelFit fit,
        ResponseTable responses, VariancePrior prior)
    {
        var table = new Sdrf.ColumnTable(ContrastColumns);
        foreach (ModeratedContrast contrast in contrasts)
        {
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                table.Add(
                    responses.FeatureIds[f],
                    contrast.Name,
                    StatusName(fit.Status[f]),
                    Finite(contrast.Estimate[f]),
                    Finite(contrast.StandardError[f]),
                    Finite(contrast.ConfidenceLow[f]),
                    Finite(contrast.ConfidenceHigh[f]),
                    Finite(contrast.T[f]),
                    Finite(contrast.DfTotal[f]),
                    Finite(contrast.PValue[f]),
                    Finite(contrast.BenjaminiHochbergAdjusted[f]),
                    Finite(contrast.PosteriorVariance[f]),
                    Finite(prior.Scale[f]),
                    Finite(fit.Sigma[f]),
                    fit.Status[f] == FeatureFitStatus.Fitted ? fit.DfResidual[f] : null,
                    fit.Observed[f],
                    Finite(fit.AverageResponse[f]));
            }
        }

        return new
        {
            confidence_level = contrasts[0].ConfidenceLevel,
            tested_contrasts = contrasts.Select(c => c.Name).ToList(),
            weights = contrasts.ToDictionary(c => c.Name, c => c.Weights.ToList(), StringComparer.Ordinal),
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
        };
    }

    /// <summary>The caveats that hold whenever contrasts were tested.</summary>
    private static IEnumerable<string> ContrastCaveats(IReadOnlyList<ModeratedContrast> contrasts)
    {
        string level = contrasts[0].ConfidenceLevel.ToString("R", CultureInfo.InvariantCulture);
        yield return
            "Contrasts are tested with the same variance prior as the coefficients (EmpiricalBayes.ModerateContrasts, as limma's contrasts.fit, eBayes and topTable(confint)), each with its own Benjamini-Hochberg family over its fitted features. " +
            $"confidence_low and confidence_high are a two-sided {level} interval per feature, estimate +/- t-quantile((1 + {level}) / 2, df_total) x standard_error, not adjusted across features.";
        yield return
            "With missing values a contrast's standard error is exact per feature, sqrt(c'(X'X)^-1 c) on the samples that feature was observed in; limma's contrasts.fit approximates it there, so on such data the contrasts differ from limma's by design. With complete data they agree, and a contrast of one coefficient gives exactly that coefficient's row.";
    }
}
