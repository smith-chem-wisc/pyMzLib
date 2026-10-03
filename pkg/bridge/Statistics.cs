using System.Globalization;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// Differential-abundance statistics over a feature-by-sample table: per-feature linear models with
/// empirical-Bayes moderated t-tests, Benjamini–Hochberg adjustment, and random-effects
/// meta-analysis. Every number is mzLib's <c>StatisticalModels</c> (mzLib #1341, #1357); the bridge
/// only parses the tables and shapes the result.
/// </summary>
/// <remarks>
/// <para>
/// <c>stats fit</c> is ONE verb for the question a differential-abundance analysis asks — "which
/// features change with this coefficient, and how sure are we?" — rather than three verbs a caller
/// must chain. It runs <see cref="LinearModel.Fit"/> once and then
/// <see cref="EmpiricalBayes.Moderate"/> once per tested coefficient. Moderation already applies
/// <see cref="MultipleTesting.BenjaminiHochberg"/> across features, so the adjusted p-values come
/// back in the same table.
/// </para>
/// <para>
/// The variance prior depends only on the fit, never on which coefficient is tested, so it is
/// reported once. That holds because <c>Moderate</c> calls <see cref="EmpiricalBayes.FitPrior"/> on
/// the fit's residual variances, residual df and (with a trend) average responses alone.
/// </para>
/// </remarks>
internal static class Statistics
{
    /// <summary>Row columns of the <c>stats fit</c> table, in order.</summary>
    internal static readonly string[] FitColumns =
    {
        "feature", "coefficient", "status", "estimate", "standard_error", "t", "df_total", "p_value",
        "bh_adjusted", "posterior_variance", "prior_variance", "sigma", "df_residual", "observed",
        "average_response",
    };

    /// <summary>Row columns of the <c>stats meta</c> table, in order.</summary>
    internal static readonly string[] MetaColumns =
    {
        "feature", "studies", "estimate", "standard_error", "confidence_low", "confidence_high",
        "p_value", "tau2", "q", "i_squared", "direction_agree", "direction_disagree",
        "leave_one_out_max_delta",
    };

    /// <summary>
    /// <c>stats fit --responses FILE --design FILE [--trend] [--spline-basis N] [--threads N]</c>,
    /// with the coefficients to test on stdin, one name per line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>responses</b> is a TSV: a header row, then one row per feature. The first column is the
    /// feature id (unique); every other column is one sample, named by its header. Cells are used as
    /// given, so log-transform first. A blank cell, <c>NA</c> or <c>NaN</c> is missing and is omitted
    /// from that feature's fit, never imputed. <b>A 0 is an observation</b>, as in mzLib.
    /// </para>
    /// <para>
    /// <b>design</b> is a TSV: a header row, then one row per sample. The first column is the sample
    /// name, matching the responses header exactly; every other column is one coefficient, named by
    /// its header, and every cell must be a finite number. Rows are matched to response columns by
    /// name, so their order does not matter. Include an intercept column yourself if you want one.
    /// </para>
    /// </remarks>
    public static object Fit(Program.Arguments arguments)
    {
        string responsesPath = RequireFile(arguments, "responses");
        string designPath = RequireFile(arguments, "design");
        bool trend = arguments.Flag("trend");
        if (arguments.WasProvided("spline-basis") && !trend)
            throw new Program.UsageException(
                "Option --spline-basis only shapes the intensity trend; give --trend as well, or omit --spline-basis.");
        RequireValueIfProvided(arguments, "spline-basis");
        int splineBasis = arguments.OptionalInt("spline-basis", EmpiricalBayes.DefaultSplineBasisCount);
        if (splineBasis < 1)
            throw new Program.UsageException($"Option --spline-basis must be 1 or more; got {splineBasis}.");
        int threads = Proteins.ThreadsFrom(arguments);

        List<string> tested = Program.ReadStdinLines().Select(l => l.Trim()).Where(l => l.Length > 0).ToList();

        ResponseTable responses = ReadResponses(responsesPath);
        DesignTable design = ReadDesign(designPath, responses.SampleNames);

        if (tested.Count == 0)
            throw new Program.UsageException(
                "No coefficient to test was given on stdin. Give one name per line; the design has: " +
                string.Join(", ", design.CoefficientNames) + ".");
        string? duplicate = tested.GroupBy(c => c, StringComparer.Ordinal).FirstOrDefault(g => g.Count() > 1)?.Key;
        if (duplicate != null)
            throw new Program.UsageException($"Coefficient '{duplicate}' was given more than once on stdin.");
        foreach (string coefficient in tested)
        {
            if (!design.CoefficientNames.Contains(coefficient, StringComparer.Ordinal))
                throw new Program.UsageException(
                    $"No coefficient named '{coefficient}' in the design. The design has: " +
                    string.Join(", ", design.CoefficientNames) + ".");
        }

        LinearModelFit fit;
        var tests = new List<ModeratedTest>(tested.Count);
        try
        {
            fit = LinearModel.Fit(responses.Values, design.Values, design.CoefficientNames, threads);
            foreach (string coefficient in tested)
                tests.Add(EmpiricalBayes.Moderate(fit, coefficient, trend, splineBasis));
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            // mzLib's own validation, reclassified as the caller's input: a design that is not of full
            // rank, or a table too sparse to estimate a variance prior from, is a question asked wrongly,
            // not a fault, and must exit 2 like every other malformed argument.
            throw new Program.UsageException(invalid.Message.Split(" (Parameter")[0]);
        }

        VariancePrior prior = tests[0].Prior;
        bool dfDiffer = tests[0].ResidualDfDiffer;

        var table = new Sdrf.ColumnTable(FitColumns);
        foreach (ModeratedTest test in tests)
        {
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                table.Add(
                    responses.FeatureIds[f],
                    test.Coefficient,
                    StatusName(fit.Status[f]),
                    Finite(test.Estimate[f]),
                    Finite(test.StandardError[f]),
                    Finite(test.T[f]),
                    Finite(test.DfTotal[f]),
                    Finite(test.PValue[f]),
                    Finite(test.BenjaminiHochbergAdjusted[f]),
                    Finite(test.PosteriorVariance[f]),
                    Finite(prior.Scale[f]),
                    Finite(fit.Sigma[f]),
                    fit.Status[f] == FeatureFitStatus.Fitted ? fit.DfResidual[f] : null,
                    fit.Observed[f],
                    Finite(fit.AverageResponse[f]));
            }
        }

        var statusCounts = Enum.GetValues<FeatureFitStatus>().ToDictionary(
            StatusName, s => fit.Status.Count(x => x == s));

        return new
        {
            responses_file = Path.GetFullPath(responsesPath),
            design_file = Path.GetFullPath(designPath),
            feature_count = fit.FeatureCount,
            sample_count = responses.SampleNames.Count,
            sample_names = responses.SampleNames,
            coefficient_names = design.CoefficientNames,
            tested_coefficients = tested,
            trend,
            threads,
            missing_count = responses.MissingCount,
            zero_count = responses.ZeroCount,
            status_counts = statusCounts,
            residual_df_differ = dfDiffer,
            prior = new
            {
                df = Finite(prior.Df),
                df_infinite = double.IsPositiveInfinity(prior.Df),
                trended = prior.Trended,
                spline_basis_count = prior.SplineBasisCount,
                // One value without a trend; with a trend it varies per feature and is the
                // prior_variance column instead.
                scale = prior.Trended ? null : FirstFinite(prior.Scale),
            },
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = FitCaveats(fit, dfDiffer, responses.ZeroCount, prior),
        };
    }

    /// <summary>
    /// <c>stats adjust</c>: Benjamini–Hochberg adjustment of the p-values on stdin, one per line.
    /// </summary>
    /// <remarks>
    /// Blank lines are kept so output row i is input line i (as <c>sdrf parse-age</c> does); a blank,
    /// <c>NA</c> or <c>NaN</c> line is untested, stays null, and is not counted in m. A value outside
    /// [0, 1] is refused by mzLib and reported as a usage error.
    /// </remarks>
    public static object Adjust(Program.Arguments arguments)
    {
        List<string> lines = Sdrf.ReadStdinCells();
        if (lines.Count == 0)
            throw new Program.UsageException("No p-values were given on stdin. Give one per line.");

        var pValues = new double[lines.Count];
        for (int i = 0; i < lines.Count; i++)
        {
            pValues[i] = ParseCell(lines[i], out bool missing);
            if (!missing && double.IsNaN(pValues[i]))
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin is not a number: '{lines[i].Trim()}'. Leave a line blank, or write NA, for a feature that was not tested.");
        }

        double[] adjusted;
        try
        {
            adjusted = MultipleTesting.BenjaminiHochberg(pValues);
        }
        catch (ArgumentOutOfRangeException outOfRange)
        {
            throw new Program.UsageException(outOfRange.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0]);
        }

        var table = new Sdrf.ColumnTable("p_value", "bh_adjusted");
        for (int i = 0; i < pValues.Length; i++)
            table.Add(Finite(pValues[i]), Finite(adjusted[i]));

        return new
        {
            line_count = lines.Count,
            tested_count = pValues.Count(double.IsFinite),
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = new[]
            {
                "Benjamini-Hochberg adjusts a family of p-values for the false discovery rate among the features tested. It is not a target-decoy q-value for identifications.",
                "m, the family size, is the number of finite p-values: a blank or NA line is not tested and does not count.",
            },
        };
    }

    /// <summary>
    /// <c>stats meta [--confidence 0.95]</c>: random-effects pooling per feature, the studies on
    /// stdin as <c>feature\testimate\tstandard_error</c>, one study per line.
    /// </summary>
    /// <remarks>
    /// No header line. Lines of one feature need not be adjacent; features are reported in order of
    /// first appearance. Each feature is pooled by <see cref="RandomEffectsMeta.Pool"/> independently.
    /// </remarks>
    public static object Meta(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "confidence");
        double confidence = arguments.OptionalDouble("confidence", 0.95);
        if (!(confidence > 0 && confidence < 1))
            throw new Program.UsageException($"Option --confidence must be between 0 and 1, exclusive; got {confidence.ToString(CultureInfo.InvariantCulture)}.");

        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No studies were given on stdin. Give one per line: feature<TAB>estimate<TAB>standard_error.");

        var order = new List<string>();
        var byFeature = new Dictionary<string, (List<double> Estimates, List<double> Errors)>(StringComparer.Ordinal);
        for (int i = 0; i < lines.Count; i++)
        {
            string[] parts = lines[i].Split('\t');
            if (parts.Length != 3)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin has {parts.Length} tab-separated fields; expected 3: feature, estimate, standard_error.");
            string feature = parts[0].Trim();
            if (feature.Length == 0)
                throw new Program.UsageException($"Line {i + 1} of stdin has an empty feature name.");
            double estimate = ParseRequired(parts[1], $"Line {i + 1} estimate");
            double error = ParseRequired(parts[2], $"Line {i + 1} standard_error");
            if (!byFeature.TryGetValue(feature, out var studies))
            {
                studies = (new List<double>(), new List<double>());
                byFeature[feature] = studies;
                order.Add(feature);
            }
            studies.Estimates.Add(estimate);
            studies.Errors.Add(error);
        }

        var table = new Sdrf.ColumnTable(MetaColumns);
        foreach (string feature in order)
        {
            var (estimates, errors) = byFeature[feature];
            MetaAnalysisResult pooled;
            try
            {
                pooled = RandomEffectsMeta.Pool(estimates, errors, confidence);
            }
            catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
            {
                throw new Program.UsageException($"Feature '{feature}': {invalid.Message.Split(" (Parameter")[0]}");
            }
            table.Add(feature, pooled.Studies, Finite(pooled.Estimate), Finite(pooled.StandardError),
                Finite(pooled.ConfidenceLow), Finite(pooled.ConfidenceHigh), Finite(pooled.PValue),
                Finite(pooled.Tau2), Finite(pooled.Q), Finite(pooled.ISquared), pooled.DirectionAgree,
                pooled.DirectionDisagree, Finite(pooled.LeaveOneOutMaxDelta));
        }

        return new
        {
            study_count = lines.Count,
            feature_count = order.Count,
            confidence,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = new[]
            {
                "DerSimonian-Laird random effects with a normal reference: the interval and p-value do not widen for few studies (no Hartung-Knapp adjustment).",
                "A feature with one study is reported with tau2 = 0 and q, i_squared and leave_one_out_max_delta null: one study has no heterogeneity to measure.",
                "Each feature is pooled independently; p_value is not adjusted across features. Pass the p_value column to stats adjust for that.",
            },
        };
    }

    // ---- input tables ----------------------------------------------------------------------------

    private sealed record ResponseTable(List<string> FeatureIds, List<string> SampleNames, double[,] Values,
        int MissingCount, int ZeroCount);

    private sealed record DesignTable(List<string> CoefficientNames, double[,] Values);

    private static ResponseTable ReadResponses(string path)
    {
        List<string[]> rows = ReadTsv(path, "responses");
        string[] header = rows[0];
        if (header.Length < 2)
            throw new Program.UsageException(
                $"The responses file '{path}' has {header.Length} column; it needs a feature id column and at least one sample column.");
        List<string> samples = header.Skip(1).Select(h => h.Trim()).ToList();
        RequireUnique(samples, "sample column", "responses header");

        int features = rows.Count - 1;
        if (features == 0)
            throw new Program.UsageException($"The responses file '{path}' has a header but no feature rows.");
        var ids = new List<string>(features);
        var values = new double[features, samples.Count];
        int missing = 0, zeros = 0;
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the responses file has {row.Length} fields; the header has {header.Length}.");
            ids.Add(row[0].Trim());
            for (int s = 0; s < samples.Count; s++)
            {
                double value = ParseCell(row[s + 1], out bool isMissing);
                if (!isMissing && double.IsNaN(value))
                    throw new Program.UsageException(
                        $"Line {r + 1} of the responses file, sample '{samples[s]}': '{row[s + 1].Trim()}' is not a number. Leave it blank, or write NA, for a missing value.");
                if (isMissing) missing++;
                else if (value == 0) zeros++;
                values[r - 1, s] = value;
            }
        }
        RequireUnique(ids, "feature id", "responses file's first column");
        return new ResponseTable(ids, samples, values, missing, zeros);
    }

    private static DesignTable ReadDesign(string path, List<string> sampleNames)
    {
        List<string[]> rows = ReadTsv(path, "design");
        string[] header = rows[0];
        if (header.Length < 2)
            throw new Program.UsageException(
                $"The design file '{path}' has {header.Length} column; it needs a sample name column and at least one coefficient column.");
        List<string> coefficients = header.Skip(1).Select(h => h.Trim()).ToList();
        RequireUnique(coefficients, "coefficient", "design header");

        var bySample = new Dictionary<string, double[]>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the design file has {row.Length} fields; the header has {header.Length}.");
            string sample = row[0].Trim();
            if (bySample.ContainsKey(sample))
                throw new Program.UsageException($"Sample '{sample}' appears more than once in the design file.");
            var cells = new double[coefficients.Count];
            for (int c = 0; c < coefficients.Count; c++)
                cells[c] = ParseRequired(row[c + 1], $"Design sample '{sample}', coefficient '{coefficients[c]}'");
            bySample[sample] = cells;
        }

        List<string> absent = sampleNames.Where(s => !bySample.ContainsKey(s)).ToList();
        List<string> extra = bySample.Keys.Where(s => !sampleNames.Contains(s, StringComparer.Ordinal)).ToList();
        if (absent.Count > 0 || extra.Count > 0)
            throw new Program.UsageException(
                "The design's samples must be exactly the responses' sample columns." +
                (absent.Count > 0 ? " Not in the design: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? " In the design but not in the responses: " + string.Join(", ", extra) + "." : ""));

        var values = new double[sampleNames.Count, coefficients.Count];
        for (int s = 0; s < sampleNames.Count; s++)
            for (int c = 0; c < coefficients.Count; c++)
                values[s, c] = bySample[sampleNames[s]][c];
        return new DesignTable(coefficients, values);
    }

    private static List<string[]> ReadTsv(string path, string what)
    {
        var rows = File.ReadAllLines(path)
            .Select((line, i) => i == 0 ? line.TrimStart('﻿') : line)
            .Where(line => line.Trim().Length > 0)
            .Select(line => line.TrimEnd('\r').Split('\t'))
            .ToList();
        if (rows.Count == 0)
            throw new Program.UsageException($"The {what} file '{path}' is empty.");
        return rows;
    }

    private static void RequireUnique(List<string> names, string what, string where)
    {
        if (names.Any(n => n.Length == 0))
            throw new Program.UsageException($"An empty {what} name in the {where}.");
        string? duplicate = names.GroupBy(n => n, StringComparer.Ordinal).FirstOrDefault(g => g.Count() > 1)?.Key;
        if (duplicate != null)
            throw new Program.UsageException($"The {what} '{duplicate}' appears more than once in the {where}.");
    }

    /// <summary>A numeric cell. Blank, NA and NaN are missing (NaN, <paramref name="missing"/> true);
    /// anything else unparseable is NaN with <paramref name="missing"/> false, for the caller to refuse.</summary>
    private static double ParseCell(string raw, out bool missing)
    {
        string cell = raw.Trim();
        missing = cell.Length == 0
            || cell.Equals("NA", StringComparison.OrdinalIgnoreCase)
            || cell.Equals("NaN", StringComparison.OrdinalIgnoreCase);
        if (missing)
            return double.NaN;
        return double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
            ? value
            : double.NaN;
    }

    private static double ParseRequired(string raw, string what)
    {
        double value = ParseCell(raw, out _);
        if (!double.IsFinite(value))
            throw new Program.UsageException($"{what} must be a finite number; got '{raw.Trim()}'.");
        return value;
    }

    private static string RequireFile(Program.Arguments arguments, string name)
    {
        string path = arguments.Required(name);
        if (!File.Exists(path))
            throw new Program.UsageException($"The --{name} file was not found: '{path}'.");
        return path;
    }

    private static void RequireValueIfProvided(Program.Arguments arguments, string name)
    {
        if (arguments.Flag(name))
            throw new Program.UsageException($"Option --{name} was given without a value.");
    }

    // ---- output ----------------------------------------------------------------------------------

    /// <summary>mzLib's <see cref="FeatureFitStatus"/> as the wire writes it.</summary>
    internal static string StatusName(FeatureFitStatus status) => status switch
    {
        FeatureFitStatus.Fitted => "fitted",
        FeatureFitStatus.TooFewObservations => "too_few_observations",
        FeatureFitStatus.RankDeficient => "rank_deficient",
        _ => status.ToString(),
    };

    private static List<string> FitCaveats(LinearModelFit fit, bool dfDiffer, int zeroCount, VariancePrior prior)
    {
        var caveats = new List<string>
        {
            "Moderation is limma's legacy empirical-Bayes estimator (eBayes(legacy = TRUE), Smyth 2004); mzLib reproduces it to 1e-8 relative. Not implemented: robust = TRUE, contrasts, the B-statistic, observation weights.",
            "bh_adjusted is Benjamini-Hochberg over the fitted features, per coefficient. It is not a target-decoy q-value.",
            "Missing cells are omitted per feature, never imputed. A 0 is an observation: log-transform first, or blank out unmeasured values.",
        };
        if (dfDiffer)
            caveats.Add(
                "residual_df_differ: the fitted features do not all have the same residual degrees of freedom (missing values were omitted). " +
                "limma 3.61 and later then default to a different prior estimator (eBayes legacy = FALSE), so default limma would report a different prior and different moderated statistics for this input.");
        int notFitted = fit.Status.Count(s => s != FeatureFitStatus.Fitted);
        if (notFitted > 0)
        {
            int tooFew = fit.Status.Count(s => s == FeatureFitStatus.TooFewObservations);
            int rank = fit.Status.Count(s => s == FeatureFitStatus.RankDeficient);
            caveats.Add(
                $"{notFitted} of {fit.FeatureCount} features were not fitted ({tooFew} too_few_observations, {rank} rank_deficient); their statistics are null and they are not in the BH family.");
        }
        if (zeroCount > 0)
            caveats.Add(
                $"{zeroCount} response cells are exactly 0 and were fitted as observations. If 0 means 'not measured' in your table, blank those cells.");
        if (double.IsPositiveInfinity(prior.Df))
            caveats.Add(
                "The prior df is infinite: the residual variances are no more dispersed than sampling alone explains, so every feature takes the prior variance.");
        return caveats;
    }

    private static double? Finite(double value) => double.IsFinite(value) ? value : null;

    private static double? FirstFinite(IReadOnlyList<double> values)
    {
        foreach (double value in values)
            if (double.IsFinite(value))
                return value;
        return null;
    }
}
