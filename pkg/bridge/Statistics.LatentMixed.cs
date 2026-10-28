using System.Globalization;
using System.Text.RegularExpressions;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats fit-latent-mixed</c> and <c>stats fit-link</c>: mzLib's <see cref="LogitNormalMixedModel"/>
/// (mzLib #1469, 1.0.596), a random-intercept model on a latent scale whose samples are either measured
/// on it or counted through it, and the detection link the counts are read through. The bridge reads
/// the tables, makes one mzLib call each and shapes the result; every number is mzLib's.
/// </summary>
/// <remarks>
/// <para>
/// <b>Why a verb of its own rather than an option of <c>stats fit-mixed</c>.</b> mzLib has two entry
/// points, <see cref="MixedModel.Fit"/> and <see cref="LogitNormalMixedModel.Fit"/>, with different
/// optimisers, different statuses (no <c>variance_ratio_at_limit</c>; <c>separated</c>), different
/// standard errors under ML, and inputs the other lacks (counts, a link, a depth, a given residual
/// variance). One verb would have the bridge pick between them from which options were present. As
/// two verbs the caller names the model, and <c>stats fit-mixed</c> is not touched.
/// </para>
/// </remarks>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats fit-latent-mixed</c> table, in order.</summary>
    internal static readonly string[] LatentMixedColumns =
    {
        "feature", "coefficient", "status", "estimate", "standard_error", "df", "t", "p_value",
        "residual_variance", "residual_variance_given", "group_variance", "log_likelihood",
        "measured", "counted", "groups",
    };

    /// <summary>
    /// The outcomes <see cref="LogitNormalMixedModel.Fit"/> can give a feature, and so the keys of
    /// <c>status_counts</c>, in <see cref="FeatureFitStatus"/> order. Read from the method body at
    /// mzLib cc67b7a29 (FitOne): too few samples with a response of either kind, a rank-deficient
    /// design over them, fewer than two groups or none seen twice, separation (a coefficient beyond
    /// <see cref="LogitNormalMixedModel.SeparationLimit"/>, a Hessian that is not positive definite,
    /// or nothing measured and every group one-sided), a variance search that did not converge, or
    /// fitted. It never returns <see cref="FeatureFitStatus.VarianceRatioAtLimit"/>. mzLib does not
    /// yet list these per model (mzLib#1487); until it does, this list is the bridge's reading of the
    /// code, held by a test.
    /// </summary>
    internal static readonly FeatureFitStatus[] LatentMixedFitStatuses =
    {
        FeatureFitStatus.Fitted, FeatureFitStatus.TooFewObservations, FeatureFitStatus.RankDeficient,
        FeatureFitStatus.Separated, FeatureFitStatus.NotConverged, FeatureFitStatus.TooFewGroups,
    };

    /// <summary>The columns a <c>stats fit-link</c> cells file names; <c>depth</c> is optional.</summary>
    internal static readonly string[] LinkCellColumns = { "latent_mean", "latent_sd", "successes", "trials", "depth" };

    /// <summary>
    /// <c>stats fit-latent-mixed --design FILE --groups FILE --link identity|A,B[,C[,RHO]]
    /// [--responses FILE] [--counts FILE] [--depth FILE] [--residual-variance FILE]
    /// [--estimator reml|ml] [--threads N]</c>, with the coefficients to report on stdin.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>responses</b> is <c>stats fit-mixed</c>'s table, with values on the latent scale (e.g. the
    /// log-odds of a site's occupancy). <b>counts</b> has the same layout, each cell
    /// <c>successes/trials</c> or blank. At least one of the two is required; one that is absent is
    /// an all-missing matrix of the other's shape (nothing measured, or nothing counted). When both
    /// are given their features and samples must be the same, in any order.
    /// </para>
    /// <para>
    /// <b>link</b> is required because mzLib's <c>link</c> parameter has no default. <c>identity</c>
    /// is mzLib's <see cref="DetectionLink.Identity"/>; two to four numbers construct a
    /// <see cref="DetectionLink"/> with those positional arguments, so an omitted depth slope or
    /// overdispersion is the record's own default. <b>estimator</b>, <b>depth</b> and
    /// <b>residual-variance</b> are passed only when given.
    /// </para>
    /// </remarks>
    public static object FitLatentMixed(Program.Arguments arguments)
    {
        string? responsesPath = OptionalFile(arguments, "responses");
        string? countsPath = OptionalFile(arguments, "counts");
        if (responsesPath is null && countsPath is null)
            throw new Program.UsageException(
                "Give --responses (measured latent values), --counts (successes/trials), or both.");
        string designPath = RequireFile(arguments, "design");
        string groupsPath = RequireFile(arguments, "groups");
        string? depthPath = OptionalFile(arguments, "depth");
        string? sigmaPath = OptionalFile(arguments, "residual-variance");
        DetectionLink link = ParseLink(arguments);
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
        ResponseTable? measured = responsesPath is null ? null : ReadResponses(responsesPath);
        CountTable? counted = countsPath is null ? null : ReadCounts(countsPath);
        ResponseTable table = AlignLatentInputs(measured, counted);
        BinomialCount?[,] counts = counted is null
            ? new BinomialCount?[table.FeatureIds.Count, table.SampleNames.Count]
            : counted.Aligned(table.FeatureIds, table.SampleNames);
        DesignTable design = ReadDesign(designPath, table.SampleNames);
        GroupTable groups = ReadGroups(groupsPath, table.SampleNames);
        double[]? depth = depthPath is null ? null : ReadSampleValues(depthPath, table.SampleNames);
        double[]? sigma2 = sigmaPath is null ? null : ReadFeatureVariances(sigmaPath, table.FeatureIds);
        CheckTestedCoefficients(tested, design);

        LogitNormalMixedModelFit fit;
        try
        {
            // The estimator is passed only when the caller named one, so an absent one is mzLib's default.
            fit = estimator is { } chosen
                ? LogitNormalMixedModel.Fit(table.Values, counts, design.Values, groups.Labels, link, depth, sigma2,
                    design.CoefficientNames, chosen, threads)
                : LogitNormalMixedModel.Fit(table.Values, counts, design.Values, groups.Labels, link, depth, sigma2,
                    design.CoefficientNames, maxThreads: threads);
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            // A cell both measured and counted, a count that is not one, a link mzLib refuses, a depth
            // or residual variance it refuses, a design not of full rank, an empty group label: the
            // caller's input, refused by mzLib, so a usage error like every malformed argument.
            throw new Program.UsageException(LatentRefusal(invalid, table));
        }

        var rows = new Sdrf.ColumnTable(LatentMixedColumns);
        int untestable = 0;
        foreach (string coefficient in tested)
        {
            int j = fit.IndexOf(coefficient);
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                double df = fit.DegreesOfFreedom(f, j);
                if (fit.Status[f] == FeatureFitStatus.Fitted && !double.IsFinite(df))
                    untestable++;
                rows.Add(
                    table.FeatureIds[f],
                    coefficient,
                    StatusName(fit.Status[f]),
                    Finite(fit.Coefficient(f, j)),
                    Finite(fit.StandardError(f, j)),
                    Finite(df),
                    Finite(fit.TValue(f, j)),
                    Finite(fit.PValue(f, j)),
                    Finite(fit.ResidualVariance[f]),
                    fit.ResidualVarianceGiven[f],
                    Finite(fit.GroupVariance[f]),
                    Finite(fit.LogLikelihood[f]),
                    fit.Measured[f],
                    fit.Counted[f],
                    fit.Groups[f]);
            }
        }

        int cells = table.FeatureIds.Count * table.SampleNames.Count;
        int measuredCells = 0, countedCells = 0;
        for (int f = 0; f < table.FeatureIds.Count; f++)
            for (int s = 0; s < table.SampleNames.Count; s++)
            {
                if (double.IsFinite(table.Values[f, s])) measuredCells++;
                if (counts[f, s] is not null) countedCells++;
            }
        var statusCounts = LatentMixedFitStatuses.ToDictionary(StatusName, s => fit.Status.Count(x => x == s));
        List<string> groupNames = groups.Labels.Distinct(StringComparer.Ordinal).ToList();

        return new
        {
            responses_file = responsesPath is null ? null : Path.GetFullPath(responsesPath),
            counts_file = countsPath is null ? null : Path.GetFullPath(countsPath),
            design_file = Path.GetFullPath(designPath),
            groups_file = Path.GetFullPath(groupsPath),
            depth_file = depthPath is null ? null : Path.GetFullPath(depthPath),
            residual_variance_file = sigmaPath is null ? null : Path.GetFullPath(sigmaPath),
            grouping = groups.Name,
            group_names = groupNames,
            group_count = groupNames.Count,
            feature_count = fit.FeatureCount,
            sample_count = table.SampleNames.Count,
            sample_names = table.SampleNames,
            coefficient_names = design.CoefficientNames,
            tested_coefficients = tested,
            estimator = EstimatorName(fit.Estimator),
            link = LinkWire(fit.Link),
            threads,
            measured_cells = measuredCells,
            counted_cells = countedCells,
            missing_count = cells - measuredCells - countedCells,
            status_counts = statusCounts,
            row_count = rows.RowCount,
            column_names = rows.Names,
            columns = rows.Columns(),
            caveats = LatentMixedCaveats(fit, untestable, countedCells),
        };
    }

    /// <summary>
    /// <c>stats fit-link --cells FILE [--with-depth true|false] [--with-overdispersion true|false]
    /// [--detection-only true|false]</c>: <see cref="LogitNormalMixedModel.FitLink"/> over the cells
    /// file's rows. Each option is mzLib's parameter of the same name and is passed only when given;
    /// the envelope echoes the value used, mzLib's declared default when it was absent.
    /// </summary>
    public static object FitLink(Program.Arguments arguments)
    {
        string cellsPath = RequireFile(arguments, "cells");
        bool? withDepth = OptionalBool(arguments, "with-depth");
        bool? withOverdispersion = OptionalBool(arguments, "with-overdispersion");
        bool? detectionOnly = OptionalBool(arguments, "detection-only");
        var (cells, depthGiven) = ReadLinkCells(cellsPath);

        DetectionLink? link;
        try
        {
            link = LogitNormalMixedModel.FitLink(cells,
                withDepth ?? FitLinkDefault<bool>("withDepth"),
                withOverdispersion ?? FitLinkDefault<bool>("withOverdispersion"),
                detectionOnly ?? FitLinkDefault<bool>("detectionOnly"));
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            throw new Program.UsageException(invalid.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0]);
        }

        bool depthUsed = withDepth ?? FitLinkDefault<bool>("withDepth");
        bool overdispersionUsed = withOverdispersion ?? FitLinkDefault<bool>("withOverdispersion");
        bool detectionUsed = detectionOnly ?? FitLinkDefault<bool>("detectionOnly");
        return new
        {
            cells_file = Path.GetFullPath(cellsPath),
            cell_count = cells.Count,
            depth_given = depthGiven,
            with_depth = depthUsed,
            with_overdispersion = overdispersionUsed,
            detection_only = detectionUsed,
            link = link is null ? null : LinkWire(link),
            caveats = LinkCaveats(link, depthUsed, overdispersionUsed, detectionUsed, depthGiven),
        };
    }

    /// <summary>
    /// The default mzLib declares for a parameter of <see cref="LogitNormalMixedModel.FitLink"/>, so an
    /// option the caller leaves out is mzLib's value at the pin, read rather than retyped.
    /// </summary>
    internal static T FitLinkDefault<T>(string parameter) =>
        (T)typeof(LogitNormalMixedModel).GetMethod(nameof(LogitNormalMixedModel.FitLink))!.GetParameters()
            .Single(p => p.Name == parameter).DefaultValue!;

    /// <summary>A <see cref="DetectionLink"/> as the wire writes it.</summary>
    private static object LinkWire(DetectionLink link) => new
    {
        intercept = link.Intercept,
        slope = link.Slope,
        depth_slope = link.DepthSlope,
        overdispersion = link.Overdispersion,
    };

    // ---- inputs ----------------------------------------------------------------------------------

    /// <summary>A <c>--name FILE</c> that may be absent; given, it must exist and have a value.</summary>
    private static string? OptionalFile(Program.Arguments arguments, string name)
    {
        RequireValueIfProvided(arguments, name);
        string? path = arguments.Optional(name);
        if (path is null)
            return null;
        if (path.Trim().Length == 0)
            throw new Program.UsageException($"Option --{name} was given without a value.");
        if (!File.Exists(path))
            throw new Program.UsageException($"The --{name} file was not found: '{path}'.");
        return path;
    }

    /// <summary>A <c>--name true|false</c> that may be absent (null: not passed to mzLib).</summary>
    private static bool? OptionalBool(Program.Arguments arguments, string name)
    {
        RequireValueIfProvided(arguments, name);
        return arguments.Optional(name)?.Trim().ToLowerInvariant() switch
        {
            null => null,
            "true" => true,
            "false" => false,
            var other => throw new Program.UsageException($"Option --{name} must be 'true' or 'false'; got '{other}'."),
        };
    }

    /// <summary>
    /// <c>--link</c>: <c>identity</c> (mzLib's <see cref="DetectionLink.Identity"/>), or two to four
    /// comma-separated numbers, the <see cref="DetectionLink"/>'s positional arguments (intercept,
    /// slope, depth slope, overdispersion) in natural log-odds. mzLib validates the values.
    /// </summary>
    private static DetectionLink ParseLink(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "link");
        string raw = arguments.Required("link").Trim();
        if (raw.Equals("identity", StringComparison.OrdinalIgnoreCase))
            return DetectionLink.Identity;
        string[] parts = raw.Split(',');
        if (parts.Length is < 2 or > 4)
            throw new Program.UsageException(
                $"Option --link must be 'identity' or 2 to 4 comma-separated numbers: intercept,slope[,depth_slope[,overdispersion]]; got '{raw}'.");
        string[] names = { "intercept", "slope", "depth_slope", "overdispersion" };
        var v = parts.Select((p, i) =>
            double.TryParse(p.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double x)
                ? x
                : throw new Program.UsageException($"Option --link: {names[i]} must be a number; got '{p.Trim()}'.")).ToArray();
        // Only the values given are passed, so an omitted depth slope or overdispersion is the record's default.
        return v.Length switch
        {
            2 => new DetectionLink(v[0], v[1]),
            3 => new DetectionLink(v[0], v[1], v[2]),
            _ => new DetectionLink(v[0], v[1], v[2], v[3]),
        };
    }

    /// <summary>The counts table: feature ids, sample names, and each cell's count (null = not counted).</summary>
    private sealed record CountTable(List<string> FeatureIds, List<string> SampleNames, BinomialCount?[,] Values)
    {
        /// <summary>The counts in the given feature and sample order (both already checked to be this table's).</summary>
        public BinomialCount?[,] Aligned(List<string> features, List<string> samples)
        {
            if (features.SequenceEqual(FeatureIds) && samples.SequenceEqual(SampleNames))
                return Values;
            var row = FeatureIds.Select((id, i) => (id, i)).ToDictionary(x => x.id, x => x.i, StringComparer.Ordinal);
            var column = SampleNames.Select((s, i) => (s, i)).ToDictionary(x => x.s, x => x.i, StringComparer.Ordinal);
            var aligned = new BinomialCount?[features.Count, samples.Count];
            for (int f = 0; f < features.Count; f++)
                for (int s = 0; s < samples.Count; s++)
                    aligned[f, s] = Values[row[features[f]], column[samples[s]]];
            return aligned;
        }
    }

    /// <summary>
    /// A feature-by-sample table laid out as the responses table, each cell <c>successes/trials</c>
    /// (two whole numbers) or blank, NA or NaN for not counted. mzLib checks that a count is one.
    /// </summary>
    private static CountTable ReadCounts(string path)
    {
        List<string[]> rows = ReadTsv(path, "counts");
        string[] header = rows[0];
        if (header.Length < 2)
            throw new Program.UsageException(
                $"The counts file '{path}' has {header.Length} column; it needs a feature id column and at least one sample column.");
        List<string> samples = header.Skip(1).Select(h => h.Trim()).ToList();
        RequireUnique(samples, "sample column", "counts header");
        if (rows.Count == 1)
            throw new Program.UsageException($"The counts file '{path}' has a header but no feature rows.");
        var ids = new List<string>(rows.Count - 1);
        var values = new BinomialCount?[rows.Count - 1, samples.Count];
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the counts file has {row.Length} fields; the header has {header.Length}.");
            ids.Add(row[0].Trim());
            for (int s = 0; s < samples.Count; s++)
            {
                string cell = row[s + 1].Trim();
                if (cell.Length == 0 || cell.Equals("NA", StringComparison.OrdinalIgnoreCase)
                    || cell.Equals("NaN", StringComparison.OrdinalIgnoreCase))
                    continue;
                string[] kn = cell.Split('/');
                if (kn.Length != 2
                    || !int.TryParse(kn[0].Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int k)
                    || !int.TryParse(kn[1].Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int n))
                    throw new Program.UsageException(
                        $"Line {r + 1} of the counts file, sample '{samples[s]}': '{cell}' is not a count. Write successes/trials (3/8), or leave it blank for a sample not counted.");
                values[r - 1, s] = new BinomialCount(k, n);
            }
        }
        RequireUnique(ids, "feature id", "counts file's first column");
        return new CountTable(ids, samples, values);
    }

    /// <summary>
    /// The feature ids, sample names and measured matrix the fit runs on: the responses table's, or,
    /// with counts alone, the counts table's with nothing measured. When both are given they must
    /// list the same features and samples, in any order; the responses table's order is kept.
    /// </summary>
    private static ResponseTable AlignLatentInputs(ResponseTable? measured, CountTable? counted)
    {
        if (measured is null)
        {
            var empty = new double[counted!.FeatureIds.Count, counted.SampleNames.Count];
            for (int f = 0; f < empty.GetLength(0); f++)
                for (int s = 0; s < empty.GetLength(1); s++)
                    empty[f, s] = double.NaN;
            return new ResponseTable(counted.FeatureIds, counted.SampleNames, empty, empty.Length, 0);
        }
        if (counted is null)
            return measured;
        static string Difference(string what, List<string> expected, List<string> actual)
        {
            List<string> absent = expected.Where(x => !actual.Contains(x, StringComparer.Ordinal)).ToList();
            List<string> extra = actual.Where(x => !expected.Contains(x, StringComparer.Ordinal)).ToList();
            return (absent.Count > 0 ? $" {what} not in the counts: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? $" {what} in the counts but not in the responses: " + string.Join(", ", extra) + "." : "");
        }
        string features = Difference("Features", measured.FeatureIds, counted.FeatureIds);
        string samples = Difference("Samples", measured.SampleNames, counted.SampleNames);
        if (features.Length > 0 || samples.Length > 0)
            throw new Program.UsageException(
                "The counts table must list exactly the responses' features and sample columns." + features + samples);
        return measured;
    }

    /// <summary>
    /// A TSV with a header row and one row per sample: the sample name, then exactly one numeric
    /// column. The samples must be exactly the fit's, in any order.
    /// </summary>
    private static double[] ReadSampleValues(string path, List<string> sampleNames)
    {
        List<string[]> rows = ReadTsv(path, "depth");
        if (rows[0].Length != 2)
            throw new Program.UsageException(
                $"The depth file '{path}' has {rows[0].Length} column{(rows[0].Length == 1 ? "" : "s")}; it needs exactly two: the sample name, then its depth.");
        var bySample = new Dictionary<string, double>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != 2)
                throw new Program.UsageException($"Line {r + 1} of the depth file has {row.Length} fields; the header has 2.");
            string sample = row[0].Trim();
            if (bySample.ContainsKey(sample))
                throw new Program.UsageException($"Sample '{sample}' appears more than once in the depth file.");
            double value = ParseCell(row[1], out bool missing);
            if (!missing && double.IsNaN(value))
                throw new Program.UsageException($"Line {r + 1} of the depth file, sample '{sample}': '{row[1].Trim()}' is not a number.");
            bySample[sample] = value;
        }
        List<string> absent = sampleNames.Where(s => !bySample.ContainsKey(s)).ToList();
        List<string> extra = bySample.Keys.Where(s => !sampleNames.Contains(s, StringComparer.Ordinal)).ToList();
        if (absent.Count > 0 || extra.Count > 0)
            throw new Program.UsageException(
                "The depth file's samples must be exactly the fit's sample columns." +
                (absent.Count > 0 ? " Not in the depth file: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? " In the depth file but not in the samples: " + string.Join(", ", extra) + "." : ""));
        return sampleNames.Select(s => bySample[s]).ToArray();
    }

    /// <summary>
    /// A TSV with a header row and one row per feature: the feature id, then sigma² to use for it.
    /// Blank, NA or NaN, or a feature not listed, is NaN, which mzLib reads as "estimate it". An id
    /// that is not one of the fit's features is refused.
    /// </summary>
    private static double[] ReadFeatureVariances(string path, List<string> featureIds)
    {
        List<string[]> rows = ReadTsv(path, "residual-variance");
        if (rows[0].Length != 2)
            throw new Program.UsageException(
                $"The residual-variance file '{path}' has {rows[0].Length} column{(rows[0].Length == 1 ? "" : "s")}; it needs exactly two: the feature id, then its residual variance.");
        var byFeature = new Dictionary<string, double>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != 2)
                throw new Program.UsageException($"Line {r + 1} of the residual-variance file has {row.Length} fields; the header has 2.");
            string feature = row[0].Trim();
            if (byFeature.ContainsKey(feature))
                throw new Program.UsageException($"Feature '{feature}' appears more than once in the residual-variance file.");
            double value = ParseCell(row[1], out bool missing);
            if (!missing && double.IsNaN(value))
                throw new Program.UsageException(
                    $"Line {r + 1} of the residual-variance file, feature '{feature}': '{row[1].Trim()}' is not a number. Leave it blank to have it estimated.");
            byFeature[feature] = value;
        }
        List<string> unknown = byFeature.Keys.Where(f => !featureIds.Contains(f, StringComparer.Ordinal)).ToList();
        if (unknown.Count > 0)
            throw new Program.UsageException(
                "The residual-variance file names features the fit does not have: " + string.Join(", ", unknown) + ".");
        return featureIds.Select(f => byFeature.TryGetValue(f, out double v) ? v : double.NaN).ToArray();
    }

    /// <summary>
    /// The cells file: a TSV naming latent_mean, latent_sd, successes and trials, and optionally depth,
    /// in any order and no other column. A cell without a depth column is constructed without one, so
    /// its depth is <see cref="LinkCell"/>'s own default.
    /// </summary>
    private static (List<LinkCell> Cells, bool DepthGiven) ReadLinkCells(string path)
    {
        var (rows, index) = ReadCandidateTable(path, "cells", LinkCellColumns, optional: new[] { "depth" });
        if (rows.Count == 0)
            throw new Program.UsageException($"The cells file '{path}' has a header but no cells.");
        bool depthGiven = index.ContainsKey("depth");
        var cells = new List<LinkCell>(rows.Count);
        foreach (var (line, row) in rows)
        {
            string where = $"Line {line} of the cells file";
            double Number(string column)
            {
                double value = ParseCell(row[index[column]], out bool missing);
                if (!missing && double.IsNaN(value))
                    throw new Program.UsageException($"{where}: {column} '{row[index[column]].Trim()}' is not a number.");
                return value;
            }
            int Whole(string column) =>
                int.TryParse(row[index[column]].Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int v)
                    ? v
                    : throw new Program.UsageException($"{where}: {column} '{row[index[column]].Trim()}' is not a whole number.");
            double mean = Number("latent_mean"), sd = Number("latent_sd");
            int successes = Whole("successes"), trials = Whole("trials");
            cells.Add(depthGiven
                ? new LinkCell(mean, sd, successes, trials, Number("depth"))
                : new LinkCell(mean, sd, successes, trials));
        }
        return (cells, depthGiven);
    }

    /// <summary>
    /// mzLib's refusal, as <see cref="MzLibRefusal"/> words it, but matching mzLib's capitalised
    /// "Feature 3, sample 5" too (LogitNormalMixedModel's messages), so the cell is also told by the
    /// names the caller wrote.
    /// </summary>
    private static string LatentRefusal(ArgumentException invalid, ResponseTable table)
    {
        string message = invalid.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0];
        var cell = Regex.Match(message, @"[Ff]eature (\d+), sample (\d+)");
        if (cell.Success
            && int.TryParse(cell.Groups[1].Value, out int f) && f < table.FeatureIds.Count
            && int.TryParse(cell.Groups[2].Value, out int s) && s < table.SampleNames.Count)
            return $"{message} That is feature '{table.FeatureIds[f]}', sample '{table.SampleNames[s]}'.";
        var sample = Regex.Match(message, @"^Sample (\d+) ");
        if (sample.Success && int.TryParse(sample.Groups[1].Value, out int t) && t < table.SampleNames.Count)
            return $"{message} That is sample '{table.SampleNames[t]}'.";
        return message;
    }

    // ---- caveats ---------------------------------------------------------------------------------

    private static string RoundTrip(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    private static List<string> LatentMixedCaveats(LogitNormalMixedModelFit fit, int untestable, int countedCells)
    {
        bool reml = fit.Estimator == VarianceEstimator.Reml;
        DetectionLink l = fit.Link;
        var caveats = new List<string>
        {
            "Per feature, the latent is X beta + u_group + e, with u ~ N(0, group_variance) and e ~ N(0, residual_variance). A measured sample reads the latent. A counted sample has successes of trials succeed, each with chance g, where logit g = intercept + slope x latent + depth_slope x depth (the link), binomially, or beta-binomially with intra-cell correlation overdispersion. One random intercept per group: no random slopes, nested or crossed groupings.",
            $"Responses, estimates and the link are on the latent scale, in natural log-odds when the latent is a logit (e.g. a site's occupancy): send logit values, not fractions or intensities. This fit read its counts through intercept {RoundTrip(l.Intercept)}, slope {RoundTrip(l.Slope)}, depth_slope {RoundTrip(l.DepthSlope)}, overdispersion {RoundTrip(l.Overdispersion)}; the link and which cells are counted are the caller's choices, not the bridge's (stats fit-link fits a link).",
            "Counted samples are integrated by adaptive Gauss-Hermite quadrature, 40 nodes over the group intercept and 12 over each count's residual, as glmer's nAGQ; measured samples exactly. mzLib agrees with lme4's glmer (25-point quadrature, ML) to about 1e-6, with lmer when nothing is counted to about 1e-5, and with glmmTMB to its Laplace error.",
            "Fixed effects are tested by Wald t on nlme's containment degrees of freedom, over the samples with a response of either kind. standard_error is sqrt of the diagonal of the inverse Hessian of the log-likelihood in beta at the estimated variances; it does not include the variances' own uncertainty.",
            reml
                ? "reml: the fixed effects are integrated out by a Laplace approximation, as glmmTMB's REML = TRUE does; it is exactly the restricted likelihood when nothing is counted. log_likelihood is that criterion: compare it only between models with the same fixed effects."
                : "ml: log_likelihood is the full log-likelihood, comparable between designs fitted to the same samples. With nothing counted the variances are stats fit-mixed --estimator ml's, but standard errors differ: these are at the ML variances, fit-mixed's scale by RSS / (observed - coefficients) as nlme's summary does.",
            $"separated means no estimate exists: a coefficient ran beyond +-{RoundTrip(LogitNormalMixedModel.SeparationLimit)} on the latent scale (mzLib's SeparationLimit), the Hessian is not positive definite, or nothing was measured and every group's counts are one-sided (all 0, or all full). It is a status, never a p-value.",
            "p_value is not adjusted across features. Send one coefficient's p_value column to stats adjust for Benjamini-Hochberg.",
        };
        int notFitted = fit.Status.Count(s => s != FeatureFitStatus.Fitted);
        if (notFitted > 0)
        {
            string counts = string.Join(", ", LatentMixedFitStatuses.Where(s => s != FeatureFitStatus.Fitted)
                .Select(s => (Status: s, Count: fit.Status.Count(x => x == s)))
                .Where(c => c.Count > 0)
                .Select(c => $"{c.Count} {StatusName(c.Status)}"));
            caveats.Add($"{notFitted} of {fit.FeatureCount} features were not fitted ({counts}); their statistics are null.");
        }
        if (untestable > 0)
            caveats.Add(
                $"{untestable} fitted rows have no positive containment df: a group-level coefficient with too few groups to test it. Their df and p_value are null; estimate, standard_error and t are reported.");
        int boundary = Enumerable.Range(0, fit.FeatureCount)
            .Count(f => fit.Status[f] == FeatureFitStatus.Fitted && fit.GroupVariance[f] == 0);
        if (boundary > 0)
            caveats.Add($"{boundary} of the fitted features have group_variance 0 (the boundary): no between-group variance.");
        int given = fit.ResidualVarianceGiven.Count(g => g);
        if (given > 0)
            caveats.Add(
                $"{given} of {fit.FeatureCount} features were fitted at the residual_variance you gave (residual_variance_given), not an estimate; their standard errors do not include its uncertainty.");
        if (countedCells == 0)
            caveats.Add(reml
                ? "No cell was counted, so this is stats fit-mixed's model, fitted by LogitNormalMixedModel's own optimiser: its numbers agree with fit-mixed's to about 1e-5, not to the last digit."
                : "No cell was counted, so this is stats fit-mixed's model, fitted by LogitNormalMixedModel's own optimiser: estimates and variances agree with fit-mixed --estimator ml's to about 1e-5, standard errors do not (see the ml caveat).");
        return caveats;
    }

    private static List<string> LinkCaveats(DetectionLink? link, bool withDepth, bool withOverdispersion, bool detectionOnly,
        bool depthGiven)
    {
        var caveats = new List<string>
        {
            "Maximum likelihood over the cells: each cell's probability integrated over its latent normal (latent_mean, latent_sd) by 12-node adaptive Gauss-Hermite quadrature, summed over cells, maximized by Nelder-Mead. mzLib returns the point estimate only: no standard errors.",
            detectionOnly
                ? "detection_only: fitted on whether any trial succeeded (successes 0 against 1 or more), not on how many. mzLib's documentation: on PTM data, a link fitted on the full counts predicts covered-zeros worse."
                : "Fitted on the full counts, P(successes of trials), not only on whether any trial succeeded.",
            "The link is in natural log-odds: logit g = intercept + slope x latent + depth_slope x depth. Pass it to stats fit-latent-mixed --link. Which cells, and each cell's latent mean and SD (e.g. the mean and SD of its feature's measured cells), are the caller's modelling choices, not the bridge's.",
        };
        if (!withDepth)
            caveats.Add("depth_slope was not fitted, so it is 0.");
        if (!withOverdispersion)
            caveats.Add("overdispersion was not fitted, so it is 0: the trials of a cell are binomial.");
        if (withDepth && !depthGiven)
            caveats.Add("The depth slope was fitted, but the cells file has no depth column, so every cell's depth is mzLib's default 0 and the data say nothing about depth_slope.");
        if (link is null)
            caveats.Add("mzLib's search did not converge within its iteration limit (FitLink returned null), so link is null.");
        return caveats;
    }
}
