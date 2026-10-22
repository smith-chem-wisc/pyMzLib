using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats correlate</c>: Spearman's rank correlation of every feature of a feature-by-sample table
/// with one per-sample trait, optionally partial given per-sample covariates, by mzLib's
/// <see cref="SpearmanCorrelation"/> (mzLib #1367, 1.0.594; the partial form #1446, 1.0.595). The
/// bridge reads the tables and makes one mzLib call per feature; every number is mzLib's.
/// </summary>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats correlate</c> table, in order.</summary>
    internal static readonly string[] CorrelateColumns = { "feature", "rho", "p_value", "n", "has_ties", "method" };

    /// <summary>
    /// <c>stats correlate --responses FILE --trait FILE [--given FILE]</c>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>responses</b> is <c>stats fit</c>'s table: a header row, then one row per feature, the first
    /// column the feature id and every other column one sample. Blank, NA or NaN is missing.
    /// </para>
    /// <para>
    /// <b>trait</b> is a TSV with a header row and one row per sample: the sample name, then exactly
    /// one value column, whose header names the trait. <b>given</b> has the same shape with one or
    /// more covariate columns. Both must list exactly the responses' samples, in any order; a blank,
    /// NA or NaN cell is missing. mzLib drops a sample from a feature's correlation when any of its
    /// values is missing (pairwise-complete), and reports how many samples were used as <c>n</c>.
    /// </para>
    /// </remarks>
    public static object Correlate(Program.Arguments arguments)
    {
        string responsesPath = RequireFile(arguments, "responses");
        string traitPath = RequireFile(arguments, "trait");
        RequireValueIfProvided(arguments, "given");
        string? givenPath = arguments.WasProvided("given") ? RequireFile(arguments, "given") : null;

        ResponseTable responses = ReadResponses(responsesPath);
        SampleTable trait = ReadSampleTable(traitPath, "trait", responses.SampleNames, oneValueColumn: true);
        SampleTable? given = givenPath is null ? null : ReadSampleTable(givenPath, "given", responses.SampleNames, oneValueColumn: false);

        double[] traitValues = trait.Columns[0];
        var covariates = given?.Columns.Select(c => (IReadOnlyList<double>)c).ToList();

        var table = new Sdrf.ColumnTable(CorrelateColumns);
        var methodCounts = Enum.GetValues<SpearmanPValueMethod>().ToDictionary(SpearmanMethodName, _ => 0);
        int sampleCount = responses.SampleNames.Count;
        for (int f = 0; f < responses.FeatureIds.Count; f++)
        {
            var x = new double[sampleCount];
            for (int s = 0; s < sampleCount; s++)
                x[s] = responses.Values[f, s];
            SpearmanResult result = covariates is null
                ? SpearmanCorrelation.Correlate(x, traitValues)
                : SpearmanCorrelation.PartialCorrelate(x, traitValues, covariates);
            string method = SpearmanMethodName(result.Method);
            methodCounts[method]++;
            table.Add(responses.FeatureIds[f], Finite(result.Rho), Finite(result.PValue), result.N, result.HasTies, method);
        }

        return new
        {
            responses_file = Path.GetFullPath(responsesPath),
            trait_file = Path.GetFullPath(traitPath),
            given_file = givenPath is null ? null : Path.GetFullPath(givenPath),
            trait_name = trait.Names[0],
            covariate_names = given?.Names ?? new List<string>(),
            partial = given != null,
            feature_count = responses.FeatureIds.Count,
            sample_count = sampleCount,
            sample_names = responses.SampleNames,
            missing_count = responses.MissingCount,
            trait_missing_count = traitValues.Count(v => !double.IsFinite(v)),
            given_incomplete_count = given is null
                ? (int?)null
                : Enumerable.Range(0, sampleCount).Count(s => given.Columns.Any(c => !double.IsFinite(c[s]))),
            method_counts = methodCounts,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = CorrelateCaveats(given != null, methodCounts, responses.FeatureIds.Count),
        };
    }

    /// <summary>
    /// mzLib's <see cref="SpearmanPValueMethod"/> as the wire writes it. Every value is named here, so
    /// a value mzLib adds fails a test on the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string SpearmanMethodName(SpearmanPValueMethod method) => method switch
    {
        SpearmanPValueMethod.Exact => "exact",
        SpearmanPValueMethod.Asymptotic => "asymptotic",
        SpearmanPValueMethod.NotEstimable => "not_estimable",
        _ => throw new InvalidOperationException(
            $"SpearmanPValueMethod.{method} has no wire name. Name it in Statistics.SpearmanMethodName."),
    };

    /// <summary>A per-sample table: value columns by header, each in the responses' sample order.</summary>
    private sealed record SampleTable(List<string> Names, List<double[]> Columns);

    /// <summary>
    /// A TSV of one row per sample, the sample name first and then exactly one numeric column (with
    /// <paramref name="oneValueColumn"/>, the trait) or one or more (the covariates). Rows are matched
    /// to the responses' samples by name, which they must be exactly. A blank, NA or NaN cell is
    /// missing (NaN); any other cell that is not a number is refused.
    /// </summary>
    private static SampleTable ReadSampleTable(string path, string what, List<string> sampleNames, bool oneValueColumn)
    {
        List<string[]> rows = ReadTsv(path, what);
        string[] header = rows[0];
        if (header.Length < 2)
            throw new Program.UsageException(
                $"The {what} file '{path}' has {header.Length} column; it needs a sample name column and at least one value column.");
        List<string> names = header.Skip(1).Select(h => h.Trim()).ToList();
        if (oneValueColumn && names.Count != 1)
            throw new Program.UsageException(
                $"The {what} file '{path}' has {names.Count} value columns ({string.Join(", ", names)}); " +
                "it needs exactly one: the sample name, then the trait. Put covariates in --given.");
        RequireUnique(names, "column", $"{what} header");

        var bySample = new Dictionary<string, double[]>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the {what} file has {row.Length} fields; the header has {header.Length}.");
            string sample = row[0].Trim();
            if (bySample.ContainsKey(sample))
                throw new Program.UsageException($"Sample '{sample}' appears more than once in the {what} file.");
            var cells = new double[names.Count];
            for (int c = 0; c < names.Count; c++)
                cells[c] = ParseOptional(row[c + 1], $"The {what} file, sample '{sample}', column '{names[c]}'");
            bySample[sample] = cells;
        }

        List<string> absent = sampleNames.Where(s => !bySample.ContainsKey(s)).ToList();
        List<string> extra = bySample.Keys.Where(s => !sampleNames.Contains(s, StringComparer.Ordinal)).ToList();
        if (absent.Count > 0 || extra.Count > 0)
            throw new Program.UsageException(
                $"The {what} file's samples must be exactly the responses' sample columns." +
                (absent.Count > 0 ? $" Not in the {what} file: " + string.Join(", ", absent) + "." : "") +
                (extra.Count > 0 ? $" In the {what} file but not in the responses: " + string.Join(", ", extra) + "." : ""));

        var columns = names.Select((_, c) => sampleNames.Select(s => bySample[s][c]).ToArray()).ToList();
        return new SampleTable(names, columns);
    }

    private static List<string> CorrelateCaveats(bool partial, Dictionary<string, int> methodCounts, int featureCount)
    {
        var caveats = new List<string>();
        if (partial)
        {
            caveats.Add(
                "Partial Spearman (mzLib's SpearmanCorrelation.PartialCorrelate, as ppcor::pcor.test(method = 'spearman')): every variable is ranked over the samples complete in all of them, and rho is the correlation of the feature's and the trait's ranks after each is regressed on the covariates' ranks. p_value is two-sided, from t with n - 2 - k degrees of freedom for k covariates, so method is always asymptotic.");
            caveats.Add(
                "mzLib drops a covariate that the others carry entirely (one that ranks the samples as they already do), and k counts the covariates kept. mzLib does not report k, so neither does this verb.");
            caveats.Add(
                "Give only nuisance covariates (batch, run order, loading). A covariate that follows the trait takes the trait's own signal with it: rho then measures what is left within its levels.");
        }
        else
        {
            caveats.Add(
                "Spearman's rho with a two-sided test of rho = 0 (mzLib's SpearmanCorrelation.Correlate). p_value is exact, from the full permutation distribution, for n <= 9 without ties, as R's cor.test(method = 'spearman'); otherwise it is asymptotic, t with n - 2 degrees of freedom, as cor.test(exact = FALSE).");
            caveats.Add(
                "For 10 <= n < 1290 without ties, R's cor.test default is an Edgeworth series (AS 89), which mzLib does not use, so p-values there differ slightly from R's default output.");
        }
        caveats.Add(
            "A sample enters a feature's correlation only when its response" + (partial ? ", its trait and every covariate are" : " and its trait are both") +
            " finite: missing values are omitted, never imputed. n is the number of samples used. Ties take their average rank.");
        caveats.Add(
            "p_value is not adjusted across features. Send the p_value column to stats adjust for Benjamini-Hochberg.");
        int notEstimable = methodCounts[SpearmanMethodName(SpearmanPValueMethod.NotEstimable)];
        if (notEstimable > 0)
            caveats.Add(
                $"{notEstimable} of {featureCount} features are not_estimable: fewer than {(partial ? "k + 3 complete samples (k the covariates kept, so at least 4)" : "3 complete samples")}, a response{(partial ? ", trait or covariate" : " or trait")} that is constant over them" +
                (partial ? ", or covariates that carry the response or the trait entirely" : "") + ". Their rho and p_value are null.");
        return caveats;
    }
}
