using System.Globalization;
using Quantification.Differential;

namespace MzLibBridge;

/// <summary>
/// <c>diff design</c>: a sample sheet and a model, turned by mzLib's <see cref="AnalysisDesign"/>
/// (mzLib #1450) into one design per stratum: its matrix, rank, aliased terms, warnings, the
/// contrasts it cannot run, and every default contrast with its estimability.
/// </summary>
/// <remarks>
/// <para>
/// <b>The sample sheet is any table of samples.</b> Its first column is the sample id; the model on
/// stdin says which other columns are factors, covariates, strata or the batch. Columns it does not
/// name are ignored and listed in <c>unused_columns</c>, so a sheet written by <c>sdrf samples</c>
/// or a spreadsheet works unedited.
/// </para>
/// <para>
/// <b>The model is one term per stdin line,</b> tab-separated, because every name in it came from a
/// TSV cell and so cannot hold a tab. An empty optional field is "not given", which leaves mzLib's
/// default in place: the bridge constructs mzLib's <see cref="FactorSpec"/> and
/// <see cref="CovariateSpec"/> without that argument. <see cref="DesignSpec.FitBatch"/> is never
/// set, so batch enters exactly when a batch column is named (mzLib's default true).
/// </para>
/// <para>
/// <b>The pivot files.</b> With <c>--out DIR</c>, stratum i (in mzLib's order, ordinal by key) gets
/// <c>design_i.tsv</c>, <see cref="AnalysisDesign.FittedMatrix"/> under a <c>sample</c> column, which
/// <c>stats fit --design</c> reads as it is; and, when any contrast is estimable,
/// <c>contrasts_i.tsv</c>, each estimable contrast's id and its
/// <see cref="AnalysisDesign.FittedWeights"/> on the same columns: the <c>ContrastWeights</c> mzLib's
/// <c>EmpiricalBayes.ModerateContrasts</c> takes. Files are numbered rather than named by stratum
/// key because a key can hold any character a sample sheet cell can.
/// </para>
/// </remarks>
internal static partial class Differential
{
    /// <summary>Row columns of the <c>diff design</c> table: one row per (stratum, default contrast).</summary>
    internal static readonly string[] DesignColumns =
    {
        "stratum", "contrast_id", "label", "numerator", "denominator", "covariate", "covariate_unit",
        "covariate_scale", "estimability", "weights", "fitted_weights",
    };

    /// <summary>mzLib's estimability verdicts, the keys of <c>estimability_counts</c>.</summary>
    internal static readonly string[] Estimabilities =
    {
        AnalysisDesign.Estimable, AnalysisDesign.NotEstimableConfounded, AnalysisDesign.NotEstimableRankDeficient,
    };

    private static readonly string[] TermKinds = { "factor", "covariate", "interaction", "stratum", "batch" };

    /// <summary>The model read from stdin, in the order given.</summary>
    private sealed class Model
    {
        public List<FactorSpec> Factors { get; } = new();
        public List<CovariateSpec> Covariates { get; } = new();
        public List<(string First, string Second)> Interactions { get; } = new();
        public List<string> Strata { get; } = new();
        public string? Batch { get; set; }

        /// <summary>The covariates' optional fields as given (null = not given, so mzLib's default).</summary>
        public List<object> CovariateEcho { get; } = new();

        /// <summary>Every sheet column a term names, with the role it was named in.</summary>
        public Dictionary<string, string> Roles { get; } = new(StringComparer.Ordinal);
    }

    /// <summary>
    /// <c>diff design --samples FILE [--out DIR]</c>, with the model on stdin, one term per line:
    /// <c>factor NAME [REFERENCE]</c>, <c>covariate NAME UNIT [SCALE [SCALED_UNIT [CENTRE]]]</c>,
    /// <c>interaction FIRST SECOND</c>, <c>stratum NAME</c>, <c>batch NAME</c> (fields tab-separated).
    /// </summary>
    public static object Design(Program.Arguments arguments)
    {
        string samplesPath = RequireFile(arguments, "samples");
        string? outDirectory = arguments.Optional("out");
        if (arguments.WasProvided("out") && string.IsNullOrWhiteSpace(outDirectory))
            throw new Program.UsageException("Option --out needs a directory to write the design and contrast tables into.");

        Model model = ReadModel(Program.ReadStdinLines());
        (List<AnalysisSample> samples, List<string> unused) = ReadSampleSheet(samplesPath, model);

        var spec = new DesignSpec
        {
            Factors = model.Factors,
            Covariates = model.Covariates,
            Interactions = model.Interactions,
        };

        IReadOnlyList<AnalysisDesign> designs;
        var contrastsByStratum = new List<IReadOnlyList<Contrast>>();
        try
        {
            designs = AnalysisDesign.ByStratum(samples, spec);
            foreach (AnalysisDesign design in designs)
                contrastsByStratum.Add(design.DefaultContrasts());
        }
        catch (ArgumentException refusal) when (refusal is not ArgumentNullException)
        {
            throw Refused(refusal);
        }

        if (outDirectory is not null)
            Directory.CreateDirectory(outDirectory);

        var table = new Sdrf.ColumnTable(DesignColumns);
        var estimability = Estimabilities.ToDictionary(e => e, _ => 0, StringComparer.Ordinal);
        var strata = new List<object>();
        int filesWritten = 0;
        for (int i = 0; i < designs.Count; i++)
        {
            AnalysisDesign design = designs[i];
            List<string> fittedNames = design.FittedColumns.Select(j => design.ColumnNames[j]).ToList();
            var estimable = new List<(Contrast Contrast, double[] Weights)>();
            foreach (Contrast contrast in contrastsByStratum[i])
            {
                string verdict = design.Estimability(contrast);
                estimability[verdict]++;
                double[]? fitted = verdict == AnalysisDesign.Estimable ? design.FittedWeights(contrast) : null;
                if (fitted is not null)
                    estimable.Add((contrast, fitted));
                table.Add(design.Stratum, contrast.Id, contrast.Label, contrast.Numerator, contrast.Denominator,
                    contrast.Covariate, contrast.CovariateUnit, contrast.CovariateScale, verdict,
                    contrast.Weights.ToList(), fitted?.ToList());
            }

            string? designFile = null, contrastsFile = null;
            if (outDirectory is not null)
            {
                double[,] fittedMatrix = design.FittedMatrix;
                designFile = Path.GetFullPath(Path.Combine(outDirectory, $"design_{i + 1}.tsv"));
                WritePlain(designFile, new[] { "sample" }.Concat(fittedNames),
                    design.Samples.Select((s, k) =>
                        new[] { s.SampleId }.Concat(Enumerable.Range(0, fittedNames.Count).Select(j => Number(fittedMatrix[k, j])))));
                filesWritten++;
                if (estimable.Count > 0)
                {
                    contrastsFile = Path.GetFullPath(Path.Combine(outDirectory, $"contrasts_{i + 1}.tsv"));
                    WritePlain(contrastsFile, new[] { "contrast" }.Concat(fittedNames),
                        estimable.Select(e => new[] { e.Contrast.Id }.Concat(e.Weights.Select(Number))));
                    filesWritten++;
                }
            }

            strata.Add(new
            {
                stratum = design.Stratum,
                sample_count = design.Samples.Count,
                sample_ids = design.Samples.Select(s => s.SampleId).ToList(),
                column_names = design.ColumnNames.ToList(),
                fitted_column_names = fittedNames,
                rank = design.Rank,
                aliased_terms = design.AliasedTerms.ToList(),
                warnings = design.Warnings.ToList(),
                not_run = design.NotRun.Select(n => new { id = n.Id, label = n.Label, reason = n.Reason }).ToList(),
                covariate_centres = model.Covariates.ToDictionary(
                    c => c.Name, c => design.CovariateCentre(c.Name), StringComparer.Ordinal),
                matrix = Rows(design.Matrix),
                fitted_matrix = Rows(design.FittedMatrix),
                design_file = designFile,
                contrasts_file = contrastsFile,
            });
        }

        return new
        {
            samples_file = Path.GetFullPath(samplesPath),
            sample_count = samples.Count,
            factors = model.Factors.Select(f => new { name = f.Name, reference_level = f.ReferenceLevel }).ToList(),
            covariates = model.CovariateEcho,
            interactions = model.Interactions.Select(p => new { first = p.First, second = p.Second }).ToList(),
            stratum_columns = model.Strata,
            batch_column = model.Batch,
            unused_columns = unused,
            stratum_count = designs.Count,
            strata,
            estimability_counts = estimability,
            written = outDirectory is null
                ? null
                : new { directory = Path.GetFullPath(outDirectory), file_count = filesWritten },
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = DesignCaveats(estimability, designs),
        };
    }

    private static Model ReadModel(List<string> lines)
    {
        var model = new Model();
        var interactionsSeen = new HashSet<(string, string)>();
        for (int i = 0; i < lines.Count; i++)
        {
            string[] f = lines[i].TrimEnd('\r').Split('\t').Select(x => x.Trim()).ToArray();
            string where = $"Line {i + 1} of stdin";
            string kind = f[0];
            switch (kind)
            {
                case "factor":
                    Fields(f, 2, 3, where, "factor NAME [REFERENCE]");
                    Claim(model, f[1], "factor", where);
                    model.Factors.Add(f.Length == 3 && f[2].Length > 0 ? new FactorSpec(f[1], f[2]) : new FactorSpec(f[1]));
                    break;
                case "covariate":
                {
                    Fields(f, 3, 6, where, "covariate NAME UNIT [SCALE [SCALED_UNIT [CENTRE]]]");
                    Claim(model, f[1], "covariate", where);
                    if (f[2].Length == 0)
                        throw new Program.UsageException($"{where}: covariate '{f[1]}' needs its stored unit, e.g. years.");
                    double? scale = OptionalNumber(f, 3, where, "SCALE");
                    string? scaledUnit = f.Length > 4 && f[4].Length > 0 ? f[4] : null;
                    double? centre = OptionalNumber(f, 5, where, "CENTRE");
                    var covariate = new CovariateSpec(f[1], f[2]);
                    if (scale is double s) covariate = covariate with { Scale = s };
                    if (scaledUnit is not null) covariate = covariate with { ScaledUnit = scaledUnit };
                    if (centre is double c) covariate = covariate with { Centre = c };
                    model.Covariates.Add(covariate);
                    model.CovariateEcho.Add(new { name = f[1], unit = f[2], scale, scaled_unit = scaledUnit, centre });
                    break;
                }
                case "interaction":
                    Fields(f, 3, 3, where, "interaction FIRST SECOND");
                    if (!interactionsSeen.Add((f[1], f[2])))
                        throw new Program.UsageException($"{where}: the interaction {f[1]}:{f[2]} is given twice.");
                    model.Interactions.Add((f[1], f[2]));
                    break;
                case "stratum":
                    Fields(f, 2, 2, where, "stratum NAME");
                    Claim(model, f[1], "stratum", where);
                    model.Strata.Add(f[1]);
                    break;
                case "batch":
                    Fields(f, 2, 2, where, "batch NAME");
                    if (model.Batch is not null)
                        throw new Program.UsageException($"{where}: a second batch column; a design has one batch.");
                    Claim(model, f[1], "batch", where);
                    model.Batch = f[1];
                    break;
                default:
                    throw new Program.UsageException(
                        $"{where} starts with '{kind}', which is not a model term. Each line is one of: " +
                        string.Join(", ", TermKinds) + ", then its tab-separated fields.");
            }
        }
        return model;
    }

    private static void Fields(string[] fields, int min, int max, string where, string shape)
    {
        if (fields.Length < min || fields.Length > max)
            throw new Program.UsageException(
                $"{where} has {fields.Length} tab-separated fields; a {fields[0]} line is: {shape}.");
        if (fields[1].Length == 0)
            throw new Program.UsageException($"{where}: the {fields[0]} has an empty name.");
    }

    private static void Claim(Model model, string column, string role, string where)
    {
        if (model.Roles.TryGetValue(column, out string? existing))
            throw new Program.UsageException(
                existing == role
                    ? $"{where}: '{column}' is named as a {role} twice."
                    : $"{where}: '{column}' is already the {existing}; a sample sheet column has one role.");
        model.Roles[column] = role;
    }

    private static double? OptionalNumber(string[] fields, int index, string where, string what)
    {
        if (fields.Length <= index || fields[index].Length == 0)
            return null;
        if (!double.TryParse(fields[index], NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
            || !double.IsFinite(value))
            throw new Program.UsageException($"{where}: {what} must be a finite number; got '{fields[index]}'.");
        return value;
    }

    /// <summary>The sheet's rows as mzLib's samples, each carrying only the columns the model names.</summary>
    private static (List<AnalysisSample> Samples, List<string> Unused) ReadSampleSheet(string path, Model model)
    {
        List<string[]> rows = ReadTable(path, "samples");
        string[] header = rows[0];
        if (header.Any(h => h.Length == 0))
            throw new Program.UsageException($"The sample sheet '{path}' has an empty column name.");
        string? repeated = header.GroupBy(h => h, StringComparer.Ordinal).FirstOrDefault(g => g.Count() > 1)?.Key;
        if (repeated is not null)
            throw new Program.UsageException($"The sample sheet's column '{repeated}' appears more than once.");

        var index = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int c = 1; c < header.Length; c++)
            index[header[c]] = c;
        foreach (var (column, role) in model.Roles)
        {
            if (!index.ContainsKey(column))
                throw new Program.UsageException(
                    $"The sample sheet has no column '{column}' (named on stdin as a {role}). Its columns after the sample id: " +
                    (header.Length > 1 ? string.Join(", ", header.Skip(1)) : "none") + ".");
        }

        if (rows.Count == 1)
            throw new Program.UsageException($"The sample sheet '{path}' has a header but no samples.");
        var samples = new List<AnalysisSample>(rows.Count - 1);
        var ids = new HashSet<string>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != header.Length)
                throw new Program.UsageException(
                    $"Line {r + 1} of the sample sheet has {row.Length} fields; the header has {header.Length}.");
            string id = row[0];
            if (id.Length == 0)
                throw new Program.UsageException($"Line {r + 1} of the sample sheet has an empty sample id.");
            // Strata are split before mzLib sees the samples, so a sample named twice in two strata
            // would never meet its twin there.
            if (!ids.Add(id))
                throw new Program.UsageException($"Sample '{id}' is given twice in the sample sheet.");

            var covariates = new Dictionary<string, double>(StringComparer.Ordinal);
            foreach (CovariateSpec covariate in model.Covariates)
            {
                string cell = row[index[covariate.Name]];
                // Blank, NA or NaN is not a value: left out, so mzLib refuses the sample by name.
                if (cell.Length == 0 || cell.Equals("NA", StringComparison.OrdinalIgnoreCase)
                    || cell.Equals("NaN", StringComparison.OrdinalIgnoreCase))
                    continue;
                if (!double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
                    throw new Program.UsageException(
                        $"Line {r + 1} of the sample sheet, covariate '{covariate.Name}': '{cell}' is not a number.");
                covariates[covariate.Name] = value;
            }

            string? batch = model.Batch is null ? null : row[index[model.Batch]];
            samples.Add(new AnalysisSample(id)
            {
                Factors = model.Factors.ToDictionary(f => f.Name, f => row[index[f.Name]], StringComparer.Ordinal),
                Strata = model.Strata.ToDictionary(s => s, s => row[index[s]], StringComparer.Ordinal),
                Batch = string.IsNullOrWhiteSpace(batch) ? null : batch,
                Covariates = covariates,
            });
        }

        var unused = header.Skip(1).Where(h => !model.Roles.ContainsKey(h)).ToList();
        return (samples, unused);
    }

    private static List<string> DesignCaveats(Dictionary<string, int> estimability, IReadOnlyList<AnalysisDesign> designs)
    {
        var caveats = new List<string>
        {
            "Each stratum's matrix is mzLib's AnalysisDesign: an intercept; each factor treatment-coded against its reference level, its other levels in ordinal order; batch against its first level; each covariate as (value - centre) / scale; then any interactions. It equals R's stats::model.matrix for the same formula and level order.",
            "Only mzLib's default contrasts are listed (AnalysisDesign.DefaultContrasts): each level against its factor's reference, every pair of levels when a factor has no reference, and each covariate's slope. With an interaction in the model, a factor's contrast is its effect at the reference level of the other factor.",
            "Contrast ids (c1, c2, ...) are numbered once over every sample of the analysis (AnalysisDesign.ByStratum), so an id names the same comparison in every stratum. A stratum without samples on both sides of a contrast skips its id and lists it in not_run.",
            "A design file holds only the fitted columns (fitted_column_names): mzLib leaves out a column that is redundant with earlier ones, because LinearModel.Fit refuses a design without full column rank. A contrast that is not estimable has no fitted weights and is not in the contrasts file.",
            "Strata are never columns: each stratum gets its own design, centres and fit; only the default contrasts' ids are shared. The individual a sample came from is never a fixed column either; it is for a mixed model's random intercept.",
        };
        int notEstimable = estimability[AnalysisDesign.NotEstimableConfounded] + estimability[AnalysisDesign.NotEstimableRankDeficient];
        if (notEstimable > 0)
            caveats.Add(
                $"{notEstimable} contrast(s) are not estimable ({estimability[AnalysisDesign.NotEstimableConfounded]} confounded, " +
                $"{estimability[AnalysisDesign.NotEstimableRankDeficient]} rank deficient): no fit of these samples can give them an estimate. " +
                "Each stratum's aliased_terms names the terms that cannot be told apart.");
        int notRun = designs.Sum(d => d.NotRun.Count);
        if (notRun > 0)
            caveats.Add($"{notRun} contrast(s) the model asks for cannot run in a stratum; each stratum's not_run says which and why.");
        return caveats;
    }
}
