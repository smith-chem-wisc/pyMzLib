using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json;
using FlashLFQ;
using Quantification.Differential;
using Readers;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>diff design</c> and <c>diff observations</c>. Their answers must be mzLib's
/// <see cref="AnalysisDesign"/>, <see cref="FlashLfqObservations"/> and
/// <see cref="SharedPeptideNormalization"/>, which mzLib holds to base R's <c>model.matrix</c> and to
/// its own stored-table tests. These tests hold the WIRE to the same references, and hold the pivot
/// files to the readers they are written for: <c>stats fit</c> and <c>EmpiricalBayes.ModerateContrasts</c>.
/// </summary>
[TestFixture]
[NonParallelizable] // swaps Console.In
public class DifferentialTests
{
    private string _temp = "";

    [SetUp]
    public void CreateTemp()
    {
        _temp = Path.Combine(Path.GetTempPath(), "mzlib-bridge-diff-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_temp);
    }

    [TearDown]
    public void DeleteTemp()
    {
        if (Directory.Exists(_temp))
            Directory.Delete(_temp, recursive: true);
    }

    // ---- locations and running -------------------------------------------------------------------

    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    private static string MzLibTest(params string[] parts)
    {
        string path = Path.Combine(new[] { RepoRoot(), "code", "mzLib", "mzLib", "Test" }.Concat(parts).ToArray());
        if (!File.Exists(path))
            Assert.Ignore($"mzLib fixture not present in the worktree: {path}");
        return path;
    }

    private static string DesignReference(string name) => MzLibTest("Quantification", "Differential", "ReferenceData", name);

    private static string External(string name) => MzLibTest("FileReadingTests", "ExternalFileTypes", name);

    /// <summary>The inputs the Python guide and docstrings run on.</summary>
    private static string DiffFixture(string name) =>
        Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "diff", name);

    private static string StatsFixture(string name) =>
        Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "stats", name);

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

    private static JsonElement Design(string stdin, params string[] args) =>
        Run(Differential.Design, stdin, new[] { "diff", "design" }.Concat(args).ToArray());

    private static JsonElement Observations(string stdin, params string[] args) =>
        Run(Differential.Observations, stdin, new[] { "diff", "observations" }.Concat(args).ToArray());

    private static List<string?> Strings(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? null : e.GetString()).ToList();

    private static double[][] Matrix(JsonElement element) =>
        element.EnumerateArray().Select(r => r.EnumerateArray().Select(v => v.GetDouble()).ToArray()).ToArray();

    private static List<string> StringList(JsonElement element) =>
        element.EnumerateArray().Select(e => e.GetString()!).ToList();

    private static List<string[]> ReadTsv(string path) =>
        File.ReadAllLines(path).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();

    private static double Parse(string s) => double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);

    // ---- diff design: the matrix is R's model.matrix, through the wire --------------------------------

    /// <summary>
    /// The model of each of mzLib's fourteen reference designs (AnalysisDesignReferenceTests.Specs), as stdin
    /// lines. A sheet's <c>batch</c> column is named as the batch, as mzLib's reader of these sheets does, and its
    /// <c>individual</c> column is left unnamed: it is never a fixed column.
    /// </summary>
    private static readonly Dictionary<string, string[]> ReferenceModels = new()
    {
        ["balanced"] = new[] { "factor\tcondition\tA" },
        ["unbalanced"] = new[] { "factor\tcondition\tA" },
        ["lost_sample"] = new[] { "factor\tcondition\tctrl" },
        ["paired"] = new[] { "factor\ttime\tpre" },
        ["within_between"] = new[] { "factor\tgenotype\twt", "factor\ttime\td0" },
        ["batch"] = new[] { "factor\tcondition\tA" },
        ["age_slope"] = new[] { "covariate\tage\tyears\t10\tdecade\t50" },
        ["age_default_centre"] = new[] { "factor\tsex\tfemale", "covariate\tage\tyears\t10\tdecade" },
        ["two_factors"] = new[] { "factor\tage_group\tyoung", "factor\ttreatment\tvehicle" },
        ["interaction"] = new[] { "factor\tage_group\tyoung", "factor\ttreatment\tvehicle", "interaction\tage_group\ttreatment" },
        ["nonalpha_reference"] = new[] { "factor\tage_group\tyoung" },
        ["confounded"] = new[] { "factor\tcondition\tA" },
        ["no_reference"] = new[] { "factor\tcondition" },
        ["factor_batch_covariate"] = new[] { "factor\tcondition\tctrl", "covariate\tage\tyears\t10\tdecade\t50" },
    };

    private static IEnumerable<string> ReferenceDesigns => ReferenceModels.Keys;

    private static string ReferenceStdin(string name)
    {
        var lines = ReferenceModels[name].ToList();
        if (ReadTsv(DesignReference($"design_{name}_samples.tsv"))[0].Contains("batch"))
            lines.Add("batch\tbatch");
        return string.Join('\n', lines) + "\n";
    }

    [TestCaseSource(nameof(ReferenceDesigns))]
    public void Design_MatrixEqualsRsModelMatrix_ThroughTheWire(string name)
    {
        JsonElement root = Design(ReferenceStdin(name), "--samples", DesignReference($"design_{name}_samples.tsv"));
        List<string[]> expected = ReadTsv(DesignReference($"design_{name}_matrix.tsv"));
        JsonElement stratum = root.GetProperty("strata")[0];
        double[][] matrix = Matrix(stratum.GetProperty("matrix"));

        Assert.That(stratum.GetProperty("stratum").GetString(), Is.EqualTo("all"));
        Assert.That(matrix.Length, Is.EqualTo(expected.Count - 1), "rows");
        Assert.That(matrix[0].Length, Is.EqualTo(expected[0].Length),
            $"columns: R has {string.Join(", ", expected[0])}; the wire {string.Join(", ", StringList(stratum.GetProperty("column_names")))}");
        for (int i = 0; i < matrix.Length; i++)
            for (int j = 0; j < matrix[i].Length; j++)
                Assert.That(matrix[i][j], Is.EqualTo(Parse(expected[i + 1][j])).Within(1e-12), $"{name} [{i}, {j}]");
        int expectedRank = name == "confounded" ? 2 : expected[0].Length;
        Assert.That(stratum.GetProperty("rank").GetInt32(), Is.EqualTo(expectedRank), "R's qr()$rank");
    }

    [TestCaseSource(nameof(ReferenceDesigns))]
    public void Design_EveryFieldIsMzLibsAnalysisDesign(string name)
    {
        JsonElement root = Design(ReferenceStdin(name), "--samples", DesignReference($"design_{name}_samples.tsv"));
        AnalysisDesign design = MzLibDesign(name);
        JsonElement stratum = root.GetProperty("strata")[0];

        Assert.That(StringList(stratum.GetProperty("sample_ids")), Is.EqualTo(design.Samples.Select(s => s.SampleId)));
        Assert.That(StringList(stratum.GetProperty("column_names")), Is.EqualTo(design.ColumnNames));
        Assert.That(StringList(stratum.GetProperty("fitted_column_names")),
            Is.EqualTo(design.FittedColumns.Select(j => design.ColumnNames[j])));
        Assert.That(StringList(stratum.GetProperty("aliased_terms")), Is.EqualTo(design.AliasedTerms));
        Assert.That(StringList(stratum.GetProperty("warnings")), Is.EqualTo(design.Warnings));
        Assert.That(Matrix(stratum.GetProperty("fitted_matrix")), Is.EqualTo(ToJagged(design.FittedMatrix)));

        var contrasts = design.DefaultContrasts();
        Assert.That(Strings(root, "contrast_id"), Is.EqualTo(contrasts.Select(c => c.Id)));
        Assert.That(Strings(root, "label"), Is.EqualTo(contrasts.Select(c => c.Label)));
        Assert.That(Strings(root, "numerator"), Is.EqualTo(contrasts.Select(c => c.Numerator)));
        Assert.That(Strings(root, "covariate_scale"), Is.EqualTo(contrasts.Select(c => c.CovariateScale)));
        Assert.That(Strings(root, "estimability"), Is.EqualTo(contrasts.Select(design.Estimability)));
        var fitted = root.GetProperty("columns").GetProperty("fitted_weights").EnumerateArray().ToList();
        for (int k = 0; k < contrasts.Count; k++)
        {
            if (design.Estimability(contrasts[k]) == AnalysisDesign.Estimable)
                Assert.That(fitted[k].EnumerateArray().Select(v => v.GetDouble()), Is.EqualTo(design.FittedWeights(contrasts[k])));
            else
                Assert.That(fitted[k].ValueKind, Is.EqualTo(JsonValueKind.Null), "no fitted weights when not estimable");
        }
    }

    /// <summary>mzLib's own reading of a reference sheet (AnalysisDesignReferenceTests.ReadSamples), with its spec.</summary>
    private static AnalysisDesign MzLibDesign(string name)
    {
        List<string[]> lines = ReadTsv(DesignReference($"design_{name}_samples.tsv"));
        string[] header = lines[0];
        var samples = lines.Skip(1).Select(row =>
        {
            var factors = new Dictionary<string, string>();
            var covariates = new Dictionary<string, double>();
            string? batch = null;
            for (int c = 1; c < header.Length; c++)
            {
                switch (header[c])
                {
                    case "individual": break;
                    case "batch": batch = row[c]; break;
                    case "age": covariates["age"] = Parse(row[c]); break;
                    default: factors[header[c]] = row[c]; break;
                }
            }
            return new AnalysisSample(row[0]) { Factors = factors, Covariates = covariates, Batch = batch };
        }).ToList();

        var factorSpecs = new List<FactorSpec>();
        var covariateSpecs = new List<CovariateSpec>();
        var interactions = new List<(string, string)>();
        foreach (string[] f in ReferenceModels[name].Select(l => l.Split('\t')))
        {
            if (f[0] == "factor") factorSpecs.Add(new FactorSpec(f[1], f.Length > 2 ? f[2] : null));
            if (f[0] == "interaction") interactions.Add((f[1], f[2]));
            if (f[0] == "covariate")
                covariateSpecs.Add(new CovariateSpec(f[1], f[2], Parse(f[3]), f[4], f.Length > 5 ? Parse(f[5]) : null));
        }
        return AnalysisDesign.Create(samples,
            new DesignSpec { Factors = factorSpecs, Covariates = covariateSpecs, Interactions = interactions });
    }

    private static double[][] ToJagged(double[,] m) =>
        Enumerable.Range(0, m.GetLength(0)).Select(i => Enumerable.Range(0, m.GetLength(1)).Select(j => m[i, j]).ToArray()).ToArray();

    // ---- diff design: what the stdin fields leave to mzLib ----------------------------------------------

    [Test]
    public void Design_AnOmittedCovariateFieldIsMzLibsDefault()
    {
        // No SCALE, SCALED_UNIT or CENTRE: CovariateSpec's own defaults (1, none, the stratum mean).
        JsonElement root = Design("covariate\tage\tyears\n", "--samples", DesignReference("design_age_slope_samples.tsv"));
        List<string[]> sheet = ReadTsv(DesignReference("design_age_slope_samples.tsv"));
        var samples = sheet.Skip(1).Select(r => new AnalysisSample(r[0])
            { Covariates = new Dictionary<string, double> { ["age"] = Parse(r[1]) } }).ToList();
        AnalysisDesign mzLib = AnalysisDesign.Create(samples,
            new DesignSpec { Covariates = new[] { new CovariateSpec("age", "years") } });
        JsonElement stratum = root.GetProperty("strata")[0];

        Assert.That(StringList(stratum.GetProperty("column_names")), Is.EqualTo(mzLib.ColumnNames));
        Assert.That(StringList(stratum.GetProperty("column_names"))[1], Is.EqualTo("age (per years)"));
        Assert.That(stratum.GetProperty("covariate_centres").GetProperty("age").GetDouble(),
            Is.EqualTo(mzLib.CovariateCentre("age")));
        Assert.That(Strings(root, "covariate_unit"), Is.EqualTo(new[] { "years" }));
        JsonElement echo = root.GetProperty("covariates")[0];
        Assert.That(echo.GetProperty("scale").ValueKind, Is.EqualTo(JsonValueKind.Null), "not given, so not sent");
        Assert.That(echo.GetProperty("centre").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void Design_ANoReferenceFactorComparesEveryPairAndWarns()
    {
        JsonElement root = Design("factor\tcondition\n", "--samples", DesignReference("design_no_reference_samples.tsv"));
        AnalysisDesign mzLib = MzLibDesign("no_reference");
        Assert.That(Strings(root, "label"), Is.EqualTo(mzLib.DefaultContrasts().Select(c => c.Label)));
        Assert.That(Strings(root, "label").Count, Is.GreaterThan(1), "every pair of levels");
        Assert.That(StringList(root.GetProperty("strata")[0].GetProperty("warnings")), Has.Some.Contains("GR-14"));
        Assert.That(root.GetProperty("factors")[0].GetProperty("reference_level").ValueKind, Is.EqualTo(JsonValueKind.Null));
    }

    [Test]
    public void Design_BatchEntersOnlyWhenNamed()
    {
        string sheet = DesignReference("design_batch_samples.tsv");
        JsonElement without = Design("factor\tcondition\tA\n", "--samples", sheet);
        JsonElement with = Design("factor\tcondition\tA\nbatch\tbatch\n", "--samples", sheet);

        Assert.That(StringList(without.GetProperty("unused_columns")), Is.EqualTo(new[] { "batch" }));
        Assert.That(without.GetProperty("batch_column").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(StringList(without.GetProperty("strata")[0].GetProperty("column_names")), Has.None.StartsWith("batch="));
        Assert.That(StringList(with.GetProperty("strata")[0].GetProperty("column_names")), Has.Some.StartsWith("batch="));
        Assert.That(with.GetProperty("unused_columns").GetArrayLength(), Is.EqualTo(0));

        // One batch in the stratum: mzLib drops it, and says so.
        string single = Path.Combine(_temp, "one_batch.tsv");
        File.WriteAllText(single, "sample\tcondition\tbatch\ns1\tA\tb1\ns2\tB\tb1\ns3\tA\tb1\ns4\tB\tb1\n");
        JsonElement one = Design("factor\tcondition\tA\nbatch\tbatch\n", "--samples", single);
        Assert.That(StringList(one.GetProperty("strata")[0].GetProperty("column_names")), Is.EqualTo(new[] { "intercept", "condition=B" }));
        Assert.That(StringList(one.GetProperty("strata")[0].GetProperty("warnings")), Has.Some.Contains("batch 'b1'"));
    }

    [Test]
    public void Design_OutOverwritesItsOwnFiles()
    {
        string sheet = DiffFixture("malat_dilution_samples.tsv");
        Design("factor\tmalat_load\t500ng\n", "--samples", sheet, "--out", _temp);
        JsonElement again = Design("factor\tmalat_load\t250ng\n", "--samples", sheet, "--out", _temp);
        Assert.That(again.GetProperty("written").GetProperty("file_count").GetInt32(), Is.EqualTo(2));
        Assert.That(ReadTsv(Path.Combine(_temp, "design_1.tsv"))[0],
            Is.EqualTo(new[] { "sample", "intercept", "malat_load=125ng", "malat_load=500ng" }), "the second call's design");
    }

    [Test]
    public void Design_ConfoundedConditionIsNotEstimableAndHasNoContrastsFile()
    {
        JsonElement root = Design("factor\tcondition\tA\nbatch\tbatch\n", "--samples",
            DesignReference("design_confounded_samples.tsv"), "--out", _temp);
        JsonElement stratum = root.GetProperty("strata")[0];

        Assert.That(Strings(root, "estimability"), Is.EqualTo(new[] { AnalysisDesign.NotEstimableConfounded }));
        Assert.That(StringList(stratum.GetProperty("aliased_terms")), Is.EqualTo(new[] { "condition", "batch" }));
        Assert.That(stratum.GetProperty("contrasts_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(File.Exists(Path.Combine(_temp, "contrasts_1.tsv")), Is.False);
        Assert.That(ReadTsv(Path.Combine(_temp, "design_1.tsv"))[0], Is.EqualTo(new[] { "sample", "intercept", "condition=B" }),
            "the redundant batch column is left out of the fitted design");
        Assert.That(root.GetProperty("estimability_counts").GetProperty(AnalysisDesign.NotEstimableConfounded).GetInt32(), Is.EqualTo(1));
        Assert.That(StringList(root.GetProperty("caveats")), Has.Some.Contains("not estimable"));
    }

    [Test]
    public void Design_StrataAreDesignedApartAndNumberedInMzLibsOrder()
    {
        JsonElement root = Design("stratum\torganism part\nfactor\tcondition\tA\n", "--samples",
            DiffFixture("stratified_samples.tsv"), "--out", _temp);
        var strata = root.GetProperty("strata").EnumerateArray().ToList();

        Assert.That(strata.Select(s => s.GetProperty("stratum").GetString()),
            Is.EqualTo(new[] { "organism part=kidney", "organism part=liver", "organism part=muscle" }), "ordinal by key");
        Assert.That(Strings(root, "stratum"), Is.EqualTo(new[] { "organism part=liver", "organism part=muscle" }));
        Assert.That(Strings(root, "contrast_id"), Is.EqualTo(new[] { "c1", "c1" }), "ids are numbered within each stratum");

        JsonElement kidney = strata[0];
        Assert.That(kidney.GetProperty("not_run")[0].GetProperty("reason").GetString(), Is.EqualTo(AnalysisDesign.NotInStratum));
        Assert.That(kidney.GetProperty("contrasts_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Path.GetFileName(kidney.GetProperty("design_file").GetString()), Is.EqualTo("design_1.tsv"));
        Assert.That(Path.GetFileName(strata[2].GetProperty("contrasts_file").GetString()), Is.EqualTo("contrasts_3.tsv"));
        Assert.That(ReadTsv(Path.Combine(_temp, "design_1.tsv")).Skip(1).Select(r => r[0]), Is.EqualTo(new[] { "s7", "s8" }));
        Assert.That(root.GetProperty("written").GetProperty("file_count").GetInt32(), Is.EqualTo(5));
        Assert.That(StringList(root.GetProperty("caveats")), Has.Some.Contains("not_run"));
    }

    [Test]
    public void Design_AnEmptyInteractionCellLeavesItsColumnOutOfTheFittedDesign()
    {
        // mzLib's AnEmptyInteractionCellIsRankDeficientButMainEffectsStayEstimable: nobody is old AND on drug.
        string sheet = Path.Combine(_temp, "samples.tsv");
        File.WriteAllText(sheet, "sample\tage_group\ttreatment\ns1\tyoung\tvehicle\ns2\tyoung\tdrug\ns3\told\tvehicle\n" +
                                 "s4\tyoung\tvehicle\ns5\tyoung\tdrug\ns6\told\tvehicle\n");
        JsonElement root = Design("factor\tage_group\tyoung\nfactor\ttreatment\tvehicle\ninteraction\tage_group\ttreatment\n",
            "--samples", sheet);
        JsonElement stratum = root.GetProperty("strata")[0];

        Assert.That(StringList(stratum.GetProperty("column_names")), Has.Member("age_group=old:treatment=drug"));
        Assert.That(StringList(stratum.GetProperty("fitted_column_names")), Has.No.Member("age_group=old:treatment=drug"));
        Assert.That(Strings(root, "estimability"), Is.All.EqualTo(AnalysisDesign.Estimable));
    }

    // ---- diff design: the pivot files are what stats fit and ModerateContrasts read ---------------------

    /// <summary>
    /// The MALAT1 dilution series, described by a sample sheet instead of a hand-made design: diff design's file,
    /// tested by stats fit, gives the numbers stats fit gives on the committed 0/1 design, coefficient by
    /// coefficient. It is the same model with its columns in another order, so the fit agrees to rounding.
    /// </summary>
    [Test]
    public void Design_TheDesignFileIsWhatStatsFitReads()
    {
        JsonElement design = Design("factor\tmalat_load\t500ng\n", "--samples", DiffFixture("malat_dilution_samples.tsv"),
            "--out", _temp);
        string designFile = design.GetProperty("strata")[0].GetProperty("design_file").GetString()!;
        string responses = StatsFixture("malat_dilution_log2_vs_fluc.tsv");

        JsonElement ours = Run(Statistics.Fit, "malat_load=250ng\nmalat_load=125ng\n",
            "stats", "fit", "--responses", responses, "--design", designFile);
        JsonElement committed = Run(Statistics.Fit, "malat_250ng\nmalat_125ng\n",
            "stats", "fit", "--responses", responses, "--design", StatsFixture("malat_dilution_design.tsv"));

        Assert.That(Strings(ours, "status"), Is.EqualTo(Strings(committed, "status")));
        foreach (string column in new[] { "estimate", "standard_error", "t", "p_value", "bh_adjusted" })
        {
            var a = ours.GetProperty("columns").GetProperty(column).EnumerateArray().ToList();
            var b = committed.GetProperty("columns").GetProperty(column).EnumerateArray().ToList();
            Assert.That(a.Count, Is.EqualTo(b.Count), column);
            for (int k = 0; k < a.Count; k++)
            {
                Assert.That(a[k].ValueKind, Is.EqualTo(b[k].ValueKind), $"{column} row {k}");
                if (b[k].ValueKind == JsonValueKind.Number)
                    Assert.That(a[k].GetDouble(), Is.EqualTo(b[k].GetDouble()).Within(1e-9 * Math.Max(1, Math.Abs(b[k].GetDouble()))),
                        $"{column} row {k}");
            }
        }
    }

    /// <summary>
    /// The contrasts file is the ContrastWeights EmpiricalBayes.ModerateContrasts takes, on the design file's columns:
    /// each row's moderated contrast equals the moderated coefficient it selects (mzLib #1452 pins that a
    /// one-coefficient contrast gives Moderate's numbers).
    /// </summary>
    [Test]
    public void Design_TheContrastsFileIsWhatModerateContrastsTakes()
    {
        Design("factor\tmalat_load\t500ng\n", "--samples", DiffFixture("malat_dilution_samples.tsv"), "--out", _temp);
        List<string[]> designRows = ReadTsv(Path.Combine(_temp, "design_1.tsv"));
        List<string[]> contrastRows = ReadTsv(Path.Combine(_temp, "contrasts_1.tsv"));
        Assert.That(contrastRows[0].Skip(1), Is.EqualTo(designRows[0].Skip(1)), "the same columns, in the same order");

        // The responses in design row order.
        List<string[]> table = ReadTsv(StatsFixture("malat_dilution_log2_vs_fluc.tsv"));
        string[] samples = table[0].Skip(1).ToArray();
        var order = designRows.Skip(1).Select(r => Array.IndexOf(samples, r[0])).ToArray();
        var y = new double[table.Count - 1, order.Length];
        for (int f = 0; f < table.Count - 1; f++)
            for (int s = 0; s < order.Length; s++)
            {
                string cell = table[f + 1][order[s] + 1];
                y[f, s] = cell.Length == 0 ? double.NaN : Parse(cell);
            }
        var x = new double[order.Length, designRows[0].Length - 1];
        for (int s = 0; s < order.Length; s++)
            for (int j = 0; j < x.GetLength(1); j++)
                x[s, j] = Parse(designRows[s + 1][j + 1]);
        string[] coefficients = designRows[0].Skip(1).ToArray();

        LinearModelFit fit = LinearModel.Fit(y, x, coefficients, 1);
        var weights = contrastRows.Skip(1).Select(r => new ContrastWeights(r[0], r.Skip(1).Select(Parse).ToArray())).ToList();
        IReadOnlyList<ModeratedContrast> moderated = EmpiricalBayes.ModerateContrasts(fit, weights, trend: false);

        Assert.That(moderated.Select(m => m.Name), Is.EqualTo(new[] { "c1", "c2" }));
        foreach (ModeratedContrast contrast in moderated)
        {
            string coefficient = coefficients[Array.IndexOf(contrast.Weights.ToArray(), 1.0)];
            ModeratedTest test = EmpiricalBayes.Moderate(fit, coefficient, trend: false);
            for (int f = 0; f < fit.FeatureCount; f++)
            {
                if (fit.Status[f] != FeatureFitStatus.Fitted) continue;
                Assert.That(contrast.Estimate[f], Is.EqualTo(test.Estimate[f]).Within(1e-12), $"{contrast.Name} feature {f}");
                Assert.That(contrast.T[f], Is.EqualTo(test.T[f]).Within(1e-9), $"{contrast.Name} feature {f}");
            }
        }
    }

    /// <summary>
    /// The whole chain over the wire: diff design --out, then stats fit --design design_1.tsv --contrasts
    /// contrasts_1.tsv, with nothing on stdin. With no reference level every pair of loads is a contrast, and
    /// one of them (500 ng vs 250 ng) is a difference of two coefficients. Each contrast row must be mzLib's
    /// EmpiricalBayes.ModerateContrasts on the same files, called directly.
    /// </summary>
    [Test]
    public void Design_ThenStatsFitContrasts_IsModerateContrastsOnThoseFiles()
    {
        JsonElement design = Design("factor\tmalat_load\n", "--samples", DiffFixture("malat_dilution_samples.tsv"), "--out", _temp);
        JsonElement stratum = design.GetProperty("strata")[0];
        string designFile = stratum.GetProperty("design_file").GetString()!;
        string contrastsFile = stratum.GetProperty("contrasts_file").GetString()!;
        string responses = StatsFixture("malat_dilution_log2_vs_fluc.tsv");

        JsonElement fit = Run(Statistics.Fit, "", "stats", "fit", "--responses", responses, "--design", designFile,
            "--contrasts", contrastsFile);

        JsonElement contrasts = fit.GetProperty("contrasts");
        Assert.That(StringList(contrasts.GetProperty("tested_contrasts")), Is.EqualTo(new[] { "c1", "c2", "c3" }));
        Assert.That(Strings(design, "label"), Is.EqualTo(new[]
        {
            "malat_load=250ng vs malat_load=125ng", "malat_load=500ng vs malat_load=125ng", "malat_load=500ng vs malat_load=250ng",
        }), "the labels a reader joins on contrast_id");
        Assert.That(contrasts.GetProperty("weights").GetProperty("c3").EnumerateArray().Select(e => e.GetDouble()),
            Is.EqualTo(new[] { 0.0, -1, 1 }), "500 ng minus 250 ng, on intercept, malat_load=250ng, malat_load=500ng");
        Assert.That(fit.GetProperty("row_count").GetInt32(), Is.EqualTo(0));
        Assert.That(contrasts.GetProperty("row_count").GetInt32(), Is.EqualTo(3 * fit.GetProperty("feature_count").GetInt32()));

        // mzLib, directly, on the same three files.
        List<string[]> designRows = ReadTsv(designFile);
        List<string[]> table = ReadTsv(responses);
        string[] samples = table[0].Skip(1).ToArray();
        var bySample = designRows.Skip(1).ToDictionary(r => r[0], r => r.Skip(1).Select(Parse).ToArray());
        string[] coefficients = designRows[0].Skip(1).ToArray();
        var y = new double[table.Count - 1, samples.Length];
        for (int f = 0; f < table.Count - 1; f++)
            for (int s = 0; s < samples.Length; s++)
            {
                string cell = table[f + 1][s + 1];
                y[f, s] = cell.Length == 0 ? double.NaN : Parse(cell);
            }
        var x = new double[samples.Length, coefficients.Length];
        for (int s = 0; s < samples.Length; s++)
            for (int j = 0; j < coefficients.Length; j++)
                x[s, j] = bySample[samples[s]][j];
        var weights = ReadTsv(contrastsFile).Skip(1).Select(r => new ContrastWeights(r[0], r.Skip(1).Select(Parse).ToArray())).ToList();
        IReadOnlyList<ModeratedContrast> direct = EmpiricalBayes.ModerateContrasts(LinearModel.Fit(y, x, coefficients, 1), weights, trend: false);

        JsonElement columns = contrasts.GetProperty("columns");
        List<string?> names = columns.GetProperty("contrast").EnumerateArray().Select(e => e.GetString()).ToList();
        foreach (ModeratedContrast contrast in direct)
        {
            foreach (var (column, values) in new[]
            {
                ("estimate", contrast.Estimate), ("standard_error", contrast.StandardError), ("t", contrast.T),
                ("p_value", contrast.PValue), ("bh_adjusted", contrast.BenjaminiHochbergAdjusted),
                ("confidence_low", contrast.ConfidenceLow), ("confidence_high", contrast.ConfidenceHigh),
            })
            {
                double[] wire = columns.GetProperty(column).EnumerateArray()
                    .Where((_, i) => names[i] == contrast.Name)
                    .Select(e => e.ValueKind == JsonValueKind.Null ? double.NaN : e.GetDouble()).ToArray();
                Assert.That(wire, Is.EqualTo(values.ToArray()), $"{contrast.Name} {column}");
            }
        }
    }

    // ---- diff design: refusals ----------------------------------------------------------------------

    [TestCase("factor\tcondition\tA\nfactor\tcondition\n", "named as a factor twice")]
    [TestCase("factor\tcondition\tA\nstratum\tcondition\n", "has one role")]
    [TestCase("level\tcondition\n", "not a model term")]
    [TestCase("factor\tnope\n", "no column 'nope'")]
    [TestCase("factor\tcondition\tA\textra\n", "a factor line is")]
    [TestCase("covariate\tage\n", "a covariate line is")]
    [TestCase("covariate\tcondition\tyears\tten\n", "SCALE must be a finite number")]
    [TestCase("batch\tcondition\nbatch\tcondition\n", "a second batch")]
    [TestCase("factor\tcondition\tA\ninteraction\tcondition\tsex\n", "names a factor that is not in the design")]
    public void Design_RefusesAMalformedModel(string stdin, string message)
    {
        string text = Usage(Differential.Design, stdin, "diff", "design", "--samples", DesignReference("design_balanced_samples.tsv"));
        Assert.That(text, Does.Contain(message));
    }

    [Test]
    public void Design_RefusesAMalformedSampleSheet()
    {
        string sheet = Path.Combine(_temp, "sheet.tsv");
        File.WriteAllText(sheet, "sample\ttissue\tcondition\ns1\tliver\tA\ns1\tmuscle\tB\n");
        Assert.That(Usage(Differential.Design, "stratum\ttissue\nfactor\tcondition\tA\n", "diff", "design", "--samples", sheet),
            Does.Contain("Sample 's1' is given twice"), "a duplicate across strata, which mzLib never sees together");

        File.WriteAllText(sheet, "sample\tcondition\tage\ns1\tA\tforty\ns2\tB\t50\n");
        Assert.That(Usage(Differential.Design, "covariate\tage\tyears\n", "diff", "design", "--samples", sheet),
            Does.Contain("'forty' is not a number"));

        File.WriteAllText(sheet, "sample\tcondition\ns1\tA\ns2\t\n");
        Assert.That(Usage(Differential.Design, "factor\tcondition\tA\n", "diff", "design", "--samples", sheet),
            Does.Contain("has no value for factor 'condition'"), "mzLib's own refusal, as a usage error");

        File.WriteAllText(sheet, "sample\tcondition\tage\ns1\tA\tNA\ns2\tB\t50\n");
        Assert.That(Usage(Differential.Design, "covariate\tage\tyears\n", "diff", "design", "--samples", sheet),
            Does.Contain("needs a finite value for covariate 'age'"), "NA is no value, so mzLib refuses the sample");

        Assert.That(Usage(Differential.Design, "", "diff", "design", "--samples", Path.Combine(_temp, "absent.tsv")),
            Does.Contain("not found"));
        Assert.That(Usage(Differential.Design, "", "diff", "design", "--samples", sheet, "--out"),
            Does.Contain("--out needs a directory"));
    }

    // ---- diff observations: the table is mzLib's, through the wire --------------------------------------

    private static string MetaMorpheusRuns()
    {
        var stored = new QuantifiedPeptideFile(External("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"));
        return string.Join('\n', stored.First().Samples.Keys
            .Select(l => l.EndsWith("-calib", StringComparison.Ordinal) ? l[..^"-calib".Length] : l)
            .Select(n => $"{n}\t{n}")) + "\n";
    }

    [TestCase("msms_only", QuantBasis.MsmsOnly)]
    [TestCase("mbr_kept", QuantBasis.MbrKept)]
    public void Observations_AreMzLibsTableBitForBit(string wire, QuantBasis basis)
    {
        string peptides = External("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv");
        string groups = External("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv");
        JsonElement root = Observations(MetaMorpheusRuns(), "--peptides", peptides, "--protein-groups", groups, "--basis", wire);

        var stored = new QuantifiedPeptideFile(peptides);
        var runs = stored.First().Samples.Keys.Select(l => l[..^"-calib".Length]).Select(n => new ObservationRun(n, n)).ToArray();
        ObservationTable table = FlashLfqObservations.FromPeptideTable(stored, runs,
            FlashLfqObservations.ReadProteinGroups(new ProteinGroupFromTsvFile(groups)), basis);
        var cells = table.Observations().ToList();

        Assert.That(Strings(root, "full_sequence"), Is.EqualTo(cells.Select(c => c.Peptide.FullSequence)));
        Assert.That(Strings(root, "sample"), Is.EqualTo(cells.Select(c => c.SampleId)));
        Assert.That(Strings(root, "state"), Is.EqualTo(cells.Select(c => Differential.StateName(c.State))));
        Assert.That(Strings(root, "unique_protein_group"), Is.EqualTo(cells.Select(c => c.Peptide.UniqueProteinGroup)));
        var log2 = root.GetProperty("columns").GetProperty("log2_intensity").EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? double.NaN : e.GetDouble()).ToList();
        Assert.That(log2.Select(BitConverter.DoubleToInt64Bits), Is.EqualTo(cells.Select(c => BitConverter.DoubleToInt64Bits(c.Log2Intensity))));
        Assert.That(StringList(root.GetProperty("samples_without_values")), Is.EqualTo(table.SamplesWithoutValues));
        Assert.That(StringList(root.GetProperty("warnings")), Is.EqualTo(table.Warnings));
        Assert.That(root.GetProperty("basis").GetString(), Is.EqualTo(wire));

        // mzLib's ReadsAMetaMorpheusPeptideTable: CACASHVAK's transfer in QE-002118_GM6_a.
        int i = cells.FindIndex(c => c.Peptide.BaseSequence == "CACASHVAK" && c.SampleId == "QE-002118_GM6_a");
        Assert.That(Strings(root, "state")[i], Is.EqualTo("mbr_transferred"));
        if (basis == QuantBasis.MbrKept)
            Assert.That(log2[i], Is.EqualTo(Math.Log2(143944.34375)));
        else
            Assert.That(log2[i], Is.NaN, "msms_only: the transferred value enters no sum, and the cell still says why");
    }

    [Test]
    public void Observations_CarryTheProteinGroupsFlags()
    {
        JsonElement root = Observations(MetaMorpheusRuns(),
            "--peptides", External("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"),
            "--protein-groups", External("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"), "--basis", "mbr_kept");
        var groups = root.GetProperty("protein_groups").EnumerateArray().ToList();
        Assert.That(groups.Where(g => g.GetProperty("is_contaminant").GetBoolean()).Select(g => g.GetProperty("name").GetString()),
            Is.EqualTo(new[] { "P02769" }));
        Assert.That(groups.Where(g => g.GetProperty("is_decoy").GetBoolean()).Select(g => g.GetProperty("name").GetString()),
            Is.EqualTo(new[] { "DECOY_P62750" }));
        Assert.That(root.GetProperty("protein_group_count").GetInt32(), Is.EqualTo(6));
        Assert.That(StringList(root.GetProperty("warnings")).Single(), Does.Contain("Q16643"),
            "the peptides' other groups are not in the protein-group table");
    }

    /// <summary>
    /// Fractions of one technical replicate are summed on the linear scale, then the replicates of a sample are
    /// averaged in log2: the ObservationTable's own arithmetic, through the stdin runs. A fraction or technical
    /// replicate left off a line is ObservationRun's default, 0.
    /// </summary>
    [Test]
    public void Observations_SumFractionsThenAverageReplicates()
    {
        string peptides = Path.Combine(_temp, "QuantifiedPeptides.tsv");
        File.WriteAllLines(peptides, new[]
        {
            "Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism\tIntensity_a\tIntensity_b\tIntensity_c\tIntensity_d\t" +
            "Detection Type_a\tDetection Type_b\tDetection Type_c\tDetection Type_d",
            "PEPA\tPEPA\tP1\t\t\t100\t300\t1000\t50\tMSMS\tMSMS\tMSMS\tMBR",
            "PEPB\tPEPB\tP1;P2\t\t\t0\t0\t0\t0\tNotDetected\tMSMSIdentifiedButNotQuantified\tNotDetected\tNotDetected",
        });
        // S1: replicate 0 = fractions a + b; replicate 1 = c. S2: d alone (a transfer).
        JsonElement root = Observations("a\tS1\nb\tS1\t1\nc\tS1\t0\t1\nd\tS2\n", "--peptides", peptides, "--basis", "msms_only");

        var log2 = root.GetProperty("columns").GetProperty("log2_intensity").EnumerateArray().ToList();
        Assert.That(log2[0].GetDouble(), Is.EqualTo((Math.Log2(100 + 300) + Math.Log2(1000)) / 2).Within(1e-12));
        Assert.That(log2[1].ValueKind, Is.EqualTo(JsonValueKind.Null), "a transfer, left out under msms_only");
        // PEPB in S1: no value in any run, so the most informative of its runs' states.
        Assert.That(log2[2].ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Strings(root, "state"), Is.EqualTo(new[] { "quantified", "mbr_transferred", "identified_not_quantified", "not_detected" }));
        Assert.That(Strings(root, "unique_protein_group"), Is.EqualTo(new[] { "P1", "P1", null, null }), "PEPB is shared");
        var runs = root.GetProperty("runs").EnumerateArray().ToList();
        Assert.That(runs.Select(r => (r.GetProperty("file_name").GetString(), r.GetProperty("fraction").GetInt32(),
                r.GetProperty("technical_replicate").GetInt32())),
            Is.EqualTo(new[] { ("a", 0, 0), ("b", 1, 0), ("c", 0, 1), ("d", 0, 0) }));
        Assert.That(StringList(root.GetProperty("caveats")), Has.Some.Contains("No --protein-groups"));
    }

    // ---- diff observations: normalization -----------------------------------------------------------

    private static string Vignette => DiffFixture("EcoliSpikeIn_QuantifiedPeptides.tsv");

    private const string VignetteRuns = "09-04-18_EcoliSpikeInSingleShot1x\tspike_1x\n09-04-18_EcoliSpikeInSingleShot2x\tspike_2x\n";

    private static ObservationTable VignetteTable() =>
        FlashLfqObservations.FromPeptideTable(new QuantifiedPeptideFile(Vignette),
            new[] { new ObservationRun("09-04-18_EcoliSpikeInSingleShot1x", "spike_1x"), new ObservationRun("09-04-18_EcoliSpikeInSingleShot2x", "spike_2x") },
            Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly);

    [Test]
    public void Observations_NormalizationIsMeasuredWithMzLibsFloorWhenNoneIsGiven()
    {
        JsonElement root = Observations(VignetteRuns, "--peptides", Vignette, "--basis", "msms_only");
        NormalizationResult mzLib = SharedPeptideNormalization.Measure(VignetteTable());
        JsonElement n = root.GetProperty("normalization");

        Assert.That(n.GetProperty("setting").GetString(), Is.EqualTo(mzLib.Setting).And.EqualTo(SharedPeptideNormalization.None));
        Assert.That(n.GetProperty("minimum_reference_set_size").GetInt32(), Is.EqualTo(SharedPeptideNormalization.DefaultMinimumPeptides));
        Assert.That(n.GetProperty("reference_set_size").GetInt32(), Is.EqualTo(mzLib.ReferenceSetSize));
        Assert.That(n.GetProperty("shifts").EnumerateObject().Select(p => p.Value.ValueKind), Is.All.EqualTo(JsonValueKind.Null));
        Assert.That(StringList(n.GetProperty("warnings")), Is.EqualTo(mzLib.Warnings));
        Assert.That(root.GetProperty("min_peptides").ValueKind, Is.EqualTo(JsonValueKind.Null), "not given, so not sent");
    }

    [Test]
    public void Observations_NormalizeSubtractsMzLibsShifts()
    {
        string responses = Path.Combine(_temp, "responses.tsv");
        JsonElement plain = Observations(VignetteRuns, "--peptides", Vignette, "--basis", "msms_only", "--min-peptides", "50");
        JsonElement shifted = Observations(VignetteRuns, "--peptides", Vignette, "--basis", "msms_only", "--min-peptides", "50",
            "--normalize", "--responses-out", responses);
        ObservationTable table = VignetteTable();
        NormalizationResult mzLib = SharedPeptideNormalization.Measure(table, 50);
        ObservationTable applied = SharedPeptideNormalization.Apply(table, mzLib);

        JsonElement shifts = shifted.GetProperty("normalization").GetProperty("shifts");
        Assert.That(shifted.GetProperty("normalization").GetProperty("setting").GetString(), Is.EqualTo(SharedPeptideNormalization.SharedPeptideMedian));
        foreach (string sample in table.Samples)
            Assert.That(shifts.GetProperty(sample).GetDouble(), Is.EqualTo(mzLib.PerSampleShift[sample]));
        Assert.That(plain.GetProperty("normalization").GetProperty("shifts").GetProperty("spike_2x").GetDouble(),
            Is.EqualTo(mzLib.PerSampleShift["spike_2x"]), "measured without --normalize too");

        var expected = applied.Observations().Select(o => BitConverter.DoubleToInt64Bits(o.Log2Intensity)).ToList();
        var actual = shifted.GetProperty("columns").GetProperty("log2_intensity").EnumerateArray()
            .Select(e => BitConverter.DoubleToInt64Bits(e.ValueKind == JsonValueKind.Null ? double.NaN : e.GetDouble())).ToList();
        Assert.That(actual, Is.EqualTo(expected));

        // The responses table is the same values, peptide by sample, blank for no value.
        List<string[]> rows = File.ReadAllLines(responses).Select(l => l.Split('\t')).ToList();
        Assert.That(rows[0], Is.EqualTo(new[] { "full_sequence", "spike_1x", "spike_2x" }));
        Assert.That(rows.Count - 1, Is.EqualTo(applied.Peptides.Count));
        for (int p = 0; p < applied.Peptides.Count; p++)
        {
            Assert.That(rows[p + 1][0], Is.EqualTo(applied.Peptides[p].FullSequence));
            for (int s = 0; s < applied.Samples.Count; s++)
            {
                double v = applied.Log2Intensity(p, s);
                Assert.That(rows[p + 1][s + 1], Is.EqualTo(double.IsNaN(v) ? "" : v.ToString("R", CultureInfo.InvariantCulture)));
            }
        }
        Assert.That(shifted.GetProperty("responses_written").GetProperty("row_count").GetInt32(), Is.EqualTo(applied.Peptides.Count));
    }

    [Test]
    public void Observations_TheResponsesTableIsWhatStatsFitReads()
    {
        string responses = Path.Combine(_temp, "responses.tsv");
        Observations(MetaMorpheusRuns(), "--peptides", External("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"),
            "--basis", "mbr_kept", "--responses-out", responses);
        string[] samples = File.ReadLines(responses).First().Split('\t').Skip(1).ToArray();
        string design = Path.Combine(_temp, "design.tsv");
        File.WriteAllLines(design, new[] { "sample\tintercept" }.Concat(samples.Select(s => s + "\t1")));

        JsonElement fit = Run(Statistics.Fit, "intercept\n", "stats", "fit", "--responses", responses, "--design", design);
        Assert.That(fit.GetProperty("feature_count").GetInt32(), Is.EqualTo(3));
        Assert.That(fit.GetProperty("sample_count").GetInt32(), Is.EqualTo(18));
    }

    [Test]
    public void Observations_SamplesNarrowsTheTableBeforeMeasuring()
    {
        string keep = Path.Combine(_temp, "keep.tsv");
        File.WriteAllText(keep, "sample\nspike_2x\n");
        JsonElement root = Observations(VignetteRuns, "--peptides", Vignette, "--basis", "msms_only", "--samples", keep, "--min-peptides", "1");
        ObservationTable narrowed = VignetteTable().ForSamples(new[] { "spike_2x" });
        NormalizationResult mzLib = SharedPeptideNormalization.Measure(narrowed, 1);

        Assert.That(StringList(root.GetProperty("samples")), Is.EqualTo(new[] { "spike_2x" }));
        Assert.That(root.GetProperty("run_count").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("normalization").GetProperty("setting").GetString(), Is.EqualTo(mzLib.Setting));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(narrowed.Peptides.Count));
    }

    [Test]
    public void Observations_OutWritesTheTableInsteadOfReturningIt()
    {
        string path = Path.Combine(_temp, "observations.tsv");
        JsonElement root = Observations(VignetteRuns, "--peptides", Vignette, "--basis", "msms_only", "--out", path);
        Assert.That(root.GetProperty("columns").GetProperty("state").GetArrayLength(), Is.EqualTo(0));
        Assert.That(root.GetProperty("written").GetProperty("row_count").GetInt32(), Is.EqualTo(486));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(486));
        string[] lines = File.ReadAllLines(path);
        Assert.That(lines[0].Split('\t'), Is.EqualTo(Differential.ObservationColumns));
        Assert.That(lines.Length, Is.EqualTo(487));
    }

    // ---- diff observations: refusals and names -----------------------------------------------------

    [Test]
    public void Observations_EveryStateHasAWireName()
    {
        foreach (ObservationState state in Enum.GetValues<ObservationState>())
            Assert.DoesNotThrow(() => Differential.StateName(state), state.ToString());
        Assert.That(Differential.States, Is.EquivalentTo(Enum.GetValues<ObservationState>()));
    }

    [Test]
    public void Observations_RefuseWhatTheyCannotRead()
    {
        string[] Args(params string[] extra) =>
            new[] { "diff", "observations", "--peptides", Vignette }.Concat(extra).ToArray();

        Assert.That(Usage(Differential.Observations, VignetteRuns, Args()), Does.Contain("--basis"));
        Assert.That(Usage(Differential.Observations, VignetteRuns, Args("--basis", "mbr")), Does.Contain("msms_only, mbr_kept"));
        Assert.That(Usage(Differential.Observations, "", Args("--basis", "msms_only")), Does.Contain("No runs"));
        Assert.That(Usage(Differential.Observations, "only-a-file\n", Args("--basis", "msms_only")), Does.Contain("tab-separated fields"));
        Assert.That(Usage(Differential.Observations, "a\tS1\tone\n", Args("--basis", "msms_only")), Does.Contain("fraction must be an integer"));
        Assert.That(Usage(Differential.Observations, VignetteRuns, Args("--basis", "msms_only", "--min-peptides", "0")),
            Does.Contain("1 or more"));
        Assert.That(Usage(Differential.Observations, VignetteRuns, Args("--basis", "msms_only", "--out", "x.csv")), Does.Contain(".tsv"));
        Assert.That(Usage(Differential.Observations, VignetteRuns, Args("--basis", "msms_only", "--protein-groups")),
            Does.Contain("without a value"));
        // mzLib's refusal: a column the runs do not name.
        Assert.That(Usage(Differential.Observations, "09-04-18_EcoliSpikeInSingleShot1x\tS1\n", Args("--basis", "msms_only")),
            Does.Contain("not in the experimental design"));
        Assert.That(Usage(Differential.Observations, VignetteRuns,
            new[] { "diff", "observations", "--peptides", External("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"), "--basis", "msms_only" }),
            Does.Contain("cannot be read by diff observations"), "a protein-group table given as the peptides");
    }

    // ---- the committed Python inputs are what they say they are ----------------------------------------

    [Test]
    public void PythonInputs_AreMzLibsFilesOrTheirRule()
    {
        static string Text(string path) => File.ReadAllText(path).Replace("\r\n", "\n");

        foreach (string name in new[] { "design_factor_batch_covariate_samples.tsv", "design_factor_batch_covariate_matrix.tsv", "design_confounded_samples.tsv" })
            Assert.That(Text(DiffFixture(name)), Is.EqualTo(Text(DesignReference(name))), name);
        foreach (string name in new[] { "MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv", "MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv" })
            Assert.That(Text(DiffFixture(name)), Is.EqualTo(Text(External(name))), name);
        Assert.That(Text(Vignette), Is.EqualTo(Text(Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "flashlfq_vignette_peptides.tsv"))),
            "the FlashLFQ vignette slice, under a name mzLib reads as a peptide table");

        // The MALAT1 sheet: each run's load, from the committed 0/1 design.
        List<string[]> design = ReadTsv(StatsFixture("malat_dilution_design.tsv"));
        int i250 = Array.IndexOf(design[0], "malat_250ng"), i125 = Array.IndexOf(design[0], "malat_125ng");
        var expected = new[] { "run\tmalat_load" }.Concat(design.Skip(1).Select(r =>
            $"{r[0]}\t{(r[i250] == "1" ? "250ng" : r[i125] == "1" ? "125ng" : "500ng")}"));
        Assert.That(Text(DiffFixture("malat_dilution_samples.tsv")), Is.EqualTo(string.Join("\n", expected) + "\n"));

        // mzLib's SamplesAreSplitByStratum samples, plus two kidney samples of condition A only.
        Assert.That(Text(DiffFixture("stratified_samples.tsv")), Is.EqualTo(
            "sample_id\torganism part\tcondition\ns1\tliver\tA\ns2\tmuscle\tA\ns3\tliver\tB\ns4\tmuscle\tB\n" +
            "s5\tliver\tA\ns6\tmuscle\tB\ns7\tkidney\tA\ns8\tkidney\tA\n"));
    }
}
