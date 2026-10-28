using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Chromatography.RetentionTimeCalibration;
using MzLibUtil;

namespace MzLibBridge.Tests;

/// <summary>
/// <c>library calibrate-rt</c>. Every number must be mzLib's: each test fits the same anchors with
/// <see cref="IrtCalibration.Fit"/> directly and holds the wire to that model, so a parsing slip, a
/// group mixed into another or an option that did not reach mzLib shows as a disagreement.
/// The anchors are real: two DIA-NN reports in mzLib's test data, read by <c>readers read-records</c>
/// and cut to six columns (pkg/python/tests/fixtures/calibration/PROVENANCE.txt).
/// </summary>
[TestFixture]
[NonParallelizable] // swaps Console.In
public class LibraryCalibrateRtTests
{
    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    private static string CalibrationFixture(string name) =>
        Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "calibration", name);

    /// <summary>The rows of a calibration fixture TSV, as name-to-text dictionaries.</summary>
    private static List<Dictionary<string, string>> Rows(string name)
    {
        string[] lines = File.ReadAllLines(CalibrationFixture(name));
        string[] header = lines[0].Split('\t');
        return lines.Skip(1).Where(l => l.Length > 0)
            .Select(l => header.Zip(l.Split('\t')).ToDictionary(p => p.First, p => p.Second))
            .ToList();
    }

    private static double D(string text) => double.Parse(text, CultureInfo.InvariantCulture);

    private static JsonElement Run(string stdin, params string[] options)
    {
        TextReader original = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            object result = Library.CalibrateRt(new Program.Arguments(new[] { "library", "calibrate-rt" }.Concat(options).ToArray()));
            return JsonDocument.Parse(JsonSerializer.Serialize(result, Program.JsonOptions)).RootElement;
        }
        finally
        {
            Console.SetIn(original);
        }
    }

    private static string Usage(string stdin, params string[] options)
    {
        TextReader original = Console.In;
        Console.SetIn(new StringReader(stdin));
        try
        {
            return Assert.Throws<Program.UsageException>(() =>
                Library.CalibrateRt(new Program.Arguments(new[] { "library", "calibrate-rt" }.Concat(options).ToArray())))!.Message;
        }
        finally
        {
            Console.SetIn(original);
        }
    }

    private static double?[] Numbers(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? (double?)null : e.GetDouble()).ToArray();

    private static string?[] Texts(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? null : e.GetString()).ToArray();

    private static bool?[] Flags(JsonElement root, string column) =>
        root.GetProperty("columns").GetProperty(column).EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? (bool?)null : e.GetBoolean()).ToArray();

    private static List<string> Caveats(JsonElement root) =>
        root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();

    /// <summary>The HeLa run: q-value 0.01 or below are anchors; the rest are queried both ways.</summary>
    private static (List<Dictionary<string, string>> Anchors, List<Dictionary<string, string>> Others, string Stdin) HeLa()
    {
        var rows = Rows("diann2_hela_rt.tsv");
        var anchors = rows.Where(r => D(r["q_value"]) <= 0.01).ToList();
        var others = rows.Where(r => D(r["q_value"]) > 0.01).ToList();
        var lines = anchors.Select(r => $"{r["run"]}\t{r["retention_time"]}\t{r["indexed_retention_time"]}")
            .Append("--")
            .Concat(others.Select(r => $"{r["run"]}\t\t{r["indexed_retention_time"]}"))
            .Concat(others.Select(r => $"{r["run"]}\t{r["retention_time"]}\t"));
        return (anchors, others, string.Join('\n', lines) + "\n");
    }

    private static IrtCalibrationModel FitDirectly(IEnumerable<Dictionary<string, string>> anchors, IrtCalibrationOptions? options = null) =>
        IrtCalibration.Fit(anchors.Select(r => (new RtMinutes(D(r["retention_time"])), new Irt(D(r["indexed_retention_time"])))).ToList(), options);

    // ---- agreement with mzLib --------------------------------------------------------------------

    [Test]
    public void ARealRunIsMzLibsFitAndEveryConversionIsMzLibs()
    {
        var (anchors, others, stdin) = HeLa();
        JsonElement root = Run(stdin);
        IrtCalibrationModel model = FitDirectly(anchors);

        JsonElement group = root.GetProperty("groups")[0];
        Assert.That(group.GetProperty("kind").GetString(), Is.EqualTo("lowess"), "mzLib's default kind");
        Assert.That(group.GetProperty("anchor_count").GetInt32(), Is.EqualTo(model.AnchorCount));
        Assert.That(group.GetProperty("residual_sd").GetDouble(), Is.EqualTo(model.ResidualSd));
        Assert.That(group.GetProperty("first_anchor_minutes").GetDouble(), Is.EqualTo(model.FirstAnchor.Value));
        Assert.That(group.GetProperty("last_anchor_minutes").GetDouble(), Is.EqualTo(model.LastAnchor.Value));
        Assert.That(group.GetProperty("query_count").GetInt32(), Is.EqualTo(2 * others.Count));
        Assert.That(group.GetProperty("error").ValueKind, Is.EqualTo(JsonValueKind.Null));

        int n = anchors.Count, m = others.Count;
        string?[] roles = Texts(root, "role");
        double?[] rt = Numbers(root, "rt_minutes"), irt = Numbers(root, "irt"), fitted = Numbers(root, "fitted_irt");
        bool?[] outside = Flags(root, "is_extrapolated");
        Assert.That(roles, Is.EqualTo(Enumerable.Repeat("anchor", n).Concat(Enumerable.Repeat("to_minutes", m)).Concat(Enumerable.Repeat("to_irt", m))));
        for (int i = 0; i < n; i++)
        {
            var at = new RtMinutes(D(anchors[i]["retention_time"]));
            Assert.That(rt[i], Is.EqualTo(at.Value));
            Assert.That(irt[i], Is.EqualTo(D(anchors[i]["indexed_retention_time"])));
            Assert.That(fitted[i], Is.EqualTo(model.ToIrt(at).Value), $"anchor {i}");
            Assert.That(outside[i], Is.EqualTo(model.IsExtrapolated(at)));
        }
        for (int j = 0; j < m; j++)
        {
            var libraryIrt = new Irt(D(others[j]["indexed_retention_time"]));
            RtMinutes predicted = model.ToRtMinutes(libraryIrt);
            Assert.That(rt[n + j], Is.EqualTo(predicted.Value), $"to_minutes {j}");
            Assert.That(irt[n + j], Is.EqualTo(libraryIrt.Value));
            Assert.That(fitted[n + j], Is.Null, "a query carries no fitted_irt");
            Assert.That(outside[n + j], Is.EqualTo(model.IsExtrapolated(predicted)));

            var observed = new RtMinutes(D(others[j]["retention_time"]));
            Assert.That(rt[n + m + j], Is.EqualTo(observed.Value));
            Assert.That(irt[n + m + j], Is.EqualTo(model.ToIrt(observed).Value), $"to_irt {j}");
            Assert.That(outside[n + m + j], Is.EqualTo(model.IsExtrapolated(observed)));
        }

        int extrapolated = outside.Skip(n).Count(x => x == true);
        Assert.That(extrapolated, Is.GreaterThan(0), "this run has queries outside its anchors");
        Assert.That(Caveats(root), Has.Some.StartWith($"{extrapolated} quer"));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(n + 2 * m));
        Assert.That(root.GetProperty("anchor_count").GetInt32(), Is.EqualTo(n));
        Assert.That(root.GetProperty("query_count").GetInt32(), Is.EqualTo(2 * m));
    }

    [Test]
    public void AbsentOptionsAreMzLibsDefaultsAndGivenOnesReachMzLib()
    {
        var (anchors, _, stdin) = HeLa();
        var defaults = new IrtCalibrationOptions();
        JsonElement plain = Run(stdin);
        Assert.That(plain.GetProperty("kind").GetString(), Is.EqualTo(Library.KindName(defaults.Kind)));
        Assert.That(plain.GetProperty("minimum_anchors").GetInt32(), Is.EqualTo(defaults.MinimumAnchors));
        Assert.That(plain.GetProperty("bandwidth").GetDouble(), Is.EqualTo(defaults.Bandwidth));
        Assert.That(plain.GetProperty("knots").GetInt32(), Is.EqualTo(defaults.Knots));
        Assert.That(plain.GetProperty("robustness_iterations").GetInt32(), Is.EqualTo(defaults.RobustnessIterations));
        Assert.That(plain.GetProperty("on_error").GetString(), Is.EqualTo("fail"));

        // Each option alone, against mzLib given the same option: each must change the answer the
        // way mzLib's does, so an option the bridge dropped would fail here.
        var cases = new (string[] Options, IrtCalibrationOptions Expected)[]
        {
            (new[] { "--kind", "linear" }, defaults with { Kind = IrtCalibrationKind.Linear }),
            (new[] { "--bandwidth", "0.1" }, defaults with { Bandwidth = 0.1 }),
            (new[] { "--knots", "7" }, defaults with { Knots = 7 }),
            (new[] { "--robustness-iterations", "0" }, defaults with { RobustnessIterations = 0 }),
            (new[] { "--minimum-anchors", "400" }, defaults with { MinimumAnchors = 400 }),
        };
        double plainSd = plain.GetProperty("groups")[0].GetProperty("residual_sd").GetDouble();
        foreach (var (options, expected) in cases)
        {
            JsonElement root = Run(stdin, options);
            IrtCalibrationModel model = FitDirectly(anchors, expected);
            JsonElement group = root.GetProperty("groups")[0];
            Assert.That(group.GetProperty("residual_sd").GetDouble(), Is.EqualTo(model.ResidualSd), string.Join(' ', options));
            Assert.That(group.GetProperty("kind").GetString(), Is.EqualTo(Library.KindName(model.Kind)));
            Assert.That(root.GetProperty("knots").GetInt32(), Is.EqualTo(expected.Knots));
            Assert.That(root.GetProperty("bandwidth").GetDouble(), Is.EqualTo(expected.Bandwidth));
            Assert.That(root.GetProperty("minimum_anchors").GetInt32(), Is.EqualTo(expected.MinimumAnchors));
            Assert.That(root.GetProperty("robustness_iterations").GetInt32(), Is.EqualTo(expected.RobustnessIterations));
            Assert.That(Numbers(root, "fitted_irt")[0],
                Is.EqualTo(model.ToIrt(new RtMinutes(D(anchors[0]["retention_time"]))).Value));
            if (options[0] != "--minimum-anchors")
                Assert.That(group.GetProperty("residual_sd").GetDouble(), Is.Not.EqualTo(plainSd),
                    $"{options[0]} changed nothing, so the test cannot tell whether it reached mzLib");
        }
        Assert.That(Usage(stdin, "--minimum-anchors", "419"), Does.Contain("needs at least 419 anchors; got 418"),
            "--minimum-anchors reaches mzLib: one more than the anchors is refused");
    }

    [Test]
    public void EachRunIsFittedOnItsOwnAnchorsInOrderOfFirstAppearance()
    {
        var rows = Rows("diann_pbmc_rt.tsv");
        // Interleave the runs, so a bridge that fitted lines in blocks would mix them.
        var interleaved = rows.OrderBy(r => r["precursor_id"]).ThenBy(r => r["run"]).ToList();
        string stdin = string.Join('\n', interleaved.Select(r => $"{r["run"]}\t{r["retention_time"]}\t{r["indexed_retention_time"]}")
            .Append("--").Append($"{interleaved[0]["run"]}\t\t50")) + "\n";
        JsonElement root = Run(stdin, "--kind", "linear", "--minimum-anchors", "10");

        string[] order = interleaved.Select(r => r["run"]).Distinct().ToArray();
        JsonElement[] groups = root.GetProperty("groups").EnumerateArray().ToArray();
        Assert.That(groups.Select(g => g.GetProperty("group").GetString()), Is.EqualTo(order));
        var options = new IrtCalibrationOptions(Kind: IrtCalibrationKind.Linear, MinimumAnchors: 10);
        foreach (JsonElement group in groups)
        {
            string run = group.GetProperty("group").GetString()!;
            IrtCalibrationModel model = FitDirectly(rows.Where(r => r["run"] == run), options);
            Assert.That(group.GetProperty("residual_sd").GetDouble(), Is.EqualTo(model.ResidualSd), run);
            Assert.That(group.GetProperty("anchor_count").GetInt32(), Is.EqualTo(rows.Count(r => r["run"] == run)));
            Assert.That(group.GetProperty("query_count").GetInt32(), Is.EqualTo(run == order[0] ? 1 : 0));
            Assert.That(group.GetProperty("kind").GetString(), Is.EqualTo("linear"));
        }
        IrtCalibrationModel first = FitDirectly(rows.Where(r => r["run"] == order[0]), options);
        Assert.That(Numbers(root, "rt_minutes")[^1], Is.EqualTo(first.ToRtMinutes(new Irt(50)).Value));
        Assert.That(Texts(root, "group"), Is.EqualTo(interleaved.Select(r => r["run"]).Append(order[0])), "rows in input order");
    }

    // ---- refusals --------------------------------------------------------------------------------

    [Test]
    public void ARunMzLibRefusesFailsTheCallOrIsSkipped()
    {
        var rows = Rows("diann_pbmc_rt.tsv");
        string stdin = string.Join('\n', rows.Select(r => $"{r["run"]}\t{r["retention_time"]}\t{r["indexed_retention_time"]}")
            .Append("--").Append($"{rows[0]["run"]}\t5\t")) + "\n";
        int firstCount = rows.Count(r => r["run"] == rows[0]["run"]);

        Assert.That(Usage(stdin), Is.EqualTo($"Group '{rows[0]["run"]}': Calibration needs at least 20 anchors; got {firstCount}."));

        JsonElement root = Run(stdin, "--on-error", "skip");
        Assert.That(root.GetProperty("failed_count").GetInt32(), Is.EqualTo(3));
        Assert.That(root.GetProperty("fitted_count").GetInt32(), Is.EqualTo(0));
        Assert.That(root.GetProperty("on_error").GetString(), Is.EqualTo("skip"));
        JsonElement group = root.GetProperty("groups")[0];
        Assert.That(group.GetProperty("kind").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(group.GetProperty("residual_sd").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(group.GetProperty("anchor_count").GetInt32(), Is.EqualTo(firstCount), "the anchors sent");
        JsonElement error = group.GetProperty("error");
        Assert.That(error.GetProperty("kind").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("ArgumentException"));
        Assert.That(error.GetProperty("message").GetString(), Is.EqualTo($"Calibration needs at least 20 anchors; got {firstCount}."));
        Assert.That(Numbers(root, "rt_minutes")[0], Is.EqualTo(D(rows[0]["retention_time"])), "the times given still cross");
        Assert.That(Numbers(root, "fitted_irt"), Is.All.Null);
        Assert.That(Flags(root, "is_extrapolated"), Is.All.Null);
        Assert.That(Numbers(root, "irt")[^1], Is.Null, "a to_irt query in a refused group has no iRT");
        Assert.That(Numbers(root, "rt_minutes")[^1], Is.EqualTo(5));
        Assert.That(Caveats(root), Has.Some.StartWith("3 groups were not fitted"));
    }

    [Test]
    public void AnOptionMzLibRefusesFailsTheCallEvenUnderSkip()
    {
        string stdin = string.Join('\n', Rows("diann_pbmc_rt.tsv").Select(r => $"{r["run"]}\t{r["retention_time"]}\t{r["indexed_retention_time"]}")) + "\n";
        Assert.That(Usage(stdin, "--knots", "2", "--on-error", "skip"), Is.EqualTo("Knots must be at least 3."));
        Assert.That(Usage(stdin, "--bandwidth", "1.5"), Is.EqualTo("Bandwidth must be in (0, 1]."));
        Assert.That(Usage(stdin, "--minimum-anchors", "2"), Is.EqualTo("MinimumAnchors must be at least 3."));
        Assert.That(Usage(stdin, "--robustness-iterations", "-1"), Is.EqualTo("RobustnessIterations cannot be negative."));
    }

    [Test]
    public void ABlankAnchorReachesMzLibAndMzLibRefusesIt()
    {
        string stdin = "r\t1\t1\nr\t2\tNA\nr\t3\t3\n";
        Assert.That(Usage(stdin, "--minimum-anchors", "3"), Is.EqualTo("Group 'r': Every anchor's retention time and iRT must be finite."));
        JsonElement root = Run(stdin, "--minimum-anchors", "3", "--on-error", "skip");
        Assert.That(Numbers(root, "irt")[1], Is.Null, "a blank anchor value crosses back as null");
    }

    [Test]
    public void AGroupWithQueriesAndNoAnchorsIsMzLibsRefusal()
    {
        string stdin = "r\t1\t1\nr\t2\t2\nr\t3\t3\n--\nother\t\t2\n";
        Assert.That(Usage(stdin, "--minimum-anchors", "3"), Is.EqualTo("Group 'other': Calibration needs at least 3 anchors; got 0."));
    }

    [TestCase("", "No anchors were given")]
    [TestCase("--\nr\t1\t\n", "No anchors were given")]
    [TestCase("r\t1\n", "has 2 tab-separated fields")]
    [TestCase("r\t1\t2\t3\n", "has 4 tab-separated fields")]
    [TestCase(" \t1\t2\n", "empty group name")]
    [TestCase("r\tearly\t2\n", "Line 1 rt_minutes is not a number: 'early'")]
    [TestCase("r\t1\t2\n--\n--\n", "second '--'")]
    [TestCase("r\t1\t2\n--\nr\t1\t2\n", "with both rt_minutes and irt")]
    [TestCase("r\t1\t2\n--\nr\t\tNA\n", "with neither rt_minutes nor irt")]
    [TestCase("r\t1\t2\n--\nr\tInfinity\t\n", "rt_minutes is not finite")]
    public void MalformedInputIsAUsageError(string stdin, string expected) =>
        Assert.That(Usage(stdin), Does.Contain(expected));

    [TestCase(new[] { "--kind", "cubic" }, "must be 'linear' or 'lowess'; got 'cubic'")]
    [TestCase(new[] { "--kind" }, "--kind was given without a value")]
    [TestCase(new[] { "--knots" }, "--knots was given without a value")]
    [TestCase(new[] { "--bandwidth", "wide" }, "--bandwidth must be a number")]
    [TestCase(new[] { "--minimum-anchors", "1.5" }, "--minimum-anchors must be an integer")]
    [TestCase(new[] { "--on-error", "ignore" }, "--on-error must be 'fail' or 'skip'")]
    public void AMalformedOptionIsAUsageError(string[] options, string expected) =>
        Assert.That(Usage("r\t1\t2\n", options), Does.Contain(expected));

    [Test]
    public void EveryCalibrationKindHasAWireName()
    {
        foreach (IrtCalibrationKind kind in Enum.GetValues<IrtCalibrationKind>())
            Assert.That(Library.KindName(kind), Is.EqualTo(kind.ToString().ToLowerInvariant()));
        Assert.Throws<InvalidOperationException>(() => Library.KindName((IrtCalibrationKind)99));
    }

    [Test]
    public void TheVerbIsDispatched()
    {
        Assert.That(Program.Verbs, Does.Contain("library calibrate-rt"));
    }
}
