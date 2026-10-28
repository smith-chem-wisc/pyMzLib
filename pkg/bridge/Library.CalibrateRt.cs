using System.Globalization;
using Chromatography.RetentionTimeCalibration;
using MzLibUtil;

namespace MzLibBridge;

/// <summary>
/// <c>library calibrate-rt</c>: one run's retention times in minutes mapped onto a spectral library's
/// iRT scale, and back, by mzLib's <see cref="IrtCalibration"/> (mzLib #1387, 1.0.596).
/// </summary>
/// <remarks>
/// <para>
/// <b>Anchors and queries in one call.</b> mzLib's fitted <see cref="IrtCalibrationModel"/> keeps its
/// knots in private fields and has an internal constructor, so the curve cannot be returned and rebuilt
/// by a later call. The caller therefore sends the anchors that fit each group (typically one run)
/// together with the times it wants converted, and the bridge makes one <see cref="IrtCalibration.Fit"/>
/// per group and converts that group's queries with the model it returned. Fitting is deterministic
/// (mzLib sorts the anchors first), so splitting a job across calls gives the same map.
/// </para>
/// <para>
/// <b>stdin</b> is two lists separated by a line holding only <c>--</c> (bridge <c>design/BULK.md</c>
/// section 2): anchors first, then queries. Both use the same three tab-separated fields,
/// <c>group</c>, <c>rt_minutes</c>, <c>irt</c>. An anchor gives both numbers; a query gives exactly one
/// and leaves the other blank, and the blank one is what mzLib computes. Keeping the two lists apart
/// means a missing value can never turn an anchor into a query: an anchor with a blank or NA field
/// reaches mzLib as NaN, and mzLib refuses it.
/// </para>
/// <para>
/// <b>Options are mzLib's.</b> Each of <c>--kind</c>, <c>--minimum-anchors</c>, <c>--bandwidth</c>,
/// <c>--knots</c> and <c>--robustness-iterations</c> sets one member of
/// <see cref="IrtCalibrationOptions"/> and is applied with <c>with</c> to mzLib's own defaults only when
/// given. mzLib validates them (<see cref="ArgumentOutOfRangeException"/>), and that refusal fails the
/// call whatever <c>--on-error</c> says, because it is the same for every group.
/// </para>
/// </remarks>
internal static class Library
{
    /// <summary>Row columns of the <c>library calibrate-rt</c> table, in order.</summary>
    internal static readonly string[] CalibrateColumns =
        { "group", "role", "rt_minutes", "irt", "fitted_irt", "is_extrapolated" };

    private sealed class Group(string name)
    {
        public string Name { get; } = name;
        public List<(RtMinutes Rt, Irt Irt)> Anchors { get; } = new();
    }

    private readonly record struct Line(Group Group, string Role, double Rt, double Irt);

    /// <summary>
    /// <c>library calibrate-rt [--kind linear|lowess] [--minimum-anchors N] [--bandwidth F] [--knots N]
    /// [--robustness-iterations N] [--on-error fail|skip]</c>, with the anchors and queries on stdin.
    /// </summary>
    public static object CalibrateRt(Program.Arguments arguments)
    {
        IrtCalibrationOptions options = OptionsFrom(arguments);
        Proteins.OnError onError = Proteins.OnErrorFrom(arguments, allowSkip: true);

        List<string> stdin = Program.ReadStdinLines();
        var groups = new Dictionary<string, Group>(StringComparer.Ordinal);
        var order = new List<Group>();
        var lines = new List<Line>();
        bool queries = false;
        for (int i = 0; i < stdin.Count; i++)
        {
            if (stdin[i].Trim() == "--")
            {
                if (queries)
                    throw new Program.UsageException(
                        $"Line {i + 1} of stdin is a second '--'. Give the anchors, one '--' line, then the queries.");
                queries = true;
                continue;
            }
            string[] fields = stdin[i].Split('\t');
            if (fields.Length != 3)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin has {fields.Length} tab-separated field{(fields.Length == 1 ? "" : "s")}; " +
                    "expected 3: group, rt_minutes, irt.");
            string name = fields[0].Trim();
            if (name.Length == 0)
                throw new Program.UsageException($"Line {i + 1} of stdin has an empty group name.");
            double rt = Number(fields[1], $"Line {i + 1} rt_minutes", out bool rtGiven);
            double irt = Number(fields[2], $"Line {i + 1} irt", out bool irtGiven);

            if (!groups.TryGetValue(name, out Group? group))
            {
                group = new Group(name);
                groups[name] = group;
                order.Add(group);
            }

            if (!queries)
            {
                // Blank or NA crosses as NaN: mzLib refuses a non-finite anchor, and its message says so.
                group.Anchors.Add((new RtMinutes(rt), new Irt(irt)));
                lines.Add(new Line(group, "anchor", rt, irt));
                continue;
            }
            if (rtGiven == irtGiven)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin is a query with {(rtGiven ? "both" : "neither")} rt_minutes {(rtGiven ? "and" : "nor")} irt. " +
                    "Give one and leave the other blank: the blank one is what the calibration computes.");
            if (!double.IsFinite(rtGiven ? rt : irt))
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin is a query whose {(rtGiven ? "rt_minutes" : "irt")} is not finite.");
            lines.Add(new Line(group, rtGiven ? "to_irt" : "to_minutes", rt, irt));
        }
        if (!lines.Any(l => l.Role == "anchor"))
            throw new Program.UsageException(
                "No anchors were given on stdin. Give one per line, group<TAB>rt_minutes<TAB>irt, then optionally " +
                "a '--' line and the queries.");

        var models = new Dictionary<Group, IrtCalibrationModel>();
        var errors = new Dictionary<Group, Exception>();
        foreach (Group group in order)
        {
            try
            {
                models[group] = IrtCalibration.Fit(group.Anchors, options);
            }
            catch (ArgumentOutOfRangeException invalid)
            {
                // An option out of range: the same for every group, so never a per-group skip.
                throw new Program.UsageException(MzLibMessage(invalid));
            }
            catch (ArgumentException refused) when (refused is not ArgumentNullException)
            {
                // Too few anchors, a non-finite one, no time span, or a fit that does not increase:
                // mzLib's refusal of the caller's anchors for this group.
                if (onError == Proteins.OnError.Fail)
                    throw new Program.UsageException($"Group '{group.Name}': {MzLibMessage(refused)}");
                errors[group] = refused;
            }
        }

        var table = new Sdrf.ColumnTable(CalibrateColumns);
        int extrapolated = 0;
        foreach (Line line in lines)
        {
            if (!models.TryGetValue(line.Group, out IrtCalibrationModel? model))
            {
                table.Add(line.Group.Name, line.Role, Given(line.Rt), Given(line.Irt), null, null);
                continue;
            }
            if (line.Role == "anchor")
            {
                var rt = new RtMinutes(line.Rt);
                table.Add(line.Group.Name, line.Role, line.Rt, line.Irt, model.ToIrt(rt).Value, model.IsExtrapolated(rt));
                continue;
            }
            // A query: mzLib converts the time given, and judges extrapolation on the run-minutes side.
            RtMinutes minutes = line.Role == "to_irt" ? new RtMinutes(line.Rt) : model.ToRtMinutes(new Irt(line.Irt));
            double irt = line.Role == "to_irt" ? model.ToIrt(minutes).Value : line.Irt;
            bool outside = model.IsExtrapolated(minutes);
            if (outside)
                extrapolated++;
            table.Add(line.Group.Name, line.Role, minutes.Value, irt, null, outside);
        }

        var groupBlocks = order.Select(group =>
        {
            models.TryGetValue(group, out IrtCalibrationModel? model);
            errors.TryGetValue(group, out Exception? error);
            return new
            {
                group = group.Name,
                kind = model is null ? null : KindName(model.Kind),
                anchor_count = model?.AnchorCount ?? group.Anchors.Count,
                query_count = lines.Count(l => l.Group == group && l.Role != "anchor"),
                residual_sd = model is null ? null : Finite(model.ResidualSd),
                first_anchor_minutes = model?.FirstAnchor.Value,
                last_anchor_minutes = model?.LastAnchor.Value,
                error = error is null
                    ? null
                    : new { kind = "usage", type = error.GetType().Name, message = MzLibMessage(error) },
            };
        }).ToList();

        return new
        {
            kind = KindName(options.Kind),
            minimum_anchors = options.MinimumAnchors,
            bandwidth = options.Bandwidth,
            knots = options.Knots,
            robustness_iterations = options.RobustnessIterations,
            on_error = onError == Proteins.OnError.Fail ? "fail" : "skip",
            group_count = order.Count,
            fitted_count = models.Count,
            failed_count = errors.Count,
            anchor_count = lines.Count(l => l.Role == "anchor"),
            query_count = lines.Count(l => l.Role != "anchor"),
            row_count = table.RowCount,
            groups = groupBlocks,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = CalibrateCaveats(errors.Count, extrapolated),
        };
    }

    /// <summary>
    /// mzLib's <see cref="IrtCalibrationOptions"/>, each member set only when its option was given, so
    /// every other member keeps mzLib's own default.
    /// </summary>
    private static IrtCalibrationOptions OptionsFrom(Program.Arguments arguments)
    {
        var options = new IrtCalibrationOptions();
        if (arguments.WasProvided("kind"))
        {
            string? kind = arguments.Optional("kind");
            options = options with
            {
                Kind = kind switch
                {
                    "linear" => IrtCalibrationKind.Linear,
                    "lowess" => IrtCalibrationKind.Lowess,
                    null or "" => throw new Program.UsageException("Option --kind was given without a value."),
                    _ => throw new Program.UsageException($"Option --kind must be 'linear' or 'lowess'; got '{kind}'."),
                },
            };
        }
        if (arguments.WasProvided("minimum-anchors"))
            options = options with { MinimumAnchors = arguments.OptionalInt("minimum-anchors", 0) };
        if (arguments.WasProvided("bandwidth"))
            options = options with { Bandwidth = arguments.OptionalDouble("bandwidth", 0) };
        if (arguments.WasProvided("knots"))
            options = options with { Knots = arguments.OptionalInt("knots", 0) };
        if (arguments.WasProvided("robustness-iterations"))
            options = options with { RobustnessIterations = arguments.OptionalInt("robustness-iterations", 0) };
        return options;
    }

    /// <summary>
    /// mzLib's <see cref="IrtCalibrationKind"/> as the wire writes it. Every value is named here, so a
    /// value mzLib adds fails a test on the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string KindName(IrtCalibrationKind kind) => kind switch
    {
        IrtCalibrationKind.Linear => "linear",
        IrtCalibrationKind.Lowess => "lowess",
        _ => throw new InvalidOperationException(
            $"IrtCalibrationKind.{kind} has no wire name. Name it in Library.KindName."),
    };

    /// <summary>A numeric field: blank, NA and NaN are not given (NaN); anything else unparseable is refused.</summary>
    private static double Number(string raw, string what, out bool given)
    {
        string cell = raw.Trim();
        given = !(cell.Length == 0
            || cell.Equals("NA", StringComparison.OrdinalIgnoreCase)
            || cell.Equals("NaN", StringComparison.OrdinalIgnoreCase));
        if (!given)
            return double.NaN;
        if (!double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
            throw new Program.UsageException(
                $"{what} is not a number: '{cell}'. Leave it blank, or write NA, for a missing value.");
        return value;
    }

    private static double? Given(double value) => double.IsNaN(value) ? null : value;

    private static double? Finite(double value) => double.IsFinite(value) ? value : null;

    /// <summary>mzLib's message without .NET's "(Parameter 'x')" suffix.</summary>
    private static string MzLibMessage(Exception exception) =>
        exception.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0];

    private static List<string> CalibrateCaveats(int failed, int extrapolated)
    {
        var caveats = new List<string>
        {
            "iRT is the scale of the library the anchors' iRT values came from. A calibration maps one run onto that library's scale, so iRT values from two libraries are not comparable, and the same run calibrated against two libraries gives two maps.",
            "The anchors are the caller's choice: mzLib asks for the observed apex retention time, in minutes, of confidently identified precursors paired with their library iRT (IrtCalibration.cs:29-31). The bridge applies no q-value cutoff and removes no duplicate. The fit is robust to some wrong anchors (binned medians, then Tukey bisquare reweighting), not to a wrong choice of them.",
            "residual_sd is mzLib's robust residual SD, in iRT units, of the anchors that fit: residuals beyond 3 robust SDs are set aside before it is computed (IrtCalibration.InlierSd). It understates the spread when many anchors are wrong, and a LOWESS fit to few anchors can pass through nearly all of them, so a residual_sd near 0 there means an overfit, not a perfect run. irt minus fitted_irt on the anchor rows shows every residual.",
            "Outside each group's first and last anchor the map continues as a straight line, not a fit (IrtCalibrationModel.IsExtrapolated). Read is_extrapolated before trusting a converted time there.",
            "The fitted curve itself does not cross the wire: mzLib keeps its knots private and its constructor internal (IrtCalibrationModel), so a calibration cannot be returned and reused. Send the queries in the same call as their anchors; refitting the same anchors gives the same map, because mzLib sorts them first.",
        };
        if (extrapolated > 0)
            caveats.Add($"{extrapolated} quer{(extrapolated == 1 ? "y falls" : "ies fall")} outside the anchors' span, where the map is extrapolated (is_extrapolated).");
        if (failed > 0)
            caveats.Add($"{failed} group{(failed == 1 ? " was" : "s were")} not fitted (on-error skip): mzLib refused the anchors, as each group's error says. Their rows carry the times given and null for the rest.");
        return caveats;
    }
}
