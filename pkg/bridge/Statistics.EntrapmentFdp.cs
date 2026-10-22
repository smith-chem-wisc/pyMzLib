using System.Globalization;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats entrapment-fdp</c>: entrapment estimates of the false discovery proportion among a
/// search's discoveries at each of several q-value thresholds, by mzLib's
/// <see cref="EntrapmentFdp.Sweep"/> (mzLib #1385, 1.0.595). The bridge reads the two candidate
/// tables and makes one mzLib call; every number is mzLib's.
/// </summary>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats entrapment-fdp</c> table, in order.</summary>
    internal static readonly string[] EntrapmentColumns =
    {
        "q_value_threshold", "target_count", "entrapment_count", "lower_bound", "combined", "paired",
    };

    private static readonly string[] TargetFileColumns = { "score", "q_value" };
    private static readonly string[] EntrapmentFileColumns = { "score", "q_value", "partner_score", "partner_q_value", "foreign" };

    /// <summary>
    /// <c>stats entrapment-fdp --targets FILE --entrapments FILE --ratio R [--thresholds Q1,Q2,...]</c>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>targets</b> is a TSV with a header row naming exactly the columns <c>score</c> and
    /// <c>q_value</c>, in any order: one row per original-target candidate
    /// (<see cref="ScoredIdentification"/>). <b>entrapments</b> names <c>score</c>, <c>q_value</c>,
    /// <c>partner_score</c> and <c>partner_q_value</c>, and optionally <c>foreign</c>: one row per
    /// entrapment candidate (<see cref="ScoredEntrapment"/>). Both partner cells blank (or NA) is no
    /// partner; <c>foreign</c> is 0/1 or true/false, and a blank cell or no column is mzLib's default
    /// for <see cref="ScoredEntrapment.IsForeign"/>. Any other column is refused, so a misspelt
    /// optional column cannot silently fall back to that default.
    /// </para>
    /// <para>
    /// <b>ratio</b> has no default in mzLib, so the caller names it. <b>thresholds</b> absent is
    /// <see cref="EntrapmentFdp.DefaultQValueThresholds"/>, read from mzLib and echoed.
    /// </para>
    /// </remarks>
    public static object Entrapment(Program.Arguments arguments)
    {
        string targetsPath = RequireFile(arguments, "targets");
        string entrapmentsPath = RequireFile(arguments, "entrapments");
        RequireValueIfProvided(arguments, "ratio");
        arguments.Required("ratio");
        double ratio = arguments.OptionalDouble("ratio", double.NaN);
        RequireValueIfProvided(arguments, "thresholds");
        bool thresholdsGiven = arguments.WasProvided("thresholds");
        IReadOnlyList<double> thresholds = thresholdsGiven
            ? ParseThresholds(arguments.Optional("thresholds") ?? "")
            : EntrapmentFdp.DefaultQValueThresholds;

        List<ScoredIdentification> targets = ReadTargets(targetsPath);
        List<ScoredEntrapment> entrapments = ReadEntrapments(entrapmentsPath);

        IReadOnlyList<EntrapmentFdpPoint> points;
        try
        {
            points = EntrapmentFdp.Sweep(targets, entrapments, ratio, thresholds);
        }
        catch (ArgumentOutOfRangeException outOfRange)
        {
            string message = outOfRange.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0];
            throw new Program.UsageException(
                $"Option --ratio {ratio.ToString("R", CultureInfo.InvariantCulture)}: {message}");
        }

        var table = new Sdrf.ColumnTable(EntrapmentColumns);
        foreach (EntrapmentFdpPoint point in points)
            table.Add(point.QValueThreshold, point.TargetCount, point.EntrapmentCount,
                Finite(point.LowerBound), Finite(point.Combined), point.Paired is { } paired ? Finite(paired) : null);

        int foreign = entrapments.Count(e => e.IsForeign);
        int partnered = entrapments.Count(e => e.Partner.HasValue);
        return new
        {
            targets_file = Path.GetFullPath(targetsPath),
            entrapments_file = Path.GetFullPath(entrapmentsPath),
            ratio,
            thresholds,
            thresholds_given = thresholdsGiven,
            target_candidate_count = targets.Count,
            entrapment_candidate_count = entrapments.Count,
            foreign_count = foreign,
            partnered_count = partnered,
            threshold_count = thresholds.Count,
            empty_threshold_count = thresholds.Count - points.Count,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = EntrapmentCaveats(points, thresholds.Count, ratio, entrapments.Count, foreign, partnered),
        };
    }

    /// <summary>Comma-separated q-values, each a finite number, in the order given.</summary>
    private static List<double> ParseThresholds(string raw)
    {
        string[] parts = raw.Split(',');
        if (raw.Trim().Length == 0)
            throw new Program.UsageException("Option --thresholds names no q-value. Give them comma-separated, e.g. 0.01,0.05.");
        return parts.Select((p, i) => ParseRequired(p, $"Option --thresholds entry {i + 1}")).ToList();
    }

    private static List<ScoredIdentification> ReadTargets(string path)
    {
        var (rows, index) = ReadCandidateTable(path, "targets", TargetFileColumns, optional: Array.Empty<string>());
        var targets = new List<ScoredIdentification>(rows.Count);
        foreach (var (line, cells) in rows)
            targets.Add(new ScoredIdentification(
                ParseRequired(cells[index["score"]], $"Line {line} of the targets file: score"),
                ParseRequired(cells[index["q_value"]], $"Line {line} of the targets file: q_value")));
        return targets;
    }

    private static List<ScoredEntrapment> ReadEntrapments(string path)
    {
        var (rows, index) = ReadCandidateTable(path, "entrapments", EntrapmentFileColumns, optional: new[] { "foreign" });
        var entrapments = new List<ScoredEntrapment>(rows.Count);
        foreach (var (line, cells) in rows)
        {
            string where = $"Line {line} of the entrapments file";
            double score = ParseRequired(cells[index["score"]], $"{where}: score");
            double q = ParseRequired(cells[index["q_value"]], $"{where}: q_value");
            double partnerScore = ParseOptional(cells[index["partner_score"]], $"{where}: partner_score");
            double partnerQ = ParseOptional(cells[index["partner_q_value"]], $"{where}: partner_q_value");
            if (double.IsNaN(partnerScore) != double.IsNaN(partnerQ))
                throw new Program.UsageException(
                    $"{where} gives a partner {(double.IsNaN(partnerScore) ? "q_value but no score" : "score but no q_value")}. " +
                    "Give both, or leave both blank for an entrapment with no partner.");
            if (double.IsInfinity(partnerScore) || double.IsInfinity(partnerQ))
                throw new Program.UsageException($"{where}: the partner's score and q_value must be finite numbers.");
            ScoredIdentification? partner = double.IsNaN(partnerScore) ? null : new ScoredIdentification(partnerScore, partnerQ);

            string? foreignCell = index.TryGetValue("foreign", out int f) ? cells[f].Trim() : null;
            // A blank cell, or no column, constructs the record without IsForeign: mzLib's default.
            entrapments.Add(string.IsNullOrEmpty(foreignCell)
                ? new ScoredEntrapment(score, q, partner)
                : new ScoredEntrapment(score, q, partner, ParseForeign(foreignCell, where)));
        }
        return entrapments;
    }

    private static bool ParseForeign(string cell, string where) => cell.ToLowerInvariant() switch
    {
        "1" or "true" => true,
        "0" or "false" => false,
        _ => throw new Program.UsageException(
            $"{where}: foreign is '{cell}'; write 1 or true for a foreign entrapment, 0 or false (or blank) for one made from a target."),
    };

    /// <summary>
    /// A TSV whose header names its columns: every one of <paramref name="columns"/> not in
    /// <paramref name="optional"/> must be there, and no other. Returns each data row with its 1-based
    /// line number, and each column's index.
    /// </summary>
    private static (List<(int Line, string[] Cells)> Rows, Dictionary<string, int> Index) ReadCandidateTable(
        string path, string what, string[] columns, string[] optional)
    {
        List<string[]> rows = ReadTsv(path, what);
        List<string> header = rows[0].Select(h => h.Trim()).ToList();
        RequireUnique(header, "column", $"{what} header");
        List<string> unknown = header.Where(h => !columns.Contains(h, StringComparer.Ordinal)).ToList();
        List<string> missing = columns.Where(c => !optional.Contains(c) && !header.Contains(c, StringComparer.Ordinal)).ToList();
        if (unknown.Count > 0 || missing.Count > 0)
            throw new Program.UsageException(
                $"The {what} file '{path}' must name the columns " + string.Join(", ", columns.Where(c => !optional.Contains(c))) +
                (optional.Length > 0 ? ", and optionally " + string.Join(", ", optional) : "") + ", in any order." +
                (missing.Count > 0 ? " Missing: " + string.Join(", ", missing) + "." : "") +
                (unknown.Count > 0 ? " Not a column this file takes: " + string.Join(", ", unknown) + "." : ""));
        var index = header.Select((h, i) => (h, i)).ToDictionary(x => x.h, x => x.i, StringComparer.Ordinal);
        var data = new List<(int, string[])>(rows.Count - 1);
        for (int r = 1; r < rows.Count; r++)
        {
            if (rows[r].Length != header.Count)
                throw new Program.UsageException(
                    $"Line {r + 1} of the {what} file has {rows[r].Length} fields; the header has {header.Count}.");
            data.Add((r + 1, rows[r]));
        }
        return (data, index);
    }

    private static List<string> EntrapmentCaveats(IReadOnlyList<EntrapmentFdpPoint> points, int thresholdCount,
        double ratio, int entrapmentCount, int foreign, int partnered)
    {
        var caveats = new List<string>
        {
            "Wen et al. 2025 (Nat Methods 22:1454), as FDRBench implements it. lower_bound, N_E / (N_T + N_E), can only show that a search fails to control its FDR. combined, N_E (1 + 1/r) / (N_T + N_E), is an upper bound: evidence that it does. paired is a tighter upper bound for 1:1 pairing, each target with exactly one entrapment made from it. None is capped at 1.",
            "A candidate is discovered when its q-value is at or below the threshold. A partner counts as discovered by its own q-value. When an entrapment and its partner are both discovered they are ranked as FDRBench ranks discoveries, lower q-value first and then higher score, and an exact tie counts the entrapment ahead, the conservative choice for an upper bound.",
            "Give one row per unit estimated (precursor, peptide or protein group), decoys removed, each with its best q-value: mzLib does not deduplicate. A partner's score and q-value are copies of its target's row; mzLib does not check them against the targets file.",
            "paired is reported only when the ratio is exactly 1 and some entrapment is not foreign. For 1:1 pairing give 1 itself, not a ratio computed from database sizes after exclusions. The sample estimator N_E / (r N_T), invalid in both directions, is deliberately absent, and the k-fold paired estimator for k entrapments per target is not implemented.",
            "An entrapment generator's exclusions change the inputs in opposite directions (mzLib's EntrapmentFdp.ClassifyExclusion): an entrapment that is really a target peptide comes out of the entrapments file; a target that cannot be traced to one target stays in the targets file but is nobody's partner.",
        };
        int empty = thresholdCount - points.Count;
        if (empty > 0)
            caveats.Add($"{empty} of {thresholdCount} thresholds had no discoveries, where every estimate is 0/0; they are left out of the table.");
        if (points.Count > 0 && points.All(p => p.Paired is null))
            caveats.Add(
                $"paired is null at every threshold. Here the ratio is {ratio.ToString("R", CultureInfo.InvariantCulture)} and {foreign} of {entrapmentCount} entrapments are foreign.");
        if (partnered == 0 && points.Any(p => p.Paired is not null))
            caveats.Add(
                "No entrapment has a partner, so paired is 2 N_E / (N_T + N_E): every discovered entrapment counts as beating an undiscovered partner.");
        return caveats;
    }
}
