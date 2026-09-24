using System.Globalization;
using Readers;

namespace MzLibBridge;

/// <summary>
/// The sample half of an SDRF, per <c>source name</c> (<c>sdrf samples</c>), and ages read into
/// years (<c>sdrf parse-age</c>, and inline in <c>sdrf samples</c>).
/// </summary>
/// <remarks>
/// <para>
/// <b>Why ages are parsed here and not in <c>sdrf read</c>.</b> <c>read</c> promises raw cells,
/// never interpreted, and that promise is what makes it safe for a file nobody has looked at. An
/// age is a sample fact, so it belongs with the samples: <c>sdrf samples</c> carries the parsed
/// age on the same row as the <c>characteristics[age]</c> cell it came from. <c>sdrf parse-age</c>
/// is the same parse over arbitrary cells, for ages that did not come from an SDRF file at all.
/// </para>
/// <para>
/// <b>Two encodings are the bridge's, and both are forced by JSON.</b> A lower bound
/// (<c>&gt;=90Y</c>) has <see cref="SdrfAge.MaxYears"/> = +∞, which JSON cannot carry, so it
/// crosses as <c>null</c> — and <c>precision</c> = <c>LowerBound</c> is what says that null means
/// "unbounded" rather than "unknown". And <see cref="SdrfAge.TryParse"/> answers only yes or no;
/// the <c>refusal</c> reason is a classification of mzLib's "no", taken from mzLib's own reserved
/// words, so that a caller can tell "the curator said not available" from "63, with no unit" without
/// re-reading the cell. The parse itself is mzLib's, untouched.
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    private const string AgeColumn = "characteristics[age]";

    private static readonly string[] AgeFields =
        { "years", "min_years", "max_years", "precision", "follows_specification", "refusal" };

    private static readonly string[] SampleColumns = new[]
        {
            "source_name", "sample_row_count", "column_name", "column_kind", "position", "status", "value",
        }
        .Concat(AgeFields.Select(f => "age_" + f))
        .ToArray();

    // ---- sdrf samples ---------------------------------------------------------------------------

    /// <summary>
    /// <c>sdrf samples (--path FILE | --paths-stdin) [--threads N] [--on-error fail|skip]</c> —
    /// every sample (<c>source name</c>) with the characteristics and factor values its rows agree
    /// on, from mzLib's <see cref="SdrfSampleBlock.BySourceName"/>, as a long table.
    /// </summary>
    /// <remarks>
    /// One row per sample × column × position. A repeated column (<c>characteristics[organism
    /// part]</c> twice) gives two rows with positions 0 and 1 rather than one joined string. A column
    /// the sample's rows DISAGREE about is withheld by mzLib; it still appears, once, with
    /// <c>status</c> = <c>conflicting</c> and a null value, so a caller filtering on
    /// <c>agreed</c> gets values and a caller reading everything cannot mistake a withheld column for
    /// an absent one.
    /// </remarks>
    public static object Samples(Program.Arguments arguments)
    {
        BulkInput input = BulkFrom(arguments);
        List<BulkOutcome<SampledDocument>> outcomes = RunEach(input, SampleOne);

        if (!input.IsBulk)
        {
            SampledDocument only = outcomes[0].Value!;
            var table = new ColumnTable(SampleColumns);
            foreach (object?[] row in only.Rows)
                table.Add(row);

            return new
            {
                path = only.Path,
                sample_count = only.SampleCount,
                row_count = only.RowCount,
                conflict_count = only.ConflictCount,
                problems = only.Problems,
                column_names = table.Names,
                columns = table.Columns(),
                caveats = SamplesCaveats(),
            };
        }

        var bulk = ColumnTable.WithSource(SampleColumns);
        var files = new List<object>(outcomes.Count);
        for (int i = 0; i < outcomes.Count; i++)
        {
            BulkOutcome<SampledDocument> outcome = outcomes[i];
            SampledDocument? s = outcome.Value;
            if (s is not null)
                foreach (object?[] row in s.Rows)
                    bulk.Add(new object?[] { i, outcome.Path }.Concat(row).ToArray());

            files.Add(new
            {
                path = outcome.Path,
                sample_count = s?.SampleCount,
                row_count = s?.RowCount,
                conflict_count = s?.ConflictCount,
                problems = s?.Problems,
                error = outcome.Error,
            });
        }

        var read = outcomes.Where(o => o.Value is not null).Select(o => o.Value!).ToList();
        return new
        {
            file_count = outcomes.Count,
            read_count = read.Count,
            failed_count = outcomes.Count - read.Count,
            sample_count = read.Sum(s => s.SampleCount),
            files,
            column_names = bulk.Names,
            columns = bulk.Columns(),
            caveats = SamplesCaveats(),
        };
    }

    private sealed record SampledDocument(
        string Path, int RowCount, int SampleCount, int ConflictCount, List<string> Problems, List<object?[]> Rows);

    private static SampledDocument SampleOne(string path)
    {
        SdrfDocument document = Load(path);
        IReadOnlyDictionary<string, SdrfSampleBlock> blocks =
            SdrfSampleBlock.BySourceName(document, out IReadOnlyList<string> problems);

        var rows = new List<object?[]>();
        int conflicts = 0;
        foreach (SdrfSampleBlock block in InDocumentOrder(document, blocks))
        {
            // The header's own spelling of "source name", which is how the block keys its cell.
            string? sourceColumn = block.Cells.Keys.FirstOrDefault(k =>
                string.Equals(k, "source name", StringComparison.OrdinalIgnoreCase));

            IEnumerable<(string Column, string Kind)> agreed =
                (sourceColumn is null ? Enumerable.Empty<(string, string)>() : new[] { (sourceColumn, "source_name") })
                .Concat(block.CharacteristicColumns.Select(c => (c, "characteristic")))
                .Concat(block.FactorValueColumns.Select(c => (c, "factor_value")));

            foreach ((string column, string kind) in agreed)
            {
                IReadOnlyList<string> values = block.All(column);
                for (int p = 0; p < values.Count; p++)
                {
                    object?[] age = IsAgeColumn(column) ? AgeCells(values[p]) : new object?[AgeFields.Length];
                    rows.Add(new object?[]
                        { block.SourceName, block.RowCount, column, kind, p, "agreed", values[p] }
                        .Concat(age).ToArray());
                }
            }

            foreach (string column in block.ConflictingColumns)
            {
                conflicts++;
                rows.Add(new object?[]
                    { block.SourceName, block.RowCount, column, KindOf(column), null, "conflicting", null }
                    .Concat(new object?[AgeFields.Length]).ToArray());
            }
        }

        return new SampledDocument(path, document.Results.Count, blocks.Count, conflicts, problems.ToList(), rows);
    }

    /// <summary>
    /// The blocks in the order their source names first appear in the document. mzLib returns a
    /// dictionary, whose enumeration order is not part of its contract; the lookup below uses that
    /// dictionary's own key comparison, so no matching rule is restated here.
    /// </summary>
    private static IEnumerable<SdrfSampleBlock> InDocumentOrder(
        SdrfDocument document, IReadOnlyDictionary<string, SdrfSampleBlock> blocks)
    {
        string? sourceColumn = document.Header.FirstOrDefault(n =>
            string.Equals(n, "source name", StringComparison.OrdinalIgnoreCase));
        if (sourceColumn is null)
            yield break;

        var emitted = new HashSet<SdrfSampleBlock>(ReferenceEqualityComparer.Instance);
        foreach (SdrfRow row in document.Results)
        {
            string? name = row[sourceColumn]?.Trim();
            if (!string.IsNullOrEmpty(name) && blocks.TryGetValue(name, out SdrfSampleBlock? block) && emitted.Add(block))
                yield return block;
        }
    }

    private static string KindOf(string column) =>
        column.StartsWith("factor value[", StringComparison.OrdinalIgnoreCase) ? "factor_value"
        : column.StartsWith("characteristics[", StringComparison.OrdinalIgnoreCase) ? "characteristic"
        : "source_name";

    /// <summary>
    /// <c>characteristics[age]</c>, matched ignoring case like mzLib's sample-column prefixes, so a
    /// mis-cased <c>characteristics[Age]</c> still has its ages read. The validator reports the
    /// casing; refusing to read the age would only lose it.
    /// </summary>
    private static bool IsAgeColumn(string column) =>
        string.Equals(column, AgeColumn, StringComparison.OrdinalIgnoreCase);

    private static List<string> SamplesCaveats() => new()
    {
        "A sample is a source name, matched ignoring case and surrounding whitespace, exactly as mzLib " +
        "keys it: one sample measured in three runs (or fractions) is ONE sample, and sample_row_count " +
        "says how many rows carried it. Across documents a source name is not unique; key on " +
        "(source_index, source_name).",

        "status 'conflicting' means the sample's rows disagree about that column, so mzLib WITHHELD " +
        "it rather than guess which row is right. value is null there. It is not the same as a column " +
        "the document lacks, which has no row at all.",

        "Values are verbatim: reserved words (\"not available\") are real answers and cross as " +
        "themselves. Column names keep the header's own spelling, so characteristics[age] and " +
        "characteristics[Age] are two columns, per the specification.",

        "age_* is filled only on characteristics[age] rows (matched ignoring case). age_years is in " +
        "YEARS (a month is 1/12 year, a week 7/365.25, a day 1/365.25). On an age row, a null " +
        "age_years means mzLib refused the cell and age_refusal says why; on any other row every " +
        "age_* field is null.",
    };

    // ---- sdrf parse-age -------------------------------------------------------------------------

    /// <summary>
    /// <c>sdrf parse-age</c> — each stdin line read as one <c>characteristics[age]</c> cell by
    /// mzLib's <see cref="SdrfAge.TryParse"/>, one output row per input line, in order.
    /// </summary>
    /// <remarks>
    /// Blank lines are KEPT, unlike every other stdin verb: output row <c>i</c> must be input line
    /// <c>i</c>, or a caller zipping the answer back onto its cells would shift every age after the
    /// first blank. A single trailing newline does not make an extra cell.
    /// </remarks>
    public static object ParseAge(Program.Arguments arguments)
    {
        List<string> cells = ReadStdinCells();
        if (cells.Count == 0)
            throw new Program.UsageException(
                "No age cells were provided on stdin. Supply one cell per line, e.g. '58Y'.");

        var table = new ColumnTable(new[] { "cell" }.Concat(AgeFields).ToArray());
        int parsed = 0;
        foreach (string cell in cells)
        {
            object?[] age = AgeCells(cell);
            if (age[0] is not null) parsed++;
            table.Add(new object?[] { cell }.Concat(age).ToArray());
        }

        return new
        {
            cell_count = cells.Count,
            parsed_count = parsed,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = AgeCaveats(),
        };
    }

    /// <summary>stdin split into lines with blank lines preserved; a UTF-8 BOM is dropped.</summary>
    private static List<string> ReadStdinCells()
    {
        string text = Console.In.ReadToEnd().TrimStart('﻿');
        if (text.Length == 0)
            return new List<string>();
        if (text.EndsWith('\n'))
            text = text[..^1];
        return text.Split('\n').Select(l => l.TrimEnd('\r')).ToList();
    }

    /// <summary>
    /// The six age fields for one cell: mzLib's parse, with +∞ as null and a refusal classified.
    /// </summary>
    internal static object?[] AgeCells(string? cell)
    {
        if (SdrfAge.TryParse(cell, out SdrfAge? age))
        {
            return new object?[]
            {
                age.Years,
                age.MinYears,
                double.IsPositiveInfinity(age.MaxYears) ? null : age.MaxYears,
                age.Precision.ToString(),
                age.FollowsSpecification,
                null,
            };
        }

        return new object?[] { null, null, null, null, null, RefusalOf(cell) };
    }

    /// <summary>
    /// Why mzLib said no, from the cell alone. Coarse on purpose: it names the two refusals a curator
    /// can act on (a reserved word was written; a number was written with no unit) and calls
    /// everything else unreadable, rather than second-guessing mzLib's grammar.
    /// </summary>
    private static string RefusalOf(string? cell)
    {
        string text = cell?.Trim() ?? "";
        if (text.Length == 0)
            return "empty";
        if (SdrfValidator.ReservedWords.Any(w => string.Equals(text, w, StringComparison.OrdinalIgnoreCase)))
            return "reserved_word";
        if (double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
            return "no_unit";
        return "unreadable";
    }

    private static List<string> AgeCaveats() => new()
    {
        "Every age is in YEARS: a month is 1/12 year, a week 7/365.25, a day 1/365.25, an hour " +
        "1/(365.25 x 24). years is the age itself, a range's MIDPOINT, or a bound; for a wide cohort " +
        "range (40Y-85Y) the midpoint is a poor stand-in for any one sample - filter on precision.",

        "precision: Exact (58Y, 30Y6M, and a degenerate range such as 40Y-40Y since mzLib #1333), " +
        "Range (both ends), LowerBound (>=90Y, the usual de-identification cap: max_years is null, " +
        "meaning UNBOUNDED, because JSON cannot carry infinity), UpperBound (<1Y: min_years is 0).",

        "follows_specification is false for unambiguous words (\"3 year\", \"6-8 weeks\", " +
        "\"4 hour\"), which are read but are not the specification's nYnMnD grammar.",

        "A bare number (63) is REFUSED, not assumed to be years: 63 years and 63 days are both " +
        "plausible in one study. So are reserved words and free text. refusal says which: empty, " +
        "reserved_word, no_unit or unreadable. A refused cell has no age; never fill it with a guess.",
    };
}
