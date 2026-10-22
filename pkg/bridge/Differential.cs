using System.Globalization;
using System.Text;

namespace MzLibBridge;

/// <summary>
/// Differential abundance: the experiment's design and the observation table a test reads, both
/// built by mzLib's <c>Quantification.Differential</c> (mzLib #1450, #1455). The bridge parses the
/// wire's tables, calls mzLib, and shapes the result; every rule (coding, estimability, which
/// value enters a sum, the normalization) is mzLib's.
/// </summary>
/// <remarks>
/// <para>
/// <b>Two building blocks, one pivot format.</b> <c>diff design --out</c> writes, per stratum, the
/// design table <c>stats fit --design</c> reads and the contrast weights mzLib's
/// <c>EmpiricalBayes.ModerateContrasts</c> takes; <c>diff observations --responses-out</c> writes
/// the peptide-by-sample table <c>stats fit --responses</c> reads. So the verbs compose through
/// files, with no glue in any binding.
/// </para>
/// <para>
/// <b>What is not here, on purpose.</b> mzLib does not yet turn an observation table into a
/// feature-by-sample matrix at a chosen grain (peptide to protein), choose the rows a test uses, or
/// turn fits into <c>DifferentialResult</c> rows (QuantProject milestone M5b). That step is not
/// written in the bridge (design/PRINCIPLE.md): the responses table is the observation table's own
/// peptide-by-sample shape and nothing more.
/// </para>
/// </remarks>
internal static partial class Differential
{
    // ---- shared input helpers --------------------------------------------------------------------

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

    /// <summary>
    /// A TSV as rows of trimmed cells: a UTF-8 byte-order mark is dropped, blank lines are skipped
    /// and a trailing carriage return is removed, as <c>stats fit</c> reads its tables.
    /// </summary>
    private static List<string[]> ReadTable(string path, string what)
    {
        var rows = File.ReadAllLines(path)
            .Select((line, i) => i == 0 ? line.TrimStart('﻿') : line)
            .Where(line => line.Trim().Length > 0)
            .Select(line => line.TrimEnd('\r').Split('\t').Select(cell => cell.Trim()).ToArray())
            .ToList();
        if (rows.Count == 0)
            throw new Program.UsageException($"The {what} file '{path}' is empty.");
        return rows;
    }

    /// <summary>The first column of a TSV after its header row: the sample ids of a sample sheet or a design file.</summary>
    private static List<string> FirstColumn(string path, string what)
    {
        List<string[]> rows = ReadTable(path, what);
        var ids = rows.Skip(1).Select(r => r[0]).ToList();
        if (ids.Count == 0)
            throw new Program.UsageException($"The {what} file '{path}' has a header but no rows.");
        if (ids.Any(id => id.Length == 0))
            throw new Program.UsageException($"The {what} file '{path}' has a row with an empty first column.");
        return ids;
    }

    /// <summary>
    /// mzLib's refusal of an input it cannot build from, as the caller's mistake (exit 2): a missing
    /// factor value, a run the producer did not write, a duplicate peptide. The parameter suffix .NET
    /// appends is dropped.
    /// </summary>
    private static Program.UsageException Refused(ArgumentException refusal) =>
        new(refusal.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0]);

    /// <summary>Round-trip text for a number in a written table, so a reader gets back exactly mzLib's value.</summary>
    private static string Number(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    /// <summary>
    /// Writes a plain tab-separated table, one line per row, for the bridge's own readers
    /// (<c>stats fit</c>'s). Every cell here is a name or a number that cannot hold a tab or a line
    /// break (names came from TSV cells or stdin fields), so nothing is quoted: <c>stats fit</c>
    /// splits on tabs and would keep a quote as part of a name.
    /// </summary>
    private static void WritePlain(string path, IEnumerable<string> header, IEnumerable<IEnumerable<string>> rows)
    {
        var text = new StringBuilder();
        text.Append(string.Join('\t', header)).Append('\n');
        foreach (IEnumerable<string> row in rows)
            text.Append(string.Join('\t', row)).Append('\n');
        string? directory = Path.GetDirectoryName(Path.GetFullPath(path));
        if (!string.IsNullOrEmpty(directory))
            Directory.CreateDirectory(directory);
        File.WriteAllText(path, text.ToString(), new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
    }

    private static double? Finite(double value) => double.IsFinite(value) ? value : null;

    private static List<List<double>> Rows(double[,] matrix) =>
        Enumerable.Range(0, matrix.GetLength(0))
            .Select(i => Enumerable.Range(0, matrix.GetLength(1)).Select(j => matrix[i, j]).ToList())
            .ToList();
}
