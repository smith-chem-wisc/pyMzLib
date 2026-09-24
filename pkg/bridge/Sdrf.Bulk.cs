using System.Runtime.ExceptionServices;
using Readers;

namespace MzLibBridge;

/// <summary>
/// The one-path-or-many input every per-document SDRF verb shares (<c>validate</c>, <c>assess</c>,
/// <c>samples</c>), implemented once so the three cannot drift apart. The contract is
/// <c>bridge/design/BULK.md</c> §1–§3.
/// </summary>
/// <remarks>
/// <para>
/// <b>Why bulk at all.</b> Each bridge call pays about 120 ms of process start-up and .NET
/// initialisation before mzLib does anything. Validating a 1,236-file corpus one call at a time
/// spends two and a half minutes starting processes; one bulk call spends it once. A binding that
/// wanted the batch faster would otherwise grow its own thread pool, which is the fan-out
/// <c>bridge/design/PARALLELISM.md</c> rules out.
/// </para>
/// <para>
/// <b>Output never depends on <c>--threads</c>.</b> Results land in an array slot per input and are
/// emitted in input order, so the payload is byte-identical at any thread count. The default is
/// <c>1</c> because each worker holds a whole document in memory; a caller who wants throughput
/// says so.
/// </para>
/// <para>
/// This is deliberately private to the SDRF verbs for now. The readers verbs are gaining the same
/// options (<c>G-bulk-readers</c>); when that lands, this should fold into the shared helper rather
/// than become a second one.
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    /// <summary>What <c>--path</c> / <c>--paths-stdin</c>, <c>--threads</c> and <c>--on-error</c> resolved to.</summary>
    /// <param name="IsBulk">True for <c>--paths-stdin</c>: the result is the bulk shape with <c>files[]</c>.</param>
    /// <param name="Paths">Absolute paths, in input order. One entry when <paramref name="IsBulk"/> is false.</param>
    /// <param name="Threads">Documents processed at once; <c>-1</c> means every core.</param>
    /// <param name="SkipFailures">True for <c>--on-error skip</c>.</param>
    internal sealed record BulkInput(bool IsBulk, IReadOnlyList<string> Paths, int Threads, bool SkipFailures);

    /// <summary>One input's outcome: either a value or the failure that replaced it.</summary>
    internal sealed record BulkOutcome<T>(string Path, T? Value, object? Error) where T : class;

    /// <summary>Resolves and validates the bulk options, before any document is read.</summary>
    internal static BulkInput BulkFrom(Program.Arguments arguments)
    {
        bool single = arguments.WasProvided("path");
        bool bulk = arguments.WasProvided("paths-stdin");
        if (single && bulk)
            throw new Program.UsageException(
                "Give --path for one document or --paths-stdin for many, not both.");
        if (!single && !bulk)
            throw new Program.UsageException(
                "Missing required option --path (one document) or --paths-stdin (many, one path " +
                "per line on stdin).");

        int threads = arguments.OptionalInt("threads", 1);
        if (threads == 0 || threads < -1)
            throw new Program.UsageException(
                $"Option --threads must be a positive whole number, or -1 for every core; got {threads}.");

        string onError = arguments.Optional("on-error") ?? "fail";
        if (arguments.WasProvided("on-error") && arguments.Optional("on-error") is null)
            throw new Program.UsageException("Option --on-error needs a value: fail or skip.");
        if (onError is not ("fail" or "skip"))
            throw new Program.UsageException($"Option --on-error must be 'fail' or 'skip'; got '{onError}'.");

        if (!bulk)
        {
            if (onError == "skip")
                throw new Program.UsageException(
                    "--on-error skip applies to --paths-stdin; with one --path a failure is the answer.");
            string path = arguments.Required("path");
            if (!File.Exists(path))
                throw new Program.UsageException($"SDRF file not found: '{path}'.");
            return new BulkInput(false, new[] { Path.GetFullPath(path) }, threads, false);
        }

        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "--paths-stdin was given but stdin held no paths. Supply one SDRF path per line.");

        // A duplicate would be read twice and reported twice under two source_index values, which
        // looks like two experiments that happen to agree perfectly. Refused rather than guessed at.
        var comparer = OperatingSystem.IsWindows() ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
        var seen = new Dictionary<string, int>(comparer);
        var paths = new List<string>(lines.Count);
        bool skip = onError == "skip";
        foreach (string line in lines)
        {
            string trimmed = line.Trim();
            string full = Path.GetFullPath(trimmed);
            if (seen.TryGetValue(full, out int first))
                throw new Program.UsageException(
                    $"'{trimmed}' appears twice on stdin (lines {first + 1} and {paths.Count + 1}). " +
                    "Each document may be given once.");
            if (!skip && !File.Exists(full))
                throw new Program.UsageException($"SDRF file not found: '{trimmed}'.");
            seen[full] = paths.Count;
            paths.Add(full);
        }

        return new BulkInput(true, paths, threads, skip);
    }

    /// <summary>
    /// Runs <paramref name="work"/> over every input with the requested degree of parallelism, and
    /// returns the outcomes in INPUT order regardless of which finished first.
    /// </summary>
    /// <remarks>
    /// Under <c>--on-error fail</c> the failure reported is the one at the LOWEST input position,
    /// not the first to happen, so the same bad batch fails with the same message at any
    /// <c>--threads</c>.
    /// </remarks>
    internal static List<BulkOutcome<T>> RunEach<T>(BulkInput input, Func<string, T> work) where T : class
    {
        var outcomes = new BulkOutcome<T>[input.Paths.Count];
        var failures = new Exception?[input.Paths.Count];

        Parallel.For(0, input.Paths.Count,
            new ParallelOptions { MaxDegreeOfParallelism = input.Threads },
            i =>
            {
                string path = input.Paths[i];
                try
                {
                    if (!File.Exists(path))
                        throw new Program.UsageException($"SDRF file not found: '{path}'.");
                    outcomes[i] = new BulkOutcome<T>(path, work(path), null);
                }
                catch (Exception exception)
                {
                    failures[i] = exception;
                    outcomes[i] = new BulkOutcome<T>(path, null, ErrorFor(exception));
                }
            });

        if (!input.SkipFailures)
        {
            Exception? first = failures.FirstOrDefault(f => f is not null);
            if (first is not null)
                ExceptionDispatchInfo.Capture(first).Throw();
        }

        return outcomes.ToList();
    }

    /// <summary>
    /// A skipped file's <c>error</c>: <c>{kind, message}</c>, with the same three kinds a whole call
    /// can fail with (<c>docs/contributing/conventions.md</c> §2).
    /// </summary>
    private static object ErrorFor(Exception exception)
    {
        Exception cause = Program.Unwrap(exception);
        if (cause is Program.UsageException)
            return new { kind = "usage", message = cause.Message };
        return new { kind = "correctness", message = $"{Program.ClassifyError(cause)}: {cause.Message}" };
    }

    /// <summary>Opens one document and parses every row, as <c>sdrf read</c> does.</summary>
    internal static SdrfDocument Load(string path)
    {
        var document = new SdrfDocument(path);
        document.LoadResults();
        return document;
    }

    /// <summary>
    /// A fixed-schema table crossing as the readers verbs' <c>column_names</c> + <c>columns</c>
    /// (name to values, one value per row), so a binding reuses the columnar type it already has.
    /// </summary>
    internal sealed class ColumnTable
    {
        private readonly List<List<object?>> _values;

        public ColumnTable(params string[] names)
        {
            Names = names.ToList();
            _values = names.Select(_ => new List<object?>()).ToList();
        }

        public List<string> Names { get; }

        public int RowCount => _values.Count == 0 ? 0 : _values[0].Count;

        public void Add(params object?[] row)
        {
            if (row.Length != Names.Count)
                throw new InvalidOperationException(
                    $"A row of {row.Length} values was added to a table of {Names.Count} columns.");
            for (int i = 0; i < row.Length; i++)
                _values[i].Add(row[i]);
        }

        /// <summary>This table with <c>source_index</c> and <c>source_path</c> prepended (BULK.md §2).</summary>
        public static ColumnTable WithSource(params string[] names) =>
            new(new[] { "source_index", "source_path" }.Concat(names).ToArray());

        public Dictionary<string, List<object?>> Columns()
        {
            var columns = new Dictionary<string, List<object?>>(StringComparer.Ordinal);
            for (int i = 0; i < Names.Count; i++)
                columns[Names[i]] = _values[i];
            return columns;
        }
    }
}
