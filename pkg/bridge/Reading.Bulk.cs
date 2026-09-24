using System.Globalization;
using System.Reflection;
using CsvHelper.Configuration.Attributes;
using Readers;

namespace MzLibBridge;

/// <summary>
/// One path or many: the shared engine every table-reading verb runs through (bridge
/// <c>design/BULK.md</c>).
/// </summary>
/// <remarks>
/// <para>
/// Each reader verb supplies one thing — how to read ONE input into its view (a
/// <see cref="FileTable"/>) — and this file does the rest, identically for all of them: the
/// <c>--path</c> form, the <c>--paths-stdin</c> form, <c>--threads</c>, <c>--on-error</c>, the
/// <c>source_index</c>/<c>source_path</c> long table, the per-file <c>files[]</c> block, and
/// <c>--out</c> streaming. Written once so seven verbs cannot drift on any of it.
/// </para>
/// <para>
/// <b>Why parallelism lives here and not in a binding.</b> A binding that wanted two hundred runs
/// read quickly could fan out two hundred bridge processes, paying .NET start-up and JIT two
/// hundred times while nobody owns the total thread count. Taking the list in one process and
/// reading <c>--threads</c> files at once is the same speed-up paid for once, with the degree of
/// parallelism stated on the wire (bridge <c>design/PARALLELISM.md</c>).
/// </para>
/// <para>
/// <b>The output never depends on <c>--threads</c>.</b> Files are read concurrently but emitted
/// strictly in input order, each file's rows in the file's own order, so <c>--threads 1</c> and
/// <c>--threads 8</c> produce byte-identical envelopes and byte-identical <c>--out</c> tables — a
/// tested property. That is what makes the thread count a resource choice rather than a correctness
/// one, unlike FlashLFQ's, where it changes answers.
/// </para>
/// </remarks>
internal static partial class Reading
{
    /// <summary>The two leading columns of every multi-file table.</summary>
    private static readonly string[] SourceColumns = ["source_index", "source_path"];

    /// <summary>
    /// One input read into one view: its per-file block (bridge <c>design/BULK.md</c> section 3) and
    /// its rows.
    /// </summary>
    /// <remarks>
    /// The rows are produced on demand by <see cref="Row"/> rather than materialised, so a file
    /// written to <c>--out</c> is held once (as mzLib's records) rather than twice. Each cell is
    /// already in its wire shape.
    /// </remarks>
    private sealed class FileTable
    {
        /// <summary>The per-file facts, keyed by wire name, in wire order.</summary>
        public required Dictionary<string, object?> Block { get; init; }

        /// <summary>This view's columns, which are also in <see cref="Block"/> as <c>column_names</c>.</summary>
        public required IReadOnlyList<string> ColumnNames { get; init; }

        /// <summary>The number of rows this read produces (after the window, in single-file mode).</summary>
        public required int RowCount { get; init; }

        /// <summary>
        /// The records those rows came from, when that differs from <see cref="RowCount"/>: a long
        /// view (one row per record per sample, or per match per score) returns more rows than
        /// records. Null means one row per record.
        /// </summary>
        public int? RecordsReturned { get; init; }

        /// <summary>One row's wire-shaped cells, in <see cref="ColumnNames"/> order.</summary>
        public required Func<int, object?[]> Row { get; init; }

        /// <summary>The offset applied (single-file mode only; 0 in bulk).</summary>
        public int Offset { get; init; }

        /// <summary>Whether the window left records behind (single-file mode only).</summary>
        public bool Truncated { get; init; }

        /// <summary>Releases anything the rows still reference, such as a vendor file handle.</summary>
        public Action? Release { get; init; }
    }

    /// <summary>Reads one input into one view.</summary>
    private delegate FileTable ReadOne(string path, Window window);

    /// <summary>What a table verb is, for the shared engine.</summary>
    /// <param name="Name">The verb, as it appears on the wire, for messages.</param>
    /// <param name="Read">How to read one input.</param>
    /// <param name="Columns">The view's columns, known before any file is read, so a multi-file
    /// <c>--out</c> table can write its header first and a batch whose every file failed still says
    /// what its columns would have been.</param>
    /// <param name="ExtraKeys">The per-file keys this verb adds to the common block, so a failed
    /// input's entry carries the same keys as a read one.</param>
    /// <param name="ReportsRows">Whether the verb's rows can outnumber its records (a long view),
    /// so its single-file answer carries <c>row_count</c> beside <c>returned_count</c>. The two
    /// counts are kept apart because <c>offset</c>, <c>limit</c>, <c>record_count</c> and
    /// <c>returned_count</c> all count RECORDS, and a window check that compared records with rows
    /// would call every long read truncated.</param>
    private sealed record TableVerb(
        string Name, ReadOne Read, IReadOnlyList<string> Columns, IReadOnlyList<string> ExtraKeys, bool ReportsRows = false);

    /// <summary>Whether the caller gave a list of paths rather than one.</summary>
    private static bool IsBulk(Program.Arguments arguments) => arguments.WasProvided("paths-stdin");

    /// <summary>
    /// Runs a table verb in whichever form the caller asked for: <c>--path</c> or
    /// <c>--paths-stdin</c>.
    /// </summary>
    /// <param name="arguments">The parsed command line.</param>
    /// <param name="verb">The verb.</param>
    /// <param name="callLevel">Options that apply to the whole call rather than to one file
    /// (<c>ms_order</c>, <c>peaks_included</c>, …), echoed at the top level in both forms.</param>
    private static object RunTable(
        Program.Arguments arguments,
        TableVerb verb,
        IReadOnlyDictionary<string, object?>? callLevel = null)
    {
        callLevel ??= new Dictionary<string, object?>();
        return IsBulk(arguments)
            ? RunBulk(Batch.From(arguments), verb, callLevel)
            : RunSingle(Window.From(arguments), verb, callLevel);
    }

    /// <summary>The <c>--path</c> form: one file, its block flattened into the top level.</summary>
    private static object RunSingle(Window window, TableVerb verb, IReadOnlyDictionary<string, object?> callLevel)
    {
        FileTable table = verb.Read(window.Path, window);
        try
        {
            Func<int, object?[]> row = RowsWithAbsentNulled(table);
            object? written = null;
            if (window.OutputPath is not null)
            {
                using var sink = new TsvSink(window.OutputPath, table.ColumnNames);
                for (int i = 0; i < table.RowCount; i++)
                    sink.Write(row(i));
                written = sink.Complete();
            }

            var result = new Dictionary<string, object?>(table.Block);
            foreach ((string key, object? value) in callLevel)
                result[key] = value;
            result["returned_count"] = written is null ? table.RecordsReturned ?? table.RowCount : 0;
            if (verb.ReportsRows)
                result["row_count"] = written is null ? table.RowCount : 0;
            result["offset"] = table.Offset;
            // True whenever records were left behind, whether by --limit or by --offset. A short
            // answer and a complete one must never look alike.
            result["truncated"] = table.Truncated;
            // Omitted entirely when writing to disk: materialising both would defeat the point.
            result["columns"] = written is null ? Columnar(table.ColumnNames, table.RowCount, row) : null;
            result["output"] = written;
            return result;
        }
        finally
        {
            table.Release?.Invoke();
        }
    }

    /// <summary>
    /// The <c>--paths-stdin</c> form: every input, one long table in input order, and one
    /// <c>files[]</c> entry per input.
    /// </summary>
    private static object RunBulk(
        Batch batch,
        TableVerb verb,
        IReadOnlyDictionary<string, object?> callLevel)
    {
        IReadOnlyList<string> columns = verb.Columns;
        string[] header = [.. SourceColumns, .. columns];

        TsvSink? sink = batch.OutputPath is null ? null : new TsvSink(batch.OutputPath, header);
        List<object?>[]? inline = sink is null ? header.Select(_ => new List<object?>()).ToArray() : null;

        var files = new List<object>(batch.Paths.Count);
        int readCount = 0, failedCount = 0, rowCount = 0, recordsReturned = 0;
        long recordCount = 0;

        try
        {
            foreach ((int index, FileTable? table, Exception? error) in ReadInOrder(batch, verb.Read))
            {
                string sourcePath = Path.GetFullPath(batch.Paths[index]);
                if (table is null)
                {
                    if (!batch.SkipFailures)
                        throw new Program.SourceFailure(index, batch.Paths[index], error!);
                    files.Add(FailedBlock(sourcePath, error!, verb.ExtraKeys));
                    failedCount++;
                    continue;
                }

                try
                {
                    // Every typed view has a fixed column set, and read-records checks up front that
                    // its inputs share one record type, so this cannot fire on a supported input. It
                    // is here so that if it ever does, the table fails loudly instead of shifting one
                    // file's values under another file's column names.
                    if (!table.ColumnNames.SequenceEqual(columns, StringComparer.Ordinal))
                        throw new InvalidOperationException(
                            $"Input {index} ('{batch.Paths[index]}') produced columns that differ from the " +
                            $"batch's ({string.Join(", ", table.ColumnNames)}); a long table needs one column set.");

                    Func<int, object?[]> rowOf = RowsWithAbsentNulled(table);
                    for (int row = 0; row < table.RowCount; row++)
                    {
                        object?[] cells = rowOf(row);
                        if (sink is not null)
                        {
                            sink.Write([index, sourcePath, .. cells]);
                        }
                        else
                        {
                            inline![0].Add(index);
                            inline[1].Add(sourcePath);
                            for (int c = 0; c < cells.Length; c++)
                                inline[c + 2].Add(cells[c]);
                        }
                    }

                    rowCount += table.RowCount;
                    recordsReturned += table.RecordsReturned ?? table.RowCount;
                    recordCount += table.Block["record_count"] is int count ? count : 0;
                    readCount++;
                    files.Add(table.Block);
                }
                finally
                {
                    // Released as soon as it is written, so memory holds at most --threads files
                    // however long the list is.
                    table.Release?.Invoke();
                }
            }
        }
        catch
        {
            // A batch that stopped part-way must not leave a table on disk that looks complete.
            if (sink is not null)
            {
                sink.Dispose();
                File.Delete(sink.Path);
            }

            throw;
        }

        object? written = sink?.Complete();
        sink?.Dispose();

        var result = new Dictionary<string, object?>
        {
            ["file_count"] = batch.Paths.Count,
            ["read_count"] = readCount,
            ["failed_count"] = failedCount,
            ["record_count"] = recordCount,
            ["returned_count"] = written is null ? recordsReturned : 0,
            // Rows can outnumber records in a long view; every bulk answer states both, so a
            // caller never has to know which verbs are long to size the table.
            ["row_count"] = written is null ? rowCount : 0,
            // --threads is deliberately NOT echoed: the envelope must be byte-identical at every
            // thread count, and an echo of the count would be the one byte that is not.
            ["on_error"] = batch.SkipFailures ? "skip" : "fail",
        };
        foreach ((string key, object? value) in callLevel)
            result[key] = value;
        result["column_names"] = header;
        result["columns"] = inline is null
            ? null
            : header.Select((name, i) => (name, values: inline[i])).ToDictionary(p => p.name, p => p.values);
        result["output"] = written;
        result["files"] = files;
        return result;
    }

    /// <summary>
    /// Reads the batch <see cref="Batch.Threads"/> files at a time and yields them in input order.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A sliding window rather than <c>Parallel.ForEach</c>: file <c>i + threads</c> is not started
    /// until file <c>i</c> has been handed on, so however slow one file is, no more than
    /// <see cref="Batch.Threads"/> whole files are ever in memory — the bound <c>--out</c> promises.
    /// A <c>Parallel.ForEach</c> would race ahead and buffer every finished file behind a slow one.
    /// </para>
    /// <para>
    /// A failure is yielded, not thrown, so the caller decides between stopping and recording it.
    /// With one thread each file is read on the calling thread, in turn: the ordinary call pays for
    /// no task machinery at all.
    /// </para>
    /// </remarks>
    private static IEnumerable<(int Index, FileTable? Table, Exception? Error)> ReadInOrder(Batch batch, ReadOne read)
    {
        IReadOnlyList<string> paths = batch.Paths;
        int threads = batch.Threads;
        var pending = new Task<FileTable>?[paths.Count];
        int started = 0;
        int consumed = 0;

        void StartUpTo(int limit)
        {
            for (; started < paths.Count && started < limit; started++)
            {
                string path = paths[started];
                pending[started] = Task.Run(() => read(path, Window.Whole(path)));
            }
        }

        try
        {
            for (int i = 0; i < paths.Count; i++)
            {
                FileTable? table = null;
                Exception? error = null;

                if (threads == 1)
                {
                    try
                    {
                        table = read(paths[i], Window.Whole(paths[i]));
                    }
                    catch (Exception exception)
                    {
                        error = exception;
                    }
                }
                else
                {
                    StartUpTo(i + threads);
                    try
                    {
                        table = pending[i]!.GetAwaiter().GetResult();
                    }
                    catch (Exception exception)
                    {
                        error = exception;
                    }

                    pending[i] = null;
                }

                consumed = i + 1;
                yield return (i, table, error);
            }
        }
        finally
        {
            // Stopped early (--on-error fail): let the reads still in flight finish and release them,
            // so nothing is left holding a file handle after the envelope is written.
            for (int i = consumed; i < started; i++)
            {
                try
                {
                    pending[i]?.GetAwaiter().GetResult().Release?.Invoke();
                }
                catch (Exception)
                {
                    // Deliberately swallowed: the batch already failed on an earlier input, which is
                    // the error being reported. A second, later failure has nothing to add to it.
                }
            }
        }
    }

    /// <summary>
    /// A table's rows with every <c>absent_fields</c> column null, whatever mzLib put there.
    /// </summary>
    /// <remarks>
    /// "Absent" means the file has no column or source for the field (bridge <c>design/BULK.md</c>
    /// section 4), so any value in it is a default mzLib filled in: MsPathFinderT's missing
    /// <c>QValue</c> reads as a perfect 0, a missing <c>Apex_intensity</c> as zero signal. Enforced
    /// here, once, for every verb, so no view can list a field as absent and still hand back its
    /// fabricated values.
    /// </remarks>
    private static Func<int, object?[]> RowsWithAbsentNulled(FileTable table)
    {
        if (table.Block.GetValueOrDefault("absent_fields") is not IReadOnlyList<string> { Count: > 0 } absent)
            return table.Row;

        int[] positions = absent
            .Select(name => table.ColumnNames.ToList().IndexOf(name))
            .Where(index => index >= 0)
            .ToArray();
        if (positions.Length == 0)
            return table.Row;

        return i =>
        {
            object?[] cells = table.Row(i);
            foreach (int position in positions)
                cells[position] = null;
            return cells;
        };
    }

    /// <summary>The <c>files[]</c> entry of an input that could not be read, under <c>--on-error skip</c>.</summary>
    /// <remarks>
    /// Carries every key a read input's entry does, so a binding can use one type for both: the
    /// facts that could not be established are null or empty, and <c>error</c> says why.
    /// </remarks>
    private static Dictionary<string, object?> FailedBlock(string path, Exception error, IReadOnlyList<string> extraKeys)
    {
        var extras = extraKeys.Select(key => (key, (object?)null)).ToArray();
        Dictionary<string, object?> block = Block(
            path, fileType: FileTypeOrNull(path), reader: null, recordCount: null, rowsNotRead: null,
            retentionTimeUnit: null, caveats: [], columnNames: [], absent: [], failed: [], excluded: [],
            extras);
        block["error"] = ErrorOf(error);
        return block;
    }

    /// <summary>
    /// A skipped input's failure, as <c>{kind, type, message}</c>.
    /// </summary>
    /// <remarks>
    /// <c>kind</c> is the three-way split every error on this wire carries — <c>usage</c> (the input
    /// was wrong: not found, not a recognised type, not this view), <c>service_unavailable</c>, or
    /// <c>correctness</c> (mzLib failed to read a file it recognised) — and <c>type</c> is what the
    /// same failure would have crossed the wire as on its own: <c>usage</c>, or the .NET type name.
    /// </remarks>
    private static object ErrorOf(Exception error)
    {
        Exception cause = Program.Unwrap(error);
        string type = Program.ClassifyError(error);
        string kind = cause is Program.UsageException ? "usage"
            : type == Program.ServiceUnavailableType ? "service_unavailable"
            : "correctness";
        return new { kind, type = kind == "usage" ? "usage" : type, message = cause.Message };
    }

    /// <summary>
    /// The per-file block every table verb reports (bridge <c>design/BULK.md</c> section 3), in wire
    /// order, with the verb's own keys after <c>reader</c>.
    /// </summary>
    private static Dictionary<string, object?> Block(
        string path,
        string? fileType,
        string? reader,
        int? recordCount,
        int? rowsNotRead,
        string? retentionTimeUnit,
        IReadOnlyList<string> caveats,
        IReadOnlyList<string> columnNames,
        IReadOnlyList<string> absent,
        IReadOnlyList<string> failed,
        IReadOnlyList<object> excluded,
        params (string Key, object? Value)[] extras)
    {
        var block = new Dictionary<string, object?>
        {
            ["path"] = Path.GetFullPath(path),
            ["file_type"] = fileType,
            ["reader"] = reader,
        };
        foreach ((string key, object? value) in extras)
            block[key] = value;
        block["record_count"] = recordCount;
        block["rows_not_read"] = rowsNotRead;
        block["retention_time_unit"] = retentionTimeUnit;
        block["caveats"] = caveats;
        block["column_names"] = columnNames;
        block["absent_fields"] = absent;
        block["failed_fields"] = failed;
        block["excluded_fields"] = excluded;
        block["error"] = null;
        return block;
    }

    /// <summary>The file type of a path, or null when mzLib cannot name one.</summary>
    private static string? FileTypeOrNull(string path)
    {
        string type = FileTypeOf(path);
        return type == "unknown" ? null : type;
    }

    /// <summary>Builds the columnar payload — one array per column — from rows.</summary>
    private static Dictionary<string, List<object?>> Columnar(
        IReadOnlyList<string> names, int rowCount, Func<int, object?[]> row)
    {
        var values = names.Select(_ => new List<object?>(rowCount)).ToArray();
        for (int r = 0; r < rowCount; r++)
        {
            object?[] cells = row(r);
            for (int c = 0; c < cells.Length; c++)
                values[c].Add(cells[c]);
        }

        var built = new Dictionary<string, List<object?>>(names.Count);
        for (int c = 0; c < names.Count; c++)
            built[names[c]] = values[c];
        return built;
    }

    /// <summary>
    /// A tab-separated table written a row at a time, so a multi-file read streams to disk.
    /// </summary>
    /// <remarks>
    /// The one writer behind every <c>--out</c>, single-file or batch, so quoting, invariant
    /// formatting and the rendering of null and of lists cannot differ between verbs. See
    /// <see cref="Render"/> for why a list is <c>;</c>-joined and null is an empty cell.
    /// </remarks>
    private sealed class TsvSink : IDisposable
    {
        private readonly StreamWriter _writer;
        private readonly CsvHelper.CsvWriter _csv;
        private int _rows;
        private bool _disposed;

        public string Path { get; }

        public TsvSink(string outputPath, IReadOnlyList<string> header)
        {
            Path = System.IO.Path.GetFullPath(outputPath);
            string? directory = System.IO.Path.GetDirectoryName(Path);
            if (!string.IsNullOrEmpty(directory))
                Directory.CreateDirectory(directory);

            var configuration = new CsvHelper.Configuration.CsvConfiguration(CultureInfo.InvariantCulture)
            {
                Delimiter = "\t",
            };
            _writer = new StreamWriter(File.Create(Path));
            _csv = new CsvHelper.CsvWriter(_writer, configuration);

            foreach (string name in header)
                _csv.WriteField(name);
            _csv.NextRecord();
        }

        public void Write(object?[] cells)
        {
            foreach (object? cell in cells)
                _csv.WriteField(Render(cell));
            _csv.NextRecord();
            _rows++;
        }

        /// <summary>Flushes the table and reports where it went.</summary>
        public object Complete()
        {
            _csv.Flush();
            _writer.Flush();
            return new { path = Path, format = "tsv", row_count = _rows };
        }

        public void Dispose()
        {
            if (_disposed)
                return;
            _disposed = true;
            _csv.Dispose();
            _writer.Dispose();
        }
    }

    // ---------------------------------------------------------------------------------------
    // The --paths-stdin request
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The <c>--paths-stdin</c>/<c>--threads</c>/<c>--on-error</c>/<c>--out</c> set, parsed and
    /// validated once for every table verb.
    /// </summary>
    private sealed record Batch(IReadOnlyList<string> Paths, int Threads, bool SkipFailures, string? OutputPath)
    {
        public static Batch From(Program.Arguments arguments)
        {
            // --paths-stdin is a flag. A value after it would be a path the caller meant to be read,
            // silently ignored in favour of stdin.
            if (arguments.Optional("paths-stdin") is { } stray)
                throw new Program.UsageException(
                    $"--paths-stdin takes no value (got '{stray}'); write the paths to stdin, one per line.");
            if (arguments.WasProvided("path"))
                throw new Program.UsageException(
                    "--path and --paths-stdin are mutually exclusive: give one file with --path, or many on stdin.");
            if (arguments.WasProvided("offset") || arguments.WasProvided("limit"))
                throw new Program.UsageException(
                    "--offset and --limit window ONE file, so they cannot be combined with --paths-stdin. " +
                    "Read each file with --path to window it, or write the batch to --out.");

            (int threads, bool skip) = ConcurrencyOptions(arguments);

            RequireValueIfProvided(arguments, "out");
            string? outputPath = arguments.Optional("out");

            // Trimmed, because a path pasted with a trailing space would otherwise be "not found"
            // for a reason nobody can see; blank lines are dropped by ReadStdinLines itself.
            List<string> paths = Program.ReadStdinLines().Select(line => line.Trim()).ToList();
            if (paths.Count == 0)
                throw new Program.UsageException(
                    "--paths-stdin was given but no paths arrived on stdin; write one path per line.");

            // Duplicates are refused rather than read twice: a repeated path is almost always a
            // mistake in how the list was built, and reading it twice would double its records in
            // every total without a word.
            StringComparer comparer = OperatingSystem.IsWindows() ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
            var seen = new Dictionary<string, int>(comparer);
            for (int i = 0; i < paths.Count; i++)
            {
                string full = System.IO.Path.GetFullPath(paths[i]);
                if (seen.TryGetValue(full, out int first))
                    throw new Program.UsageException(
                        $"Input {i} ('{paths[i]}') repeats input {first}; each file may be listed once.");
                seen[full] = i;
            }

            // The same rule --path has: --out may not overwrite anything being read.
            if (outputPath is not null && seen.TryGetValue(System.IO.Path.GetFullPath(outputPath), out int clash))
                throw new Program.UsageException(
                    $"Option --out must differ from every input: writing to '{outputPath}' would overwrite input {clash}.");

            return new Batch(paths, threads, skip, string.IsNullOrWhiteSpace(outputPath) ? null : outputPath);
        }
    }

    /// <summary>
    /// <c>--threads</c> and <c>--on-error</c>, validated the same way in both forms.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <c>--threads</c> defaults to 1. Every reader materialises the whole file, so N threads means up
    /// to N whole files in memory at once, and a caller who wants throughput can say so; since the
    /// output never depends on it, the default costs speed, never correctness. <c>-1</c> means one
    /// per core. With <c>--path</c> it is accepted and has nothing to do.
    /// </para>
    /// <para>
    /// <c>--on-error skip</c> is refused with <c>--path</c>: skipping the only file would turn a
    /// failure into an empty success, which is the one outcome this bridge exists to prevent.
    /// </para>
    /// </remarks>
    private static (int Threads, bool Skip) ConcurrencyOptions(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "threads");
        int threads = arguments.OptionalInt("threads", 1);
        if (threads == 0 || threads < -1)
            throw new Program.UsageException(
                $"Option --threads must be 1 or more, or -1 for one per core; got {threads}.");
        if (threads == -1)
            threads = Environment.ProcessorCount;

        RequireValueIfProvided(arguments, "on-error");
        string onError = arguments.Optional("on-error") ?? "fail";
        bool skip = onError switch
        {
            "fail" => false,
            "skip" => true,
            _ => throw new Program.UsageException(
                $"Option --on-error must be 'fail' or 'skip'; got '{onError}'."),
        };

        return (threads, skip);
    }

    // ---------------------------------------------------------------------------------------
    // absent_fields: what a file's header does not have
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The properties of a CsvHelper-mapped record type whose column this file's header does not
    /// have — so every value of that property is mzLib's default, not something the file said.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is the general form of the <c>absent_fields</c> list (bridge <c>design/BULK.md</c>
    /// section 4), and it reads mzLib's own declaration rather than a table kept here: a property
    /// marked <c>[Optional]</c> with a <c>[Name(...)]</c> is one mzLib reads when the column exists
    /// and leaves at its default when it does not. Checking the header for those names is exactly
    /// the question CsvHelper answered while reading. The cases this catches at mzLib 1.0.592 include
    /// <c>mbr_score</c> on a current FlashLFQ peaks table, which has no <c>MBR Score</c> column
    /// (#1345); the seven columns #1345 added, on a table written before them; and MsPathFinderT's
    /// <c>q_value</c> on a targets-only file, where mzLib's non-nullable double reads 0 for a column
    /// that is not there.
    /// </para>
    /// <para>
    /// Deliberately conservative. It says nothing about a type with no <c>[Name]</c> attributes, a
    /// directory or compressed input, or a file whose first line lacks any of the type's REQUIRED
    /// columns — a preamble, say, or a header mzLib matched in a way this simple reading does not.
    /// Names compare case-insensitively after trimming, so an ambiguity can only make a column look
    /// present. An empty answer means "no basis to say", never "nothing is absent".
    /// </para>
    /// </remarks>
    private static HashSet<string> AbsentProperties(Type? recordType, string path)
    {
        var absent = new HashSet<string>(StringComparer.Ordinal);
        if (recordType is null || !File.Exists(path) || path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
            return absent;

        var mapped = new List<(PropertyInfo Property, string[] Names, bool Optional)>();
        bool anyNamed = false;
        foreach (PropertyInfo property in recordType.GetProperties(BindingFlags.Public | BindingFlags.Instance))
        {
            if (!property.CanWrite || property.GetIndexParameters().Length > 0 || property.IsDefined(typeof(IgnoreAttribute), true))
                continue;
            // Only a property with an explicit [Name] is judged: an unnamed [Optional] one may be
            // filled by the reader itself rather than from a column (MsPathFinderT's
            // FileNameWithoutExtension is set from the file path), and would read as absent here.
            NameAttribute? name = property.GetCustomAttribute<NameAttribute>(true);
            if (name is null)
                continue;
            anyNamed = true;
            mapped.Add((property, name.Names, property.IsDefined(typeof(OptionalAttribute), true)));
        }

        if (!anyNamed)
            return absent;

        string? headerLine;
        try
        {
            headerLine = File.ReadLines(path).FirstOrDefault(line => !string.IsNullOrWhiteSpace(line));
        }
        catch (IOException)
        {
            return absent;
        }

        if (headerLine is null)
            return absent;

        char delimiter = headerLine.Contains('\t') ? '\t' : ',';
        var header = new HashSet<string>(
            headerLine.Split(delimiter).Select(cell => cell.Trim().Trim('"')), StringComparer.OrdinalIgnoreCase);

        // Not the header this type reads (a preamble, or a mapping this reading does not model):
        // no basis to call anything absent.
        if (mapped.Any(m => !m.Optional && !m.Names.Any(header.Contains)))
            return absent;

        foreach ((PropertyInfo property, string[] names, bool optional) in mapped)
        {
            if (optional && !names.Any(header.Contains))
                absent.Add(property.Name);
        }

        return absent;
    }

    /// <summary>The header cells of a tab-separated file, or an empty set when it cannot be read.</summary>
    private static IReadOnlyList<string> HeaderCells(string path)
    {
        try
        {
            string? line = File.ReadLines(path).FirstOrDefault(l => !string.IsNullOrWhiteSpace(l));
            return line is null ? [] : line.Split('\t').Select(cell => cell.Trim()).ToList();
        }
        catch (IOException)
        {
            return [];
        }
    }
}
