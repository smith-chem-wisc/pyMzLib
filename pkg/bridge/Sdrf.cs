using Readers;

namespace MzLibBridge;

/// <summary>
/// SDRF-Proteomics experimental-design files: read one, or pool several into one table.
/// </summary>
/// <remarks>
/// <para>
/// SDRF is the format that says what a run <i>is</i> — which sample, which organism part, which
/// replicate, which instrument settings. Every other reader here answers "what did the search
/// find"; this one answers "what was searched", which is the half a caller needs to group results
/// across experiments.
/// </para>
/// <para>
/// <b>The payload is row-major, and that is not a stylistic choice.</b> Every other readers verb
/// returns <c>columns</c>, a name-to-values map. SDRF cannot use one: its column names are DATA
/// rather than a schema, and they REPEAT — 649 files in the curated corpus carry
/// <c>comment[modification parameters]</c> more than once, up to eight times in one file, and one
/// file repeats an EMPTY name 23 times. A map keyed by name would silently drop every occurrence
/// but one. So <c>column_names</c> is a list that may contain duplicates, <c>rows</c> is a list of
/// cell lists, and position is what links them.
/// </para>
/// <para>
/// <b>Cells cross as raw strings and are never interpreted.</b> mzLib's reader makes the same
/// promise and for the same reason: the SDRF key=value grammar (<c>NT=…;AC=…</c>) is not safely
/// detectable from a cell's shape, because <c>comment[file uri]</c> cells in the corpus carry
/// pre-signed download URLs whose query strings contain <c>Signature=</c>, <c>Expires=</c> and
/// <c>Id=</c> — 5,750 occurrences each. Anything that decoded cells "that look like CV terms"
/// would decode those. Interpretation is a projection to be layered on top, not baked in here.
/// </para>
/// <para>
/// <b>What this does not do yet is validate.</b> mzLib models the specification's structural
/// rules in <c>SdrfValidator</c> and its drift rules in <c>SdrfDriftLint</c>. Both are public as
/// of mzLib #1207, which the pinned mzLib includes, so they can be reached from here. They are not
/// yet exposed as verbs. When they are, the verbs must call those two entry points rather than
/// re-derive their rules: a second implementation — in C#, then again in Rust and R — is exactly
/// the per-binding repair this bridge exists to avoid. See <c>bridge/UPSTREAM.md</c> (U8).
/// </para>
/// </remarks>
internal static class Sdrf
{
    /// <summary>
    /// <c>sdrf read --path FILE [--limit N] [--offset N]</c> — one SDRF document as a header and
    /// its rows.
    /// </summary>
    public static object Read(Program.Arguments arguments)
    {
        string path = arguments.Required("path");
        (int offset, int limit) = WindowFrom(arguments);

        var document = new SdrfDocument(path);
        document.LoadResults();

        IReadOnlyList<SdrfRow> all = document.Results;
        IReadOnlyList<SdrfRow> selected = Slice(all, offset, limit, out bool truncated);

        return new
        {
            path,
            column_names = document.Header.ToList(),
            row_count = all.Count,
            returned_count = selected.Count,
            offset,
            truncated,
            rows = selected.Select(r => r.Cells.ToList()).ToList(),
            caveats = ReadCaveats(document, all),
        };
    }

    /// <summary>
    /// <c>sdrf pool [--out FILE] [--limit N] [--offset N]</c> — several SDRF documents pooled into
    /// one table, their paths supplied on stdin, one per line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Documents arrive on stdin, one per line, tab-separated: <c>path[\tlabel]</c>. stdin rather
    /// than repeated <c>--path</c> options because <see cref="Program.Arguments"/> holds named
    /// options in a dictionary, so a repeated option keeps only the last value — <c>--path a
    /// --path b</c> would pool <c>b</c> alone and report success. It is also what the two
    /// <c>quant</c> verbs already use for their many-file input, so this is the established shape
    /// rather than a new one.
    /// </para>
    /// <para>
    /// The optional label is the provenance value stamped into
    /// <c>comment[source document]</c>. It is worth supplying: mzLib's default is
    /// <c>containing-folder/file-stem</c> — deliberately, since a folder is what usually
    /// distinguishes one search's output from another's — but that makes the value depend on where
    /// the files happen to sit, so the same two documents pooled from a different directory
    /// produce a different table. Either every line carries a label or none does; a partial set is
    /// refused rather than silently mixing the two schemes.
    /// </para>
    /// </remarks>
    public static object Pool(Program.Arguments arguments)
    {
        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No SDRF files were provided on stdin. Supply one path per line, optionally " +
                "followed by a tab and a provenance label, e.g. 'PXD000070.sdrf.tsv'.");

        var paths = new List<string>(lines.Count);
        var labels = new List<string>(lines.Count);
        foreach (string line in lines)
        {
            string[] fields = line.Split('\t');
            string path = fields[0].Trim();
            if (path.Length == 0)
                throw new Program.UsageException("A stdin line has an empty file path.");
            if (!File.Exists(path))
                throw new Program.UsageException($"SDRF file not found: '{path}'.");

            paths.Add(path);
            if (fields.Length > 1 && fields[1].Trim().Length > 0)
                labels.Add(fields[1].Trim());
        }

        if (labels.Count != 0 && labels.Count != paths.Count)
            throw new Program.UsageException(
                $"{labels.Count} of {paths.Count} stdin lines carry a label. Label every document " +
                "or none: a partial set would stamp some rows with your name for the document and " +
                "others with a path-derived default, and the two are not comparable.");

        (int offset, int limit) = WindowFrom(arguments);
        string? outputPath = arguments.Optional("out");
        if (string.IsNullOrWhiteSpace(outputPath))
            outputPath = null;

        var collection = new SdrfCollection(
            paths.Select(p => new SdrfDocument(p)),
            labels.Count == 0 ? null : labels);
        SdrfDocument merged = collection.Merge();

        IReadOnlyList<SdrfRow> all = merged.Results;
        IReadOnlyList<SdrfRow> selected = Slice(all, offset, limit, out bool truncated);

        object? written = null;
        if (outputPath is not null)
        {
            // The WHOLE merged document, never the windowed slice. --limit and --offset are a
            // transport convenience for the payload; writing only the visible rows would produce a
            // file that silently disagrees with the one the caller asked to pool.
            merged.WriteResults(outputPath);
            written = new { path = outputPath, row_count = all.Count };
        }

        return new
        {
            document_count = collection.Count,
            paths,
            labels = collection.Labels.ToList(),
            column_names = merged.Header.ToList(),
            row_count = all.Count,
            returned_count = selected.Count,
            offset,
            truncated,
            rows = selected.Select(r => r.Cells.ToList()).ToList(),
            written,
            caveats = PoolCaveats(collection, labelsWereSupplied: labels.Count != 0),
        };
    }

    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The truths about an SDRF payload that a caller cannot see from the data and would otherwise
    /// get wrong.
    /// </summary>
    private static List<string> ReadCaveats(SdrfDocument document, IReadOnlyList<SdrfRow> rows)
    {
        var caveats = new List<string>
        {
            "Cells are RAW STRINGS, never interpreted. The SDRF key=value grammar " +
            "(\"NT=Oxidation;AC=UNIMOD:35\") is left exactly as written, because it cannot be told " +
            "apart from a cell that merely contains '=' and ';' - comment[file uri] routinely " +
            "carries pre-signed URLs that do.",

            "column_names is a LIST, not a set: names may repeat, and their order is part of the " +
            "document. comment[modification parameters] legitimately appears many times in one " +
            "file, so a name identifies a POSITION rather than a column.",

            "An SDRF reserved word - \"not available\", \"not applicable\" - is a real value " +
            "meaning the experiment stated an absence. It is NOT the same as a column the document " +
            "does not have, and the two must not be collapsed.",
        };

        int shortRows = rows.Count(r => r.Cells.Count < document.Header.Count);
        if (shortRows > 0)
        {
            caveats.Add(
                $"{shortRows} of {rows.Count} rows are SHORT: they carry fewer cells than " +
                $"column_names has entries ({document.Header.Count}). Real files are ragged and " +
                "mzLib preserves that rather than padding, so the file round-trips. Index by " +
                "position and treat a missing position as absent, not as empty.");
        }

        return caveats;
    }

    private static List<string> PoolCaveats(SdrfCollection collection, bool labelsWereSupplied)
    {
        var caveats = new List<string>
        {
        "Every cell a source document did not have is filled with the reserved word \"not " +
        "available\". Those words are therefore partly an artefact of pooling rather than " +
        "something an experiment said, and a fill-rate computed over this table is not the fill " +
        "rate of the originals.",

        $"Provenance is in the '{SdrfCollection.SourceDocumentColumn}' column. It is stamped only " +
        "where a source document had none, so pooling an already-pooled table keeps the inner " +
        "labels rather than overwriting them.",

        "This is an ANALYSIS table, not something to deposit. source name + assay name + " +
        "comment[label] is unique within one document, but two experiments may both have a " +
        "\"Sample 1\", so the pooled table will usually violate SDRF's uniqueness rule. Use the " +
        "source-document column as part of any key.",

        $"Columns are the UNION over all {collection.Count} documents, ordered by SDRF's own block " +
        "structure, and a name that repeats is carried at the highest multiplicity any single " +
        "document used, so no values are dropped.",
        };

        if (!labelsWereSupplied)
        {
            caveats.Add(
                "No labels were supplied, so provenance fell back to mzLib's default of " +
                "containing-folder/file-stem - which depends on WHERE THE FILES SIT. Pooling the " +
                "same documents from a different directory produces different provenance values, " +
                "so this table is not reproducible on another machine. Supply a label per document " +
                "for anything you intend to keep.");
        }

        return caveats;
    }

    /// <summary>Reads and validates <c>--offset</c> and <c>--limit</c>, matching the readers verbs.</summary>
    private static (int Offset, int Limit) WindowFrom(Program.Arguments arguments)
    {
        int offset = arguments.OptionalInt("offset", 0);
        if (offset < 0)
            throw new Program.UsageException($"Option --offset must be zero or greater; got {offset}.");

        bool limited = arguments.WasProvided("limit");
        int limit = arguments.OptionalInt("limit", int.MaxValue);
        if (limited && limit < 0)
            throw new Program.UsageException($"Option --limit must be zero or greater; got {limit}.");

        return (offset, limit);
    }

    /// <summary>The selected window, and whether anything was left outside it.</summary>
    private static IReadOnlyList<T> Slice<T>(IReadOnlyList<T> all, int offset, int limit, out bool truncated)
    {
        int start = Math.Min(offset, all.Count);
        int count = (int)Math.Min((long)limit, all.Count - start);

        var selected = new List<T>(count);
        for (int i = 0; i < count; i++)
            selected.Add(all[start + i]);

        // "Were any rows left behind", by either the limit or the offset - deliberately not
        // `offset + count < all.Count`, which reports a complete answer for an offset past the end.
        // Same rule as the readers verbs.
        truncated = selected.Count < all.Count;
        return selected;
    }
}
