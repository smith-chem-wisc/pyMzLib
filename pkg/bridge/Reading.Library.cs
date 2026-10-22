using Omics.BioPolymer;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace MzLibBridge;

/// <summary>
/// <c>library read</c>: an mzLib binary spectral library (<c>.msl</c>) as one row per precursor, or
/// with <c>--fragments</c> one row per fragment ion.
/// </summary>
/// <remarks>
/// <para>
/// <b>A reader verb in another module.</b> A <c>.msl</c> file is not one of mzLib's
/// <c>SupportedFileType</c>s, so <c>readers</c> cannot open it; mzLib reads it with
/// <see cref="MslLibrary"/>, outside <c>FileReader</c>. It is still a table read from a file, so it
/// runs through the same engine as every reader verb (<see cref="RunTable"/>): <c>--path</c> or
/// <c>--paths-stdin</c>, <c>--offset</c>/<c>--limit</c>/<c>--out</c>, <c>--threads</c> and
/// <c>--on-error</c> behave exactly as they do there (bridge <c>design/BULK.md</c>).
/// </para>
/// <para>
/// <b>Opening is mzLib's decision.</b> <see cref="MslFileTypeHandler.Open"/> loads a library whole
/// when it is under <see cref="MslFileTypeHandler.DefaultIndexOnlyThresholdBytes"/> (1 GiB at the
/// pin) and index-only at or above it, where fragment blocks stay on disk. The bridge passes
/// <c>--index-only-threshold</c> to it when given and nothing otherwise. A compressed library always
/// loads whole (<see cref="MslLibrary.LoadIndexOnly"/> falls back), which <c>is_index_only</c>
/// reports.
/// </para>
/// <para>
/// <b>Rows follow the precursor index, not m/z.</b> <see cref="MslLibrary.GetAllEntries"/> yields
/// entries in ascending precursor m/z and carries no index, so a fragment row could not be joined
/// back to its precursor by a key, and in index-only mode it reads every fragment block from disk.
/// The bridge instead asks <see cref="MslLibrary.GetEntry(int, bool)"/> for each precursor index in
/// turn, which is the metadata-only loop mzLib's own documentation of <c>GetAllEntries</c>
/// prescribes. Rows therefore come in the library's stored order, the same in both load modes, and
/// <c>precursor_index</c> is the key <c>GetEntry</c> takes.
/// </para>
/// <para>
/// <b>Sentinels mzLib documents cross as null</b>, each named in the spec: an ion mobility of 0.0
/// (<see cref="MslLibraryEntry.IonMobility"/>), a NaN q-value (<see cref="MslLibraryEntry.QValue"/>),
/// an NCE of 0 (<see cref="MslLibraryEntry.Nce"/>), an empty protein accession, name or gene, a
/// terminal ion's secondary fragment number of 0 (<see cref="MslFragmentIon.SecondaryFragmentNumber"/>),
/// and a residue position of 0, which is what mzLib reads back for every custom-loss fragment.
/// </para>
/// </remarks>
internal static partial class Reading
{
    /// <summary>
    /// <c>library read --path FILE [--fragments] [--index-only-threshold BYTES] [--offset N]
    /// [--limit N] [--out FILE.tsv]</c>, or the <c>--paths-stdin</c> form.
    /// </summary>
    public static object ReadLibrary(Program.Arguments arguments)
    {
        // A flag. A value after it ("--fragments true") would otherwise land among the named
        // options, the flag would read as absent, and the caller would get a table without the
        // fragments they asked for.
        if (arguments.Optional("fragments") is { } stray)
            throw new Program.UsageException(
                $"--fragments takes no value (got '{stray}'); give --fragments alone to read one row per fragment ion.");
        bool fragments = arguments.Flag("fragments");

        long? threshold = IndexOnlyThreshold(arguments);
        var columns = fragments ? LibraryFragmentColumns : LibraryPrecursorColumns;

        var verb = new TableVerb(
            "library read",
            (path, window) => ReadLibraryFile(path, window, fragments, threshold),
            columns.Select(c => c.Name).ToList(),
            LibraryExtraKeys,
            ReportsRows: true);

        return RunTable(arguments, verb, new Dictionary<string, object?>
        {
            ["fragments_included"] = fragments,
            // The threshold mzLib applied: the caller's, or its own default when none was given.
            ["index_only_threshold"] = threshold ?? MslFileTypeHandler.DefaultIndexOnlyThresholdBytes,
        });
    }

    /// <summary><c>--index-only-threshold</c>, in bytes, or null when the caller gave none.</summary>
    private static long? IndexOnlyThreshold(Program.Arguments arguments)
    {
        if (!arguments.WasProvided("index-only-threshold"))
            return null;
        RequireValueIfProvided(arguments, "index-only-threshold");
        string raw = arguments.Optional("index-only-threshold")!;
        if (!long.TryParse(raw, System.Globalization.NumberStyles.Integer,
                System.Globalization.CultureInfo.InvariantCulture, out long bytes))
            throw new Program.UsageException(
                $"Option --index-only-threshold must be a whole number of bytes; got '{raw}'.");
        if (bytes < 0)
            throw new Program.UsageException(
                $"Option --index-only-threshold must be zero or greater (bytes); got {bytes}. " +
                "0 opens every library index-only.");
        return bytes;
    }

    /// <summary>The per-file keys this verb adds to the common block, in wire order.</summary>
    private static readonly string[] LibraryExtraKeys =
    [
        "format_version", "is_compressed", "is_index_only", "target_count", "decoy_count",
        "proteoform_count", "min_precursor_mz", "max_precursor_mz",
    ];

    /// <summary>One table row: a precursor, and one of its fragment ions when fragments were asked for.</summary>
    /// <remarks><c>Fragment</c> is null on a precursor-only row, and on the single row of a
    /// precursor that has no fragment ions.</remarks>
    private sealed record LibraryRow(int PrecursorIndex, MslLibraryEntry Entry, MslFragmentIon? Fragment);

    /// <summary>The precursor columns, under <see cref="MslLibraryEntry"/>'s own names.</summary>
    private static IReadOnlyList<Column<LibraryRow>> LibraryPrecursorColumns { get; } =
    [
        // The zero-based index MslLibrary.GetEntry takes (MslPrecursorIndexEntry.PrecursorIdx).
        new("precursor_index", r => r.PrecursorIndex),
        new("full_sequence", r => r.Entry.FullSequence),
        new("base_sequence", r => r.Entry.BaseSequence),
        new("precursor_mz", r => r.Entry.PrecursorMz),
        new("charge_state", r => r.Entry.ChargeState),
        new("retention_time", r => r.Entry.RetentionTime),
        // MslLibraryEntry.RtIsCalibrated: false means RetentionTime holds iRT, true run minutes.
        new("retention_time_unit", r => r.Entry.RtIsCalibrated ? "minutes" : "iRT"),
        // MslLibraryEntry.IonMobility: "0.0 means the value is not available". It is read from a
        // stored float32, so an exact comparison is the documented sentinel, not a derived value.
        new("ion_mobility", r => r.Entry.IonMobility != 0.0 ? r.Entry.IonMobility : null),
        new("is_decoy", r => r.Entry.IsDecoy),
        new("is_entrapment", r => r.Entry.IsEntrapment),
        // mzLib's one label writer. The .msl format has no contaminant flag, so C never appears.
        new("decoy_contaminant_target", r => DecoyContaminantTargetLabel.For(r.Entry.IsDecoy, false, r.Entry.IsEntrapment)),
        // "Empty string or null when absent" (MslLibraryEntry.ProteinAccession, ProteinName, GeneName).
        new("protein_accession", r => NullIfEmpty(r.Entry.ProteinAccession)),
        new("protein_name", r => NullIfEmpty(r.Entry.ProteinName)),
        new("gene_name", r => NullIfEmpty(r.Entry.GeneName)),
        new("is_proteotypic", r => r.Entry.IsProteotypic),
        // MslLibraryEntry.QValue: "float.NaN when no q-value is available (e.g. purely predicted library)".
        new("q_value", r => float.IsNaN(r.Entry.QValue) ? null : r.Entry.QValue),
        new("elution_group_id", r => r.Entry.ElutionGroupId),
        new("molecule_type", r => r.Entry.MoleculeType.ToString()),
        // MslLibraryEntry.Source (MslFormat.SourceType). Not "source", which a bulk table would
        // set beside source_index and source_path.
        new("source_type", r => r.Entry.Source.ToString()),
        new("dissociation_type", r => r.Entry.DissociationType.ToString()),
        // MslLibraryEntry.Nce: "0 = unknown or not applicable".
        new("nce", r => r.Entry.Nce != 0 ? r.Entry.Nce : null),
    ];

    /// <summary>The precursor columns, then one fragment ion's, under <see cref="MslFragmentIon"/>'s names.</summary>
    private static IReadOnlyList<Column<LibraryRow>> LibraryFragmentColumns { get; } =
    [
        .. LibraryPrecursorColumns,
        // Counted from the fragments GetEntry(i, includeFragments: true) read, so 0 is a real zero:
        // the precursor's single row then has every fragment cell null.
        new("fragment_count", r => r.Entry.MatchedFragmentIons.Count),
        new("fragment_mz", r => r.Fragment?.Mz),
        new("fragment_intensity", r => r.Fragment?.Intensity),
        new("product_type", r => r.Fragment?.ProductType.ToString()),
        // Null for every terminal ion (MslFragmentIon.SecondaryProductType).
        new("secondary_product_type", r => r.Fragment?.SecondaryProductType?.ToString()),
        new("fragment_number", r => r.Fragment?.FragmentNumber),
        // "0 for all terminal ions" (MslFragmentIon.SecondaryFragmentNumber).
        new("secondary_fragment_number", r => r.Fragment is { IsInternalFragment: true } f ? f.SecondaryFragmentNumber : null),
        // 0 is what MslReader.ConvertFragment reads back for a custom-loss fragment ("not available
        // for custom-loss fragments"), and no terminal ion has residue 0, since positions are 1-based.
        new("residue_position", r => r.Fragment is { ResiduePosition: not 0 } f ? f.ResiduePosition : null),
        new("fragment_charge", r => r.Fragment?.Charge),
        // A loss of 0.0 Da is the value of "no loss", not a missing one, so it crosses as 0.
        new("neutral_loss", r => r.Fragment?.NeutralLoss),
        new("is_internal", r => r.Fragment?.IsInternalFragment),
        new("is_diagnostic", r => r.Fragment?.IsDiagnosticIon),
        new("exclude_from_quant", r => r.Fragment?.ExcludeFromQuant),
        new("annotation", r => r.Fragment?.Annotation),
    ];

    private static string? NullIfEmpty(string? text) => string.IsNullOrEmpty(text) ? null : text;

    /// <summary>One file through <c>library read</c>.</summary>
    private static FileTable ReadLibraryFile(string path, Window window, bool fragments, long? threshold)
    {
        if (!File.Exists(path))
            throw new Program.UsageException(
                Directory.Exists(path)
                    ? $"'{path}' is a directory; library read reads one .msl file."
                    : $"File not found: '{path}'.");
        // mzLib's own test for the format. Open does not check the extension itself ("call
        // IsMslFile first if extension validation is required").
        if (!MslFileTypeHandler.IsMslFile(path))
            throw new Program.UsageException(
                $"library read reads an mzLib spectral library, a .msl file (MslFileTypeHandler.IsMslFile); " +
                $"'{path}' is not one. Text libraries (.msp) are not read by this verb.");

        MslLibrary library = threshold is null
            ? MslFileTypeHandler.Open(path)
            : MslFileTypeHandler.Open(path, threshold.Value);

        try
        {
            return fragments
                ? LibraryFragmentTable(path, window, library)
                : LibraryPrecursorTable(path, window, library);
        }
        catch
        {
            // An index-only library holds the file open; nothing else will release it.
            library.Dispose();
            throw;
        }
    }

    /// <summary>The window over precursor indices: the first index, how many, and whether any were left out.</summary>
    private static (int Start, int Count, bool Truncated) LibraryWindow(Window window, int total)
    {
        int start = Math.Min(window.Offset, total);
        int count = (int)Math.Min((long)window.Limit, total - start);
        // Window.Apply's rule: whether any record was left behind, by either limit or offset.
        return (start, count, count < total);
    }

    /// <summary>The entry at a precursor index, or an error naming it: mzLib returns null only out of range.</summary>
    private static MslLibraryEntry LibraryEntryAt(MslLibrary library, int index, bool includeFragments) =>
        library.GetEntry(index, includeFragments)
        ?? throw new InvalidOperationException(
            $"mzLib returned no entry for precursor {index} of {library.PrecursorCount} (MslLibrary.GetEntry).");

    /// <summary>One row per precursor, read without fragments.</summary>
    private static FileTable LibraryPrecursorTable(string path, Window window, MslLibrary library)
    {
        (int start, int count, bool truncated) = LibraryWindow(window, library.PrecursorCount);

        // includeFragments: false reads nothing from disk in index-only mode. In full-load mode the
        // fragments come along anyway, and mzLib says not to rely on either, so none is read here.
        var rows = new List<LibraryRow>(count);
        for (int k = 0; k < count; k++)
            rows.Add(new LibraryRow(start + k, LibraryEntryAt(library, start + k, includeFragments: false), null));

        var columns = LibraryPrecursorColumns;
        return LibraryTable(path, window, library, columns, FailedColumns(columns, rows), fragments: false,
            rowCount: rows.Count, recordsReturned: null, row: i => Cells(columns, rows[i]), truncated);
    }

    /// <summary>One row per fragment ion, the precursor's columns repeated on each.</summary>
    /// <remarks>
    /// The engine needs the row count before the first row, and an index-only library cannot say
    /// how many fragments a precursor has without reading them (mzLib keeps the stored count
    /// private). So the fragments are read twice: once to count, and once per row as it is
    /// written. Nothing is held between the two but one count per precursor, so an index-only read
    /// written to <c>--out</c> stays within mzLib's fragment cache however large the library is.
    /// </remarks>
    private static FileTable LibraryFragmentTable(string path, Window window, MslLibrary library)
    {
        (int start, int count, bool truncated) = LibraryWindow(window, library.PrecursorCount);
        var columns = LibraryFragmentColumns;

        // firstRow[k] is the first row of the k-th selected precursor; every precursor has at least
        // one row, so the array is strictly increasing.
        var firstRow = new int[count + 1];
        var failed = new SortedSet<string>(StringComparer.Ordinal);
        for (int k = 0; k < count; k++)
        {
            MslLibraryEntry entry = LibraryEntryAt(library, start + k, includeFragments: true);
            List<LibraryRow> rows = RowsOf(start + k, entry);
            failed.UnionWith(FailedColumns(columns, rows));
            firstRow[k + 1] = checked(firstRow[k] + rows.Count);
        }

        int cachedPrecursor = -1;
        MslLibraryEntry? cached = null;
        object?[] Row(int i)
        {
            int k = Array.BinarySearch(firstRow, 0, count, i);
            if (k < 0)
                k = ~k - 1;
            if (k != cachedPrecursor)
            {
                cached = LibraryEntryAt(library, start + k, includeFragments: true);
                cachedPrecursor = k;
            }

            int ordinal = i - firstRow[k];
            MslFragmentIon? fragment = cached!.MatchedFragmentIons.Count == 0 ? null : cached.MatchedFragmentIons[ordinal];
            return Cells(columns, new LibraryRow(start + k, cached, fragment));
        }

        return LibraryTable(path, window, library, columns, failed.ToList(), fragments: true,
            rowCount: firstRow[count], recordsReturned: count, row: Row, truncated);
    }

    /// <summary>A precursor's rows: one per fragment ion, or one with no fragment when it has none.</summary>
    private static List<LibraryRow> RowsOf(int index, MslLibraryEntry entry) =>
        entry.MatchedFragmentIons.Count == 0
            ? [new LibraryRow(index, entry, null)]
            : entry.MatchedFragmentIons.Select(f => new LibraryRow(index, entry, f)).ToList();

    private static FileTable LibraryTable(
        string path,
        Window window,
        MslLibrary library,
        IReadOnlyList<Column<LibraryRow>> columns,
        IReadOnlyList<string> failed,
        bool fragments,
        int rowCount,
        int? recordsReturned,
        Func<int, object?[]> row,
        bool truncated)
    {
        List<string> names = columns.Select(c => c.Name).ToList();
        return new FileTable
        {
            Block = Block(
                path,
                // A .msl file is not an mzLib SupportedFileType, so there is no type name to give.
                fileType: null,
                reader: nameof(MslLibrary),
                recordCount: library.PrecursorCount,
                // Binary, and validated whole (magic, version, footer, count, CRC-32) before any
                // entry is returned: mzLib reads every precursor or throws.
                rowsNotRead: null,
                // The unit varies by precursor; it is the retention_time_unit column.
                retentionTimeUnit: null,
                caveats: fragments ? [.. LibraryCaveats, .. LibraryFragmentCaveats] : LibraryCaveats,
                columnNames: names,
                absent: AbsentLibraryFields(library),
                failed: failed,
                excluded: [],
                ("format_version", library.Header.FormatVersion),
                ("is_compressed", library.IsCompressed),
                ("is_index_only", library.IsIndexOnly),
                ("target_count", library.TargetCount),
                ("decoy_count", library.DecoyCount),
                ("proteoform_count", library.ProteoformCount),
                // An empty library has no range. mzLib's statistics hold +/- infinity then (its
                // property docs say 0); either way there is no m/z to report.
                ("min_precursor_mz", FiniteOrNull(library.MinPrecursorMz)),
                ("max_precursor_mz", FiniteOrNull(library.MaxPrecursorMz))),
            ColumnNames = names,
            RowCount = rowCount,
            RecordsReturned = recordsReturned,
            Row = row,
            Offset = window.Offset,
            Truncated = truncated,
            Release = library.Dispose,
        };
    }

    private static float? FiniteOrNull(float value) => float.IsFinite(value) ? value : null;

    /// <summary>
    /// The columns this library has no source for, from mzLib's file-level flags: no precursor
    /// carries an ion mobility, there is no protein table, or there are no gene names.
    /// </summary>
    private static List<string> AbsentLibraryFields(MslLibrary library)
    {
        int flags = library.Header.FileFlags;
        var absent = new List<string>();
        if ((flags & MslFormat.FileFlagHasIonMobility) == 0)
            absent.Add("ion_mobility");
        if ((flags & MslFormat.FileFlagHasProteinData) == 0)
        {
            absent.Add("protein_accession");
            absent.Add("protein_name");
        }
        if ((flags & MslFormat.FileFlagHasGeneData) == 0)
            absent.Add("gene_name");
        return absent;
    }

    private static List<string> LibraryCaveats =>
    [
        "Rows are in the library's stored order, by precursor_index, not by m/z: mzLib's own " +
        "MslLibrary.GetAllEntries enumerates in ascending precursor m/z instead. Sort on precursor_mz " +
        "for that order.",
        "retention_time is iRT for some precursors and run minutes for others: read retention_time_unit " +
        "on each row (MslLibraryEntry.RtIsCalibrated). A library can hold both, so never plot or " +
        "compare the column without it. There is no conversion: iRT becomes minutes only through a " +
        "calibration against a run.",
        "The library is unfiltered: decoys (D, ED) and entrapment precursors (ET, ED) are rows. Filter " +
        "on is_decoy, or on decoy_contaminant_target, before counting targets.",
        "A retention_time of 0 may mean the source had none: MslLibraryEntry.FromLibrarySpectrum, the " +
        "path mzLib converts a text library through, writes RetentionTime ?? 0.0, and mzLib documents " +
        "no sentinel for it, so 0 is passed through as a value.",
        "nce is a whole number: the format stores NCE x 10 and MslReader reads it back with integer " +
        "division, so an NCE of 27.5 reads as 27.",
        "A library converted from a text library (.msp) has no NCE, dissociation type, q-value, ion " +
        "mobility or protein: LibrarySpectrum has no field for them, and MslLibraryEntry.FromLibrarySpectrum " +
        "fills mzLib's not-available defaults, so nce and q_value are null and dissociation_type is Unknown.",
        "precursor_mz, retention_time and ion_mobility are stored as 32-bit floats and widened to double " +
        "by mzLib, so 465.74 reads back as 465.739990234375. The extra digits are the float, not a " +
        "measurement; round them for display.",
        "offset and limit count PRECURSORS. With fragments, each precursor gives one row per fragment ion.",
    ];

    private static List<string> LibraryFragmentCaveats =>
    [
        "residue_position is null for every fragment with a custom neutral loss: the .msl format reuses " +
        "that field to index its loss table, and MslReader reads it back as 0.",
        "fragment_intensity is relative within its precursor, 0 to 1: mzLib's writer scales each " +
        "precursor's most intense fragment to 1 (MslWriter). It is not comparable across precursors.",
    ];
}
