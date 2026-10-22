using Readers;

namespace MzLibBridge;

/// <summary>
/// <c>sdrf restrict</c>: an SDRF restricted to the files one search read, by mzLib's
/// <see cref="SdrfSearchScope.Restrict"/>, with everything the join could not place.
/// </summary>
/// <remarks>
/// <para>
/// <b>The join is mzLib's, all of it.</b> A search reads <c>X-calib.mzML</c> while the SDRF names the
/// acquired <c>X.raw</c>; mzLib joins them by stem, case-insensitively, after removing MetaMorpheus's
/// derivative suffixes, and writes the searched name into <c>comment[searched data file]</c>. The
/// bridge parses stdin, calls <see cref="SdrfSearchScope.Restrict"/> once, and shapes the restricted
/// <see cref="SdrfDocument"/> exactly as <c>sdrf read</c> shapes a document: row-major raw cells.
/// </para>
/// <para>
/// <b>It pairs with the design verbs.</b> <c>sdrf design</c> and <c>sdrf isobaric-design</c> match
/// searched files by EXACT name, by design (mzLib leaves the stem join to this type). A document
/// restricted here names the searched files in <c>comment[searched data file]</c>, which both
/// design readers key on first (MAP-12), so restrict, then design.
/// </para>
/// <para>
/// <b><c>--out</c> follows <c>sdrf pool --out</c></b> as bridge thread 010 ruled it: <c>.tsv</c>
/// only, checked before any read (PYB-4), and a missing folder created (PYB-5). The write is mzLib's
/// <see cref="SdrfDocument.WriteResults"/>, which adds <c>.sdrf.tsv</c> to a name that does not end
/// in it; <c>written.path</c> reports the file mzLib wrote, by mzLib's own test
/// (<see cref="ResultFile{TResult}.CanRead"/>), never the name as typed.
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    /// <summary>
    /// <c>sdrf restrict --path FILE [--out FILE.tsv] [--limit N] [--offset N]</c>, with the searched
    /// files on stdin, one per line: <c>searched</c>, or <c>searched&lt;TAB&gt;acquired</c> when the
    /// caller knows the acquired file a searched one came from.
    /// </summary>
    public static object Restrict(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        string? outputPath = arguments.Optional("out");
        if (arguments.WasProvided("out") && string.IsNullOrWhiteSpace(outputPath))
            throw new Program.UsageException(
                "Option --out needs a file path ending in .tsv, e.g. 'restricted.sdrf.tsv'.");
        Program.RequireTsvOutput(outputPath);

        (int offset, int limit) = WindowFrom(arguments);

        if (!File.Exists(path))
            throw new Program.UsageException($"SDRF file not found: '{path}'.");

        (List<string> searched, Dictionary<string, string>? acquiredNameOf) = ReadSearchedFiles();

        var document = new SdrfDocument(path);
        // Throws ArgumentException for an SDRF with no comment[data file] column: a fact about the
        // document, so it crosses as mzLib's own failure (correctness), not as a usage error.
        SdrfSearchScoping scoping = SdrfSearchScope.Restrict(document, searched, acquiredNameOf);
        SdrfDocument restricted = scoping.Document;

        IReadOnlyList<SdrfRow> all = restricted.Results;
        IReadOnlyList<SdrfRow> selected = Slice(all, offset, limit, out bool truncated);

        object? written = null;
        if (outputPath is not null)
        {
            // The WHOLE restricted document, never the windowed slice, as sdrf pool writes.
            CreateOutputFolder(outputPath);
            restricted.WriteResults(outputPath);
            written = new { path = WrittenSdrfPath(restricted, outputPath), row_count = all.Count };
        }

        return new
        {
            path,
            searched_files = searched,
            acquired_name_count = acquiredNameOf?.Count ?? 0,
            input_row_count = document.Results.Count,
            column_names = restricted.Header.ToList(),
            row_count = all.Count,
            returned_count = selected.Count,
            offset,
            truncated,
            rows = selected.Select(r => r.Cells.ToList()).ToList(),
            searched_without_row = scoping.SearchedWithoutRow.ToList(),
            dropped_data_files = scoping.DroppedDataFiles.ToList(),
            ambiguous = scoping.Ambiguous.ToList(),
            shared_stem = scoping.SharedStem.ToList(),
            written,
            caveats = RestrictCaveats(),
        };
    }

    /// <summary>
    /// The searched files from stdin, in order, and the acquired name given for any of them (null
    /// when none was, which is mzLib's default: the suffix rule for every file).
    /// </summary>
    private static (List<string> Searched, Dictionary<string, string>? AcquiredNameOf) ReadSearchedFiles()
    {
        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No searched files were given on stdin. Supply one per line, as a path or a file " +
                "name, optionally followed by a tab and the acquired file it came from, e.g. " +
                "'X-calib.mzML' or 'X-calib.mzML<TAB>X.raw'.");

        var searched = new List<string>(lines.Count);
        // Case-insensitive, as mzLib collapses searched files that differ only in case
        // (SdrfSearchScope.Restrict keeps the first spelling and looks it up here).
        var acquiredNameOf = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        foreach (string line in lines)
        {
            string[] fields = line.Split('\t');
            if (fields.Length > 2)
                throw new Program.UsageException(
                    $"A stdin line has {fields.Length} tab-separated fields: '{line}'. Give " +
                    "'searched' or 'searched<TAB>acquired'.");
            string file = fields[0].Trim();
            if (file.Length == 0)
                throw new Program.UsageException($"A stdin line has an empty searched file: '{line}'.");
            searched.Add(file);

            string acquired = fields.Length == 2 ? fields[1].Trim() : string.Empty;
            if (acquired.Length == 0)
                continue;
            if (acquiredNameOf.TryGetValue(file, out string? earlier) &&
                !string.Equals(earlier, acquired, StringComparison.Ordinal))
                throw new Program.UsageException(
                    $"Searched file '{file}' is given two acquired files, '{earlier}' and '{acquired}'. " +
                    "A searched file comes from one acquisition.");
            acquiredNameOf[file] = acquired;
        }

        return (searched, acquiredNameOf.Count == 0 ? null : acquiredNameOf);
    }

    /// <summary>
    /// Creates the folder an <c>--out</c> file goes into when it is missing (bridge thread 010,
    /// PYB-5: the caller named the path, so creating its folder writes nowhere they did not ask for).
    /// mzLib's <see cref="SdrfDocument.WriteResults"/> does not create it.
    /// </summary>
    private static void CreateOutputFolder(string outputPath)
    {
        string? folder = Path.GetDirectoryName(Path.GetFullPath(outputPath));
        if (!string.IsNullOrEmpty(folder))
            Directory.CreateDirectory(folder);
    }

    /// <summary>
    /// The file <see cref="SdrfDocument.WriteResults"/> wrote for <paramref name="outputPath"/>: it
    /// adds <c>.sdrf.tsv</c> to a name it would not read back as an SDRF (its own
    /// <see cref="ResultFile{TResult}.CanRead"/>), so <c>x.tsv</c> becomes <c>x.tsv.sdrf.tsv</c>.
    /// sdrf pool and sdrf restrict both report this, never the name as typed (SDRF-B3 asks mzLib to
    /// return the path itself).
    /// </summary>
    internal static string WrittenSdrfPath(SdrfDocument document, string outputPath) =>
        document.CanRead(outputPath) ? outputPath : outputPath + SupportedFileType.Sdrf.GetFileExtension();

    private static List<string> RestrictCaveats() => new()
    {
        "The join is by STEM, case-insensitively, after mzLib removes MetaMorpheus's derivative " +
        "suffixes (-calib, -averaged, in any order and number), so X-calib.mzML finds the row " +
        "naming X.raw. A searched<TAB>acquired line overrides that rule for its file. A row whose " +
        "stem joins nothing is tried once more by an earlier search's comment[searched data file].",

        "Kept rows are carried cell for cell: no replicate or fraction is renumbered. mzLib adds " +
        "comment[searched data file] (the searched file's name, no folder) directly after " +
        "comment[data file], and replaces one an earlier search wrote.",

        "A kept row shorter than the header is padded with empty cells, so this document is never " +
        "ragged, unlike what sdrf read returns for the same file.",

        "ambiguous: two or more searched files claimed one acquisition; the first in name order " +
        "kept it, and the others get no row and are not in searched_without_row. shared_stem: one " +
        "searched file joined rows naming two or more DIFFERENT acquired files (X.raw and X.mzML); " +
        "every such row is kept, so the document can describe one run twice, and only the " +
        "depositor can say which is right. A row naming no file is dropped without being listed " +
        "in dropped_data_files.",

        "The design verbs match searched files by EXACT name. Restrict first: sdrf design and sdrf " +
        "isobaric-design key on comment[searched data file] when a document has it, so the " +
        "restricted document joins the search's own file names.",

        "--out writes the WHOLE restricted document, whatever limit and offset are, through mzLib's " +
        "SdrfDocument.WriteResults: atomically, with CRLF line endings whatever the input used. " +
        "mzLib adds .sdrf.tsv to a name that does not already end in it, so written.path is the " +
        "file actually written.",
    };
}
