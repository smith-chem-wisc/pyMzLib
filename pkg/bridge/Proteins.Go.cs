using System.Security.Cryptography;
using Readers;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.GeneOntology;

namespace MzLibBridge;

/// <summary>
/// Gene Ontology on protein groups: <c>proteins annotate-go</c> and <c>proteins update-go</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>The bridge adds no biology here either.</b> Propagation up <c>is_a</c> and <c>part_of</c> is
/// <see cref="GeneOntologyGraph"/>, the per-group rows are <see cref="GoGroupAnnotator"/> (mzLib #1353),
/// reading a stored MetaMorpheus table is <see cref="ProteinGroupFromTsvExtensions.ToGoAnnotationGroups"/>
/// (#1366), and categories are <see cref="GoCategoryResolver"/>. The column names are mzLib's own
/// <see cref="GoAnnotationTsv.Schema"/> headers, the file a caller asks for is written by
/// <see cref="GoAnnotationTsv.Write"/> itself, and the provenance header on the wire is the one that
/// writer produced, read back. So a table read here, a file written here, and a file mzLib writes
/// inside a search are the same thing.
/// </para>
/// <para>
/// <b>Nothing is downloaded behind the caller's back.</b> mzLib's <see cref="Loaders.LoadGeneOntology"/>
/// fetches go.obo when the file is absent, from a PURL that moves with every GO release, so which
/// release a run used would depend on the day it first ran. <c>annotate-go</c> therefore reads only a
/// file that exists, and fetching is its own verb, <c>update-go</c>, which a caller runs on purpose.
/// </para>
/// </remarks>
internal static partial class Proteins
{
    /// <summary>
    /// The annotation table's columns: exactly mzLib's <see cref="GoAnnotationTsv.Schema"/> headers, in its
    /// order. A test holds them equal, so the wire and mzLib's own file share one vocabulary.
    /// </summary>
    internal static readonly string[] GoAnnotationColumns =
    {
        "protein_group", "accession_used", "accession_direct", "accession_inherited", "go_id", "go_name",
        "aspect", "evidence", "evidence_by_member", "inherited", "propagated", "n_members", "n_with",
        "entrapment_members", "annotation_status", "q_value", "go_release", "go_obo_sha256",
        "annotation_db_sha256",
    };

    /// <summary>The category table's columns, as <see cref="GoCategoryTsv.Write"/> writes them.</summary>
    internal static readonly string[] GoCategoryColumns = { "go_id", "category", "subcategory" };

    // ---- proteins annotate-go -----------------------------------------------------------------

    /// <summary>
    /// <c>proteins annotate-go --groups FILE --database FILE --go-obo FILE [--skip-unknown-go-ids]
    /// [--category-map FILE] [--out FILE.tsv] [--categories-out FILE.tsv] [--limit N] [--offset N]</c> —
    /// one row per (protein group, GO term) that any member of the group holds, directly or through an
    /// ancestor.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Every argument is checked before any file is read</b> (bridge thread 008, PYB-1): an
    /// <c>--out</c> that is not <c>.tsv</c> fails in milliseconds, not after a proteome has loaded.
    /// </para>
    /// <para>
    /// <b>One annotation database.</b> <see cref="GoGroupAnnotator"/> stamps every row with the sha256 of
    /// <i>the</i> database its terms came from, so two databases have no honest single value. Use the
    /// UniProt XML the search used; a FASTA is refused, because it carries no GO and would make every
    /// group read <c>no_go_terms</c>, which is a statement about the file format, not about the proteins.
    /// </para>
    /// <para>
    /// <b>The rows on the wire are a window; the file is everything.</b> As <c>sdrf pool</c> does,
    /// <c>--limit</c>/<c>--offset</c> select the rows returned, and <c>--out</c> always receives the whole
    /// table. A large run (87,513 rows from one 715-group table, mzLib #1353) belongs in <c>--out</c> with
    /// <c>--limit 0</c>.
    /// </para>
    /// </remarks>
    public static object AnnotateGo(Program.Arguments arguments)
    {
        foreach (string name in new[] { "groups", "database", "go-obo", "category-map", "out", "categories-out", "limit", "offset" })
            RequireValueIfProvided(arguments, name);

        string groupsPath = ExistingFile(arguments.Required("groups"), "Protein-group table");
        string databasePath = ExistingFile(arguments.Required("database"), "Annotation database");
        string goOboPath = arguments.Required("go-obo");
        if (!File.Exists(goOboPath))
            throw new Program.UsageException(
                $"go.obo not found: '{goOboPath}'. Nothing is downloaded here, so a run cannot silently pick up " +
                "whichever release the PURL serves today. Fetch one on purpose with 'proteins update-go --go-obo " +
                "<path>', then keep that file: it pins the release every row is computed against.");
        goOboPath = Path.GetFullPath(goOboPath);

        if (FormatOf(databasePath) == DatabaseFormat.Fasta)
            throw new Program.UsageException(
                $"'{Path.GetFileName(databasePath)}' is a FASTA, which carries no GO terms, so every group would read " +
                "no_go_terms whatever the proteins are. Give the UniProt XML (.xml or .xml.gz) of the proteome the " +
                "search used.");

        string? mapPath = arguments.Optional("category-map");
        if (mapPath is not null)
            mapPath = ExistingFile(mapPath, "Category map");

        string? outPath = TsvOutput(arguments, "out");
        string? categoriesOutPath = TsvOutput(arguments, "categories-out");
        if (categoriesOutPath is not null && mapPath is null)
            throw new Program.UsageException("Option --categories-out writes the category table, so it needs --category-map.");

        var inputs = new[] { groupsPath, databasePath, goOboPath, mapPath };
        foreach ((string? output, string option) in new[] { (outPath, "out"), (categoriesOutPath, "categories-out") })
        {
            if (output is not null && inputs.Any(i => i is not null && SamePath(i, output)))
                throw new Program.UsageException($"Option --{option} names an input file; writing there would overwrite it.");
        }
        if (outPath is not null && categoriesOutPath is not null && SamePath(outPath, categoriesOutPath))
            throw new Program.UsageException("Options --out and --categories-out name the same file.");

        (int offset, int limit) = GoWindowFrom(arguments);
        bool skipUnknown = arguments.Flag("skip-unknown-go-ids");

        // ---- read ----

        GeneOntologyGraph ontology = GeneOntologyGraph.Load(goOboPath);

        var database = new LoadedDatabase
        {
            Input = new DatabaseInput(0, databasePath, Contaminant: false),
            AbsolutePath = databasePath,
        };
        Load(database);
        database.Sha256 = DecompressedSha256(databasePath);

        var annotator = new GoGroupAnnotator(ontology, database.Proteins!, database.Sha256, skipUnknown);

        var table = new ProteinGroupFromTsvFile(groupsPath);
        List<ProteinGroupFromTsv> stored = table.Results;
        List<GoAnnotationGroup> groups = stored.ToGoAnnotationGroups().ToList();
        IReadOnlyList<GoAnnotationRow> rows = annotator.AnnotateAll(groups);
        string groupsSha256 = FileSha256(groupsPath);

        // ---- write: mzLib's own writer, to the caller's file or to memory for the header ----

        Dictionary<string, string> header;
        object? written = null;
        if (outPath is not null)
        {
            WriteAtomically(outPath, writer =>
                GoAnnotationTsv.Write(writer, rows, ontology, database.Sha256, groupsSha256, annotator.UnresolvedGoIds));
            header = ReadHeader(File.ReadLines(outPath));
            written = new { path = outPath, row_count = rows.Count };
        }
        else
        {
            using var buffer = new StringWriter();
            GoAnnotationTsv.Write(buffer, rows, ontology, database.Sha256, groupsSha256, annotator.UnresolvedGoIds);
            header = ReadHeader(Lines(buffer.ToString()));
        }

        object? categories = null;
        object? categoriesWritten = null;
        if (mapPath is not null)
        {
            GoCategoryMap map = GoCategoryMap.Load(mapPath);
            var resolver = new GoCategoryResolver(map, ontology);
            using var buffer = new StringWriter();
            GoCategoryTsv.Write(buffer, resolver, rows);
            string text = buffer.ToString();
            if (categoriesOutPath is not null)
            {
                WriteAtomically(categoriesOutPath, writer => writer.Write(text));
                categoriesWritten = new { path = categoriesOutPath };
            }
            categories = CategoryTable(map, text);
        }

        // ---- shape ----

        int start = Math.Min(offset, rows.Count);
        int count = (int)Math.Min((long)limit, rows.Count - start);
        var columns = NewColumns(GoAnnotationColumns);
        foreach (GoAnnotationRow row in rows.Skip(start).Take(count))
        {
            AddRow(columns, row.ProteinGroup, row.AccessionUsed.ToList(), row.AccessionDirect.ToList(),
                row.AccessionInherited.ToList(), row.GoId, row.GoName, AspectName(row), row.Evidence.ToList(),
                row.EvidenceByMember.ToDictionary(m => m.Key, m => m.Value.ToList(), StringComparer.Ordinal),
                row.Inherited, row.Propagated, row.NMembers, row.NWith, row.EntrapmentMembers.ToList(),
                GoAnnotationTsv.StatusName(row.Status), row.QValue, row.GoRelease, row.GoOboSha256,
                row.AnnotationDbSha256);
        }

        return new
        {
            groups_file = groupsPath,
            groups_file_sha256 = groupsSha256,
            table_row_count = stored.Count,
            decoy_group_count = stored.Count(r => r.IsDecoy),
            group_count = groups.Count,
            annotation_database = new
            {
                path = databasePath,
                file_type = database.FileType,
                reader = database.Reader,
                protein_count = database.Proteins!.Count,
                sha256 = database.Sha256,
                caveats = database.Caveats,
            },
            go = new
            {
                source_file_name = ontology.SourceFileName,
                sha256 = ontology.SourceSha256,
                release = ontology.Release,
                term_count = ontology.Count,
            },
            skip_unknown_go_ids = skipUnknown,
            unresolved_go_ids = annotator.UnresolvedGoIds.ToList(),
            header,
            row_count = rows.Count,
            returned_count = count,
            offset,
            truncated = count < rows.Count,
            caveats = GoCaveats(ontology, annotator, database),
            written,
            categories,
            categories_written = categoriesWritten,
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    // ---- proteins update-go -------------------------------------------------------------------

    /// <summary>
    /// <c>proteins update-go --go-obo FILE</c> — download the current go.obo to <c>FILE</c>, on purpose.
    /// </summary>
    /// <remarks>
    /// Wraps <see cref="Loaders.UpdateGeneOntology(string, CancellationToken)"/>: the whole file (~37 MB)
    /// is streamed from GO's PURL; a file already there is kept as a timestamped backup when the
    /// download differs, and left alone when it is the same. A transport failure leaves any existing file
    /// untouched and crosses as <c>ServiceUnavailable</c>.
    /// </remarks>
    public static object UpdateGo(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "go-obo");
        string path = Path.GetFullPath(arguments.Required("go-obo"));
        string? directory = Path.GetDirectoryName(path);
        if (directory is not null && !Directory.Exists(directory))
            throw new Program.UsageException($"The folder for --go-obo does not exist: '{directory}'.");

        bool existed = File.Exists(path);
        string? previous = existed ? FileSha256(path) : null;

        Loaders.UpdateGeneOntology(path);
        GeneOntologyGraph ontology = GeneOntologyGraph.Load(path);
        bool changed = !string.Equals(previous, ontology.SourceSha256, StringComparison.Ordinal);

        var caveats = new List<string>();
        if (existed && changed)
            caveats.Add(
                $"The file that was there (sha256 {previous}) was kept beside it as " +
                $"'{Path.GetFileName(path)}.<yyyyMMdd-HHmmss-fff>', so a run made against it can still be reproduced.");
        if (ontology.Release is null)
            caveats.Add("The downloaded file has no data-version header, so its release is unknown; its sha256 still pins it.");

        return new
        {
            go_obo_file = path,
            url = Loaders.GeneOntologyUrl,
            existed_before = existed,
            previous_sha256 = previous,
            changed,
            go = new
            {
                source_file_name = ontology.SourceFileName,
                sha256 = ontology.SourceSha256,
                release = ontology.Release,
                term_count = ontology.Count,
            },
            caveats,
        };
    }

    // ---- helpers ------------------------------------------------------------------------------

    private static string ExistingFile(string path, string what)
    {
        if (!File.Exists(path))
            throw new Program.UsageException($"{what} not found: '{path}'.");
        return Path.GetFullPath(path);
    }

    /// <summary>
    /// An output path, which must end in <c>.tsv</c> (any case): the file is mzLib's TSV, and a path that
    /// says otherwise would mislabel it (bridge thread 008, PYB-1). A path with no extension is refused
    /// too, rather than given one.
    /// </summary>
    private static string? TsvOutput(Program.Arguments arguments, string option)
    {
        string? path = arguments.Optional(option);
        if (path is null)
            return null;
        Program.RequireTsvOutput(path, option);
        string full = Path.GetFullPath(path);
        string? directory = Path.GetDirectoryName(full);
        if (directory is not null && !Directory.Exists(directory))
            throw new Program.UsageException($"The folder for --{option} does not exist: '{directory}'.");
        return full;
    }

    private static bool SamePath(string a, string b) =>
        string.Equals(Path.GetFullPath(a), Path.GetFullPath(b), StringComparison.OrdinalIgnoreCase);

    private static (int Offset, int Limit) GoWindowFrom(Program.Arguments arguments)
    {
        int offset = arguments.OptionalInt("offset", 0);
        if (offset < 0)
            throw new Program.UsageException($"Option --offset must be zero or greater; got {offset}.");
        int limit = arguments.OptionalInt("limit", int.MaxValue);
        if (limit < 0)
            throw new Program.UsageException($"Option --limit must be zero or greater; got {limit}.");
        return (offset, limit);
    }

    /// <summary>
    /// The aspect exactly as mzLib's writer spells it (<c>biological_process</c>, ...), taken from its own
    /// schema column rather than re-derived, so the wire cannot drift from the file.
    /// </summary>
    private static string? AspectName(GoAnnotationRow row) => AspectColumn.GetValue(row);

    private static readonly Omics.BioPolymerGroup.TsvColumn<GoAnnotationRow> AspectColumn =
        GoAnnotationTsv.Schema.Single(c => c.Header == "aspect");

    /// <summary>Writes through a temporary file beside the target, so a failure never leaves half a table.</summary>
    private static void WriteAtomically(string path, Action<TextWriter> write)
    {
        string temp = path + ".partial";
        try
        {
            using (var writer = new StreamWriter(temp, append: false, new System.Text.UTF8Encoding(false)))
                write(writer);
            File.Move(temp, path, overwrite: true);
        }
        finally
        {
            if (File.Exists(temp))
                File.Delete(temp);
        }
    }

    private static IEnumerable<string> Lines(string text) =>
        text.Split('\n').Select(l => l.TrimEnd('\r'));

    /// <summary>The <c>#!key value</c> lines a GO table opens with, as mzLib wrote them.</summary>
    private static Dictionary<string, string> ReadHeader(IEnumerable<string> lines)
    {
        var header = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (string line in lines)
        {
            if (!line.StartsWith("#!", StringComparison.Ordinal))
                break;
            string body = line[2..];
            int space = body.IndexOf(' ');
            header[space < 0 ? body : body[..space]] = space < 0 ? "" : body[(space + 1)..];
        }
        return header;
    }

    /// <summary>The category table, read back from the text mzLib's <see cref="GoCategoryTsv"/> wrote.</summary>
    private static object CategoryTable(GoCategoryMap map, string text)
    {
        var columns = NewColumns(GoCategoryColumns);
        int rowCount = 0;
        bool pastHeader = false;
        foreach (string line in Lines(text))
        {
            if (line.Length == 0 || line.StartsWith("#!", StringComparison.Ordinal))
                continue;
            if (!pastHeader)
            {
                pastHeader = true;
                continue;
            }
            string[] cells = line.Split('\t');
            AddRow(columns, cells[0], cells[1], cells.Length > 2 && cells[2].Length > 0 ? cells[2] : null);
            rowCount++;
        }
        return new
        {
            map_name = map.MapName,
            map_version = map.MapVersion,
            source_file_name = map.SourceFileName,
            sha256 = map.SourceSha256,
            anchor_count = map.Anchors.Count,
            row_count = rowCount,
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    private static List<string> GoCaveats(GeneOntologyGraph ontology, GoGroupAnnotator annotator, LoadedDatabase database)
    {
        var caveats = new List<string>();
        if (annotator.UnresolvedGoIds.Count > 0)
            caveats.Add(
                $"{annotator.UnresolvedGoIds.Count} GO id(s) the database cites are absent from go.obo release " +
                $"{ontology.Release ?? "(unversioned)"} and were dropped (listed in unresolved_go_ids). A protein " +
                "whose only terms were dropped reads no_go_terms. Use a go.obo at least as new as the database " +
                "to keep them.");
        if (ontology.Release is null)
            caveats.Add("go.obo has no data-version header, so go_release is null on every row; go_obo_sha256 still pins the file.");
        caveats.AddRange(database.Caveats);
        return caveats;
    }

    private static string FileSha256(string path)
    {
        using FileStream file = File.OpenRead(path);
        return Convert.ToHexString(SHA256.HashData(file)).ToLowerInvariant();
    }
}
