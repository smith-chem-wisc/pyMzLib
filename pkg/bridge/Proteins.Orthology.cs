using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge;

/// <summary>
/// <c>genes orthology</c>: an Ensembl Compara orthology snapshot for a chosen set of species, built by
/// mzLib's <see cref="OrthologySnapshot.Build"/> (mzLib #1381), optionally written as Parquet by
/// <see cref="OrthologySnapshotWriter.WriteAsync"/>, and optionally read back gene by gene with
/// <see cref="OrthologySnapshot.PairStatus"/>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Projection only.</b> Every rule here is mzLib's: which inputs fit together (every species'
/// dump, one release), what a relationship row holds, which status a gene gets, the Parquet layout and
/// the snapshot id. The bridge parses its own stdin, checks that the files exist, and maps mzLib's
/// refusals onto the wire's error kinds: an <see cref="ArgumentException"/> from <c>Build</c> is the
/// caller's inputs not fitting together, a usage error; an <see cref="InvalidDataException"/> is data
/// that breaks a reader's or the snapshot's rules, a correctness failure under its own type.
/// </para>
/// <para>
/// <b>Two checks are made before any input is read</b>, because the inputs are large (a human GTF is
/// 141 MB compressed) and a refusal after a minute of loading is a poor answer to a typo. Both apply
/// mzLib's own rule, not a bridge one: <see cref="OrthologySnapshotWriter.WriteAsync"/> refuses a
/// directory that is not empty, and <see cref="OrthologySnapshot.PairStatus"/> refuses a species that is
/// not one of the snapshot's, or the same species twice.
/// </para>
/// <para>
/// <b>No <c>--threads</c>.</b> None of these mzLib calls takes a degree of parallelism, so the inputs
/// are read one after another (bridge PARALLELISM.md: parallelism goes on the wire only where mzLib
/// takes it).
/// </para>
/// </remarks>
internal static partial class Proteins
{
    /// <summary>The pair-status table's columns, in wire order.</summary>
    internal static readonly string[] PairStatusColumns = { "gene_id", "gene_name", "gene_biotype", "status" };

    /// <summary>The two gene-set kinds a species line may name, as genes resolve spells its options.</summary>
    private const string GtfKind = "gtf";
    private const string GeneSetKind = "gene-set";

    /// <summary>One stdin species line: the production name, which reader, and the file.</summary>
    private sealed record GeneSetInput(string Species, string Kind, string Path);

    /// <summary>
    /// <c>genes orthology --release R --gene-trees FILE [--restrict-gene-trees] [--out DIR]
    /// [--pair-status FROM,TO]</c>, with the species lines, a <c>--</c> line and the dump paths on stdin.
    /// </summary>
    public static async Task<object> OrthologyAsync(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "release");
        RequireValueIfProvided(arguments, "gene-trees");
        RequireValueIfProvided(arguments, "out");
        RequireValueIfProvided(arguments, "pair-status");
        // A flag. A value after it ("--restrict-gene-trees true") would otherwise land among the named
        // options, the flag would read as absent, and the caller would get the unrestricted read.
        if (arguments.Optional("restrict-gene-trees") is { } stray)
            throw new Program.UsageException(
                $"--restrict-gene-trees takes no value (got '{stray}'); give --restrict-gene-trees alone.");

        string release = arguments.Required("release");
        string geneTreesPath = arguments.Required("gene-trees");
        bool restrict = arguments.Flag("restrict-gene-trees");
        string? outDirectory = arguments.Optional("out") is { } outText ? Path.GetFullPath(outText) : null;

        (List<GeneSetInput> geneSetInputs, List<string> dumpPaths) = OrthologyInputsFrom(Program.ReadStdinLines());
        var speciesGiven = geneSetInputs.Select(g => g.Species).ToHashSet(StringComparer.Ordinal);
        (string From, string To)? pairStatus = PairStatusFrom(arguments, speciesGiven);

        foreach (GeneSetInput input in geneSetInputs)
            RequireFile(input.Path, input.Kind == GtfKind ? $"GTF of {input.Species}" : $"Gene-set table of {input.Species}");
        RequireFile(geneTreesPath, "Gene-trees file");
        foreach (string dump in dumpPaths)
            RequireFile(dump, "Homology dump");
        if (outDirectory is not null)
            RequireWritableSnapshotDirectory(outDirectory);

        var geneSets = new Dictionary<string, EnsemblGeneSet>(StringComparer.Ordinal);
        foreach (GeneSetInput input in geneSetInputs)
            geneSets[input.Species] = input.Kind == GtfKind
                ? EnsemblGeneSet.LoadGtf(input.Path)
                : EnsemblGeneSetReader.Load(input.Path);

        ComparaGeneTreeContent geneTrees = ComparaGeneTreeContent.Load(geneTreesPath, restrict ? geneSets.Values : null);
        List<ComparaHomologyDump> dumps = dumpPaths.Select(p => ComparaHomologyDump.Load(p, geneSets.Keys)).ToList();

        OrthologySnapshot snapshot;
        try
        {
            snapshot = OrthologySnapshot.Build(release, geneSets, geneTrees, dumps);
        }
        catch (ArgumentException refused) when (refused is not ArgumentNullException)
        {
            // mzLib's "these inputs do not fit together": another release, a dump missing, repeated or
            // of a species not given. The caller's inputs, so a usage error, with mzLib's own words.
            throw new Program.UsageException(refused.Message.Split(" (Parameter")[0]);
        }

        SnapshotManifest? manifest = outDirectory is null
            ? null
            : await OrthologySnapshotWriter.WriteAsync(snapshot, outDirectory).ConfigureAwait(false);

        return new
        {
            release = snapshot.Release,
            species = snapshot.Species,
            format = OrthologySnapshotWriter.FormatName,
            format_version = OrthologySnapshotWriter.FormatVersion,
            snapshot_id = OrthologySnapshotWriter.SnapshotId(snapshot),
            identical_duplicates_dropped = snapshot.IdenticalDuplicatesDropped,
            gene_sets = snapshot.Species.ToDictionary(s => s, s => GeneSetWire(snapshot.GeneSets[s])),
            gene_trees = new
            {
                source_file_name = geneTrees.SourceFileName,
                sha256 = geneTrees.SourceSha256,
                release = geneTrees.Release,
                collection = geneTrees.Collection,
                member_count = geneTrees.Count,
                restricted_to_gene_sets = geneTrees.RestrictedToGeneSets,
            },
            dumps = snapshot.Dumps.ToDictionary(d => d.Genome, d => (object)new
            {
                source_file_name = d.SourceFileName,
                sha256 = d.SourceSha256,
                release = d.Release,
                rows_read = d.RowsRead,
                rows_kept = d.Rows.Count,
            }),
            pairs = snapshot.Pairs.Select(p => PairWire(snapshot, p.A, p.B)).ToList(),
            written = manifest is null
                ? null
                : new
                {
                    directory = outDirectory,
                    files = manifest.Files.Select(f => new { path = f.Path, rows = f.Rows, bytes = f.Bytes, sha256 = f.Sha256 }).ToList(),
                },
            pair_status = pairStatus is { } direction ? PairStatusWire(snapshot, direction.From, direction.To) : null,
            caveats = OrthologyCaveats(snapshot, geneTrees),
        };
    }

    /// <summary>
    /// The two stdin sections: species lines (<c>species TAB gtf|gene-set TAB path</c>), a line holding
    /// only <c>--</c>, then one dump path per line.
    /// </summary>
    private static (List<GeneSetInput> GeneSets, List<string> Dumps) OrthologyInputsFrom(List<string> stdin)
    {
        int separator = stdin.FindIndex(l => l.Trim() == StdinSectionSeparator);
        if (separator < 0)
            throw new Program.UsageException(
                "stdin holds one line per species ('homo_sapiens<TAB>gtf<TAB>Homo_sapiens.GRCh38.116.gtf.gz', or " +
                $"gene-set for a compact gene table), then a line holding only '{StdinSectionSeparator}', then one " +
                "Compara homology dump per line. No such line was found.");

        List<GeneSetInput> geneSets = SpeciesLinesFrom(stdin.Take(separator));
        if (geneSets.Count == 0)
            throw new Program.UsageException(
                "No species lines before the '--' line: give one per species, e.g. " +
                "'homo_sapiens<TAB>gtf<TAB>Homo_sapiens.GRCh38.116.gtf.gz'.");

        List<string> dumps = stdin.Skip(separator + 1).Select(l => l.Trim()).Where(l => l.Length > 0).ToList();
        if (dumps.Count == 0)
            throw new Program.UsageException(
                "No homology dumps after the '--' line. Give every species' dump " +
                "(homologies/<species>/Compara.<release>.protein_default.homologies.tsv.gz): Compara splits each " +
                "pair's relationships between the pair's two dumps, so mzLib requires them all.");
        return (geneSets, dumps);
    }

    /// <summary>
    /// Species lines (<c>species TAB gtf|gene-set TAB path</c>), shared by <c>genes orthology</c> and
    /// <c>genes map-residues</c>. A line that is not three fields, a kind other than gtf or gene-set, or a
    /// species given twice is a usage error. None is not an error here: each verb says whether it needs one.
    /// </summary>
    private static List<GeneSetInput> SpeciesLinesFrom(IEnumerable<string> lines)
    {
        var geneSets = new List<GeneSetInput>();
        var seen = new HashSet<string>(StringComparer.Ordinal);
        foreach (string line in lines)
        {
            string[] fields = line.Split('\t');
            if (fields.Length != 3 || fields.Any(f => f.Trim().Length == 0))
                throw new Program.UsageException(
                    $"Species line '{line}' is not three tab-separated fields: the species' Ensembl production " +
                    "name, gtf or gene-set, and the path, e.g. 'homo_sapiens<TAB>gtf<TAB>Homo_sapiens.GRCh38.116.gtf.gz'.");
            string species = fields[0].Trim();
            string kind = fields[1].Trim();
            if (kind is not (GtfKind or GeneSetKind))
                throw new Program.UsageException(
                    $"Species '{species}' names its gene set as '{kind}'; the kinds are '{GtfKind}' (an Ensembl GTF) and " +
                    $"'{GeneSetKind}' (a compact table mzLib's EnsemblGeneSetWriter made from one).");
            if (!seen.Add(species))
                throw new Program.UsageException($"Species '{species}' is listed twice; give one gene set per species.");
            geneSets.Add(new GeneSetInput(species, kind, fields[2].Trim()));
        }
        return geneSets;
    }

    /// <summary><c>--pair-status FROM,TO</c>: two different species given on stdin, or null when absent.</summary>
    private static (string From, string To)? PairStatusFrom(Program.Arguments arguments, HashSet<string> species)
    {
        string? raw = arguments.Optional("pair-status");
        if (raw is null)
            return null;
        string[] names = raw.Split(',', StringSplitOptions.TrimEntries);
        if (names.Length != 2 || names.Any(n => n.Length == 0))
            throw new Program.UsageException(
                $"Option --pair-status takes two species separated by a comma, e.g. 'homo_sapiens,mus_musculus'; got '{raw}'.");
        // OrthologySnapshot.PairStatus's own rule (RequireDistinct), applied before the inputs are read.
        foreach (string name in names.Where(n => !species.Contains(n)))
            throw new Program.UsageException(
                $"--pair-status names '{name}', which is not one of the species on stdin " +
                $"({string.Join(", ", species.OrderBy(s => s, StringComparer.Ordinal))}).");
        if (names[0] == names[1])
            throw new Program.UsageException(
                "Orthology is between two different species: --pair-status names one species twice.");
        return (names[0], names[1]);
    }

    private static void RequireFile(string path, string what)
    {
        if (!File.Exists(path))
            throw new Program.UsageException($"{what} not found: '{path}'.");
    }

    /// <summary>
    /// OrthologySnapshotWriter.WriteAsync's refusal (a directory that exists and is not empty), applied
    /// before the inputs are read; and a path that is a file, which the writer would fail on mid-write.
    /// </summary>
    private static void RequireWritableSnapshotDirectory(string directory)
    {
        if (File.Exists(directory))
            throw new Program.UsageException($"--out '{directory}' is a file. Give a directory for the snapshot.");
        if (Directory.Exists(directory) && Directory.EnumerateFileSystemEntries(directory).Any())
            throw new Program.UsageException(
                $"--out '{directory}' is not empty. mzLib writes a snapshot once and never updates one in place " +
                "(OrthologySnapshotWriter.WriteAsync); give a new directory, e.g. one named for the release.");
    }

    /// <summary>The same fields as genes resolve's gene_set, so one binding type reads both.</summary>
    private static object GeneSetWire(EnsemblGeneSet set) => new
    {
        source_file_name = set.SourceFileName,
        sha256 = set.SourceSha256,
        release = set.Release,
        genome_build = set.GenomeBuild,
        genebuild_last_updated = set.GenebuildLastUpdated,
        gene_count = set.Count,
    };

    private static object PairWire(OrthologySnapshot snapshot, string a, string b)
    {
        IReadOnlyList<ComparaHomology> rows = snapshot.Pair(a, b);
        return new
        {
            species_a = a,
            species_b = b,
            row_count = rows.Count,
            ortholog_count = rows.Count(r => r.Class == HomologyClass.Ortholog),
            paralog_count = rows.Count(r => r.Class == HomologyClass.Paralog),
            homoeolog_count = rows.Count(r => r.Class == HomologyClass.Homoeolog),
        };
    }

    private static object PairStatusWire(OrthologySnapshot snapshot, string from, string to)
    {
        List<(EnsemblGene Gene, OrthologyStatus Status)> rows;
        try
        {
            rows = snapshot.PairStatus(from, to).ToList();
        }
        catch (ArgumentException refused) when (refused is not ArgumentNullException)
        {
            // Checked before the inputs were read; kept so a change to mzLib's rule still exits 2.
            throw new Program.UsageException(refused.Message.Split(" (Parameter")[0]);
        }

        var columns = NewColumns(PairStatusColumns);
        var counts = Enum.GetValues<OrthologyStatus>().ToDictionary(OrthologySnapshot.StatusName, _ => 0);
        foreach ((EnsemblGene gene, OrthologyStatus status) in rows)
        {
            string name = OrthologySnapshot.StatusName(status);
            counts[name]++;
            AddRow(columns, gene.GeneId, gene.Symbol, gene.Biotype, name);
        }

        return new
        {
            from_species = from,
            to_species = to,
            row_count = rows.Count,
            status_counts = counts,
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    /// <summary>This call's traps, each one a fact of these inputs that the result alone would hide.</summary>
    private static List<string> OrthologyCaveats(OrthologySnapshot snapshot, ComparaGeneTreeContent geneTrees)
    {
        var caveats = new List<string>();

        var notChecked = snapshot.Dumps.Where(d => d.Release is null).Select(d => $"{d.SourceFileName} (the {d.Genome} dump)").ToList();
        if (geneTrees.Release is null)
            notChecked.Insert(0, geneTrees.SourceFileName);
        if (notChecked.Count > 0)
            caveats.Add(
                $"{string.Join(", ", notChecked)} carr{(notChecked.Count == 1 ? "ies" : "y")} no release in the file name, " +
                $"so mzLib could not hold {(notChecked.Count == 1 ? "it" : "them")} to --release {snapshot.Release}. Keep Ensembl's " +
                "file names (Compara.<release>.protein_default.homologies.tsv.gz, " +
                "<collection>.GeneTree_content.default.e<release>.txt.gz) so the release is checked.");

        foreach (IGrouping<string, ComparaHomologyDump> shared in snapshot.Dumps
                     .GroupBy(d => d.SourceFileName, StringComparer.Ordinal).Where(g => g.Count() > 1))
            caveats.Add(
                $"The dumps of {string.Join(", ", shared.Select(d => d.Genome))} share the file name {shared.Key}, " +
                "and the pairs tables' source_dump column holds only the file name, so it cannot say which species' " +
                "dump a row came from. manifest.json's inputs list each dump with its species and sha256.");

        foreach (string species in snapshot.Species)
        {
            string? setRelease = snapshot.GeneSets[species].Release;
            if (setRelease is not null && setRelease != snapshot.Release)
                caveats.Add(
                    $"The gene set of {species} ({snapshot.GeneSets[species].SourceFileName}) is named for release " +
                    $"{setRelease}, not {snapshot.Release}. mzLib does not compare a GTF's release: Ensembl publishes some species " +
                    "with their own annotation release. It held the gene set to the dumps by content instead: every gene " +
                    "a homology names is in it.");
        }

        if (snapshot.IdenticalDuplicatesDropped > 0)
            caveats.Add(
                $"{snapshot.IdenticalDuplicatesDropped} homology row(s) were in both of a pair's dumps, identical in every " +
                "value, and are counted once. Release 116's dumps partition, with none.");

        int inGeneSets = snapshot.Species.Sum(s => snapshot.GeneSets[s].GeneIds.Count(g => geneTrees.TryGetMember(g, out _)));
        if (geneTrees.RestrictedToGeneSets is null && geneTrees.Count > inGeneSets)
            caveats.Add(
                $"The gene trees were read unrestricted: {geneTrees.Count} genes were held in memory, {inGeneSets} of them " +
                "in the given gene sets. --restrict-gene-trees holds only those and gives the same snapshot.");

        return caveats;
    }
}
