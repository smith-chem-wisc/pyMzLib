using System.Collections.Concurrent;
using System.IO.Compression;
using System.Security.Cryptography;
using MzLibUtil;
using Omics.Modifications;
using Proteomics;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge;

/// <summary>
/// Protein databases: what a UniProt XML or FASTA says about each protein, which Ensembl genes those
/// proteins belong to, and whether a peptide identifies one of them.
/// </summary>
/// <remarks>
/// <para>
/// Three verbs, one loader. <c>proteins read</c>, <c>genes resolve</c> and
/// <c>proteins classify-peptides</c> all load their databases through <see cref="LoadDatabases"/>,
/// which calls mzLib's <see cref="ProteinDbLoader"/> exactly as a search would, minus the parts a
/// search adds on top: no decoys are generated and no sequence variants are expanded. Every row
/// therefore describes an entry as the file wrote it, once.
/// </para>
/// <para>
/// <b>The bridge adds no biology here.</b> GO terms are <see cref="Protein.GoTerms"/> (mzLib #1336),
/// Ensembl links are <see cref="Protein.EnsemblGeneReferences"/>, gene resolution is
/// <see cref="EnsemblGeneResolver"/> (#1338), and peptide sharing is
/// <see cref="PeptideUniquenessClassifier"/> (#1348). What the bridge owns is the shape: long tables
/// in place of nested objects, one path or many (bridge BULK.md), and the per-file facts that say
/// when an empty column means "this format cannot say" rather than "there is nothing".
/// </para>
/// </remarks>
internal static class Proteins
{
    /// <summary>The tables <c>proteins read</c> can return, in wire order.</summary>
    internal static readonly string[] TableNames = { "proteins", "go_terms", "ensembl_genes" };

    /// <summary>The protein table's columns, declared once so the header and the values cannot disagree.</summary>
    internal static readonly string[] ProteinColumns =
    {
        "source_index", "source_path", "accession", "name", "full_name", "organism", "ncbi_taxonomy_id",
        "primary_gene_name", "gene_names", "length", "monoisotopic_mass", "is_contaminant", "is_decoy",
        "is_entrapment", "ensembl_gene_ids",
    };

    internal static readonly string[] GoTermColumns =
    {
        "source_index", "source_path", "accession", "go_id", "aspect", "term_name", "evidence_codes", "projects",
    };

    internal static readonly string[] EnsemblGeneColumns =
    {
        "source_index", "source_path", "accession", "transcript_id", "protein_id", "gene_id",
        "versioned_gene_id", "gene_version",
    };

    /// <summary>
    /// The resolution table's columns. The names after the two source columns are mzLib's own
    /// <see cref="GeneResolutionTsv.Schema"/> headers, in its order; a test holds them equal, so a
    /// resolution read here and one written by mzLib's TSV writer use one vocabulary.
    /// </summary>
    internal static readonly string[] GeneResolutionColumns =
    {
        "source_index", "source_path", "accession", "entry_accession", "isoform", "namespace", "outcome",
        "n_genes", "gene_id", "versioned_gene_id", "gene_symbol", "gene_biotype", "off_primary_genes",
        "uniprot_gene_name", "source", "search_database_sha256", "gene_set_release", "gene_set_sha256",
        "ensembl_xref_agrees", "ensembl_xref_info_type", "ensembl_xref_sha256",
    };

    internal static readonly string[] ClassificationColumns =
    {
        "peptide", "sharing", "accession_count", "accessions", "shared_gene_keys",
    };

    /// <summary>The line that separates two lists sharing stdin, e.g. database paths and then accessions.</summary>
    internal const string StdinSectionSeparator = "--";

    // ---- proteins read ------------------------------------------------------------------------

    /// <summary>
    /// <c>proteins read (--path FILE [--contaminant] | --paths-stdin) [--tables proteins,go_terms,ensembl_genes]
    /// [--accessions-stdin] [--sequences] [--threads 1] [--on-error fail|skip]</c> — one row per protein,
    /// and on request one row per GO annotation and one per Ensembl transcript link.
    /// </summary>
    public static object Read(Program.Arguments arguments)
    {
        HashSet<string> tables = TablesFrom(arguments);
        bool withSequences = arguments.Flag("sequences");
        bool filterRequested = arguments.Flag("accessions-stdin");

        (List<DatabaseInput> inputs, List<string>? filterLines) = InputsFrom(arguments, "accessions-stdin");
        List<string>? filter = null;
        if (filterRequested)
        {
            filter = filterLines!.Select(a => a.Trim()).Where(a => a.Length > 0).Distinct(StringComparer.Ordinal).ToList();
            if (filter.Count == 0)
                throw new Program.UsageException(
                    "--accessions-stdin was given but no accessions were supplied. Omit it to read every protein.");
        }

        LoadedDatabase[] loaded = LoadDatabases(inputs, ThreadsFrom(arguments), OnErrorFrom(arguments, allowSkip: true));
        var wanted = filter is null ? null : new HashSet<string>(filter, StringComparer.Ordinal);

        var proteinColumns = tables.Contains("proteins")
            ? NewColumns(withSequences ? ProteinColumns.Append("sequence") : ProteinColumns)
            : null;
        var goColumns = tables.Contains("go_terms") ? NewColumns(GoTermColumns) : null;
        var ensemblColumns = tables.Contains("ensembl_genes") ? NewColumns(EnsemblGeneColumns) : null;
        var found = new HashSet<string>(StringComparer.Ordinal);
        int recordCount = 0;
        int goCount = 0;
        int ensemblCount = 0;
        int unweighable = 0;

        foreach (LoadedDatabase database in loaded)
        {
            if (database.Proteins is null)
                continue;

            foreach (Protein protein in database.Proteins)
            {
                if (wanted is not null && !wanted.Contains(protein.Accession))
                    continue;
                found.Add(protein.Accession);
                database.RecordCount++;
                recordCount++;

                if (proteinColumns is not null)
                {
                    double? mass = MonoisotopicMassOf(protein.BaseSequence);
                    if (mass is null)
                        unweighable++;
                    AddRow(proteinColumns, database.Input.Index, database.AbsolutePath, protein.Accession,
                        NullIfEmpty(protein.Name), NullIfEmpty(protein.FullName), NullIfEmpty(protein.Organism),
                        NullIfEmpty(protein.NcbiTaxonomyId), PrimaryGeneName(protein), GeneNames(protein),
                        protein.Length, mass, protein.IsContaminant, protein.IsDecoy, protein.IsEntrapment,
                        protein.EnsemblGeneIds.ToList());
                    if (withSequences)
                        proteinColumns["sequence"].Add(protein.BaseSequence);
                }

                if (goColumns is not null)
                {
                    foreach (GoTerm term in protein.GoTerms)
                    {
                        goCount++;
                        AddRow(goColumns, database.Input.Index, database.AbsolutePath, protein.Accession, term.Id,
                            term.Aspect.ToString(), NullIfEmpty(term.Name), Sorted(term.EvidenceCodes),
                            Sorted(term.Projects));
                    }
                }

                if (ensemblColumns is not null)
                {
                    foreach (EnsemblGeneReference link in protein.EnsemblGeneReferences)
                    {
                        ensemblCount++;
                        AddRow(ensemblColumns, database.Input.Index, database.AbsolutePath, protein.Accession,
                            NullIfEmpty(link.TranscriptId), NullIfEmpty(link.ProteinId), link.GeneId,
                            link.VersionedGeneId, link.GeneVersion);
                    }
                }
            }
        }

        var caveats = new List<string>();
        if (unweighable > 0)
            caveats.Add(
                $"{unweighable} protein(s) have no monoisotopic_mass: the sequence holds a letter with no " +
                "defined residue mass (X, B, Z, J or similar), so any mass would be a guess.");

        return new
        {
            file_count = loaded.Length,
            read_count = loaded.Count(d => d.Error is null),
            failed_count = loaded.Count(d => d.Error is not null),
            record_count = recordCount,
            tables = TableNames.Where(tables.Contains).ToList(),
            accession_filter_count = filter?.Count,
            accessions_not_found = filter?.Where(a => !found.Contains(a)).ToList(),
            caveats,
            files = loaded.Select(d => d.ToWire(extra: null)).ToList(),
            column_names = proteinColumns?.Keys.ToList(),
            columns = proteinColumns,
            go_terms = goColumns is null ? null : SubTable(goColumns, goCount),
            ensembl_genes = ensemblColumns is null ? null : SubTable(ensemblColumns, ensemblCount),
        };
    }

    // ---- genes resolve ------------------------------------------------------------------------

    /// <summary>
    /// <c>genes resolve (--path FILE [--contaminant] | --paths-stdin) (--gtf FILE | --gene-set FILE)
    /// [--xref FILE] [--threads 1] [--on-error fail|skip]</c> — every protein resolved to stable Ensembl
    /// gene ids, counted against one release's gene set.
    /// </summary>
    public static object ResolveGenes(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "gtf");
        RequireValueIfProvided(arguments, "gene-set");
        RequireValueIfProvided(arguments, "xref");
        string? gtf = arguments.Optional("gtf");
        string? geneSetTable = arguments.Optional("gene-set");
        string? xrefPath = arguments.Optional("xref");

        if ((gtf is null) == (geneSetTable is null))
            throw new Program.UsageException(
                "Give exactly one of --gtf (an Ensembl GTF, e.g. Homo_sapiens.GRCh38.116.gtf.gz) or " +
                "--gene-set (a compact table mzLib's EnsemblGeneSetWriter made from one). The gene set " +
                "is what every resolution is counted against, so it is required and never defaulted.");

        string referencePath = gtf ?? geneSetTable!;
        if (!File.Exists(referencePath))
            throw new Program.UsageException($"{(gtf is null ? "Gene set" : "GTF")} not found: '{referencePath}'.");
        if (xrefPath is not null && !File.Exists(xrefPath))
            throw new Program.UsageException($"Ensembl xref table not found: '{xrefPath}'.");

        (List<DatabaseInput> inputs, _) = InputsFrom(arguments, secondSection: null);
        int threads = ThreadsFrom(arguments);
        OnError onError = OnErrorFrom(arguments, allowSkip: true);

        EnsemblGeneSet geneSet = gtf is not null ? EnsemblGeneSet.LoadGtf(gtf) : EnsemblGeneSetReader.Load(geneSetTable!);
        EnsemblXrefTable? xrefs = xrefPath is null ? null : EnsemblXrefTable.Load(xrefPath);
        var resolver = new EnsemblGeneResolver(geneSet, xrefs);

        LoadedDatabase[] loaded = LoadDatabases(inputs, threads, onError, hashDatabase: true);

        var columns = NewColumns(GeneResolutionColumns);
        var outcomeCounts = Enum.GetValues<GeneResolutionOutcome>()
            .ToDictionary(o => GeneResolutionTsv.OutcomeName(o), _ => 0);
        int proteinCount = 0;
        int rowCount = 0;

        foreach (LoadedDatabase database in loaded)
        {
            if (database.Proteins is null)
                continue;

            foreach (Protein protein in database.Proteins)
            {
                IReadOnlyList<GeneResolution> rows = resolver.Resolve(protein, database.Sha256);
                proteinCount++;
                database.RecordCount += rows.Count;
                outcomeCounts[GeneResolutionTsv.OutcomeName(rows[0].Outcome)]++;

                foreach (GeneResolution row in rows)
                {
                    rowCount++;
                    AddRow(columns, database.Input.Index, database.AbsolutePath, row.Accession, row.EntryAccession,
                        row.Isoform, row.Namespace.ToString().ToLowerInvariant(),
                        GeneResolutionTsv.OutcomeName(row.Outcome), row.GeneCount, row.GeneId, row.VersionedGeneId,
                        row.GeneSymbol, row.GeneBiotype, row.OffPrimaryGenes, row.UniProtGeneName, row.Source,
                        row.SearchDatabaseSha256, row.GeneSetRelease, row.GeneSetSha256, row.EnsemblXrefAgrees,
                        row.EnsemblXrefInfoType, row.EnsemblXrefSha256);
                }
            }
        }

        return new
        {
            file_count = loaded.Length,
            read_count = loaded.Count(d => d.Error is null),
            failed_count = loaded.Count(d => d.Error is not null),
            protein_count = proteinCount,
            record_count = rowCount,
            gene_set = new
            {
                source_file_name = geneSet.SourceFileName,
                sha256 = geneSet.SourceSha256,
                release = geneSet.Release,
                genome_build = geneSet.GenomeBuild,
                genebuild_last_updated = geneSet.GenebuildLastUpdated,
                gene_count = geneSet.Count,
            },
            xref = xrefs is null
                ? null
                : new
                {
                    source_file_name = xrefs.SourceFileName,
                    sha256 = xrefs.SourceSha256,
                    release = xrefs.Release,
                    accession_count = xrefs.AccessionCount,
                },
            outcome_counts = outcomeCounts,
            caveats = GeneCaveats(geneSet, xrefs, loaded),
            files = loaded.Select(d => d.ToWire(extra: d.Sha256)).ToList(),
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    /// <summary>Per-call traps a caller would otherwise only find by reading mzLib's source.</summary>
    private static List<string> GeneCaveats(EnsemblGeneSet geneSet, EnsemblXrefTable? xrefs, LoadedDatabase[] loaded)
    {
        var caveats = new List<string>();
        int fromFasta = loaded.Where(d => d.FileType == "Fasta" && d.Proteins is not null)
            .Sum(d => d.Proteins!.Count(p => !p.IsContaminant));
        if (fromFasta > 0)
            caveats.Add(
                $"{fromFasta} protein(s) came from a FASTA, which carries no Ensembl links, so each is " +
                "not_in_source whatever Ensembl knows about it. Resolve against the UniProt XML of the same " +
                "proteome, or pass --xref to add Ensembl's own links as ensembl_xref rows.");
        if (geneSet.Release is null)
            caveats.Add(
                $"The gene set's file name '{geneSet.SourceFileName}' carries no Ensembl release number, so " +
                "gene_set_release is null. Keep Ensembl's own file name (Species.Assembly.Release.gtf.gz) so the " +
                "release travels with every row; the sha256 still pins the exact file.");
        if (geneSet.GenomeBuild is null)
            caveats.Add(
                "The gene set has no '#!genome-build' header, so the assembly it covers is not recorded.");
        if (xrefs is null)
            caveats.Add(
                "No --xref was given, so ensembl_xref_agrees is null on every row: Ensembl's agreement is " +
                "unknown, not false.");
        return caveats;
    }

    // ---- proteins classify-peptides -----------------------------------------------------------

    /// <summary>
    /// <c>proteins classify-peptides (--path FILE [--contaminant] | --paths-stdin) [--threads 1]</c>, with the
    /// peptides on stdin (after the paths and a <c>--</c> line when <c>--paths-stdin</c> is given) — one row
    /// per peptide saying whether it is unique to one sequence, shared within a gene, or shared across genes.
    /// </summary>
    public static object ClassifyPeptides(Program.Arguments arguments)
    {
        // A database that failed to load would silently shrink the search space and turn "shared" into
        // "unique" — the one direction of error this verb exists to prevent. So there is no skip mode.
        if (arguments.WasProvided("on-error") && arguments.Optional("on-error") != "fail")
            throw new Program.UsageException(
                "proteins classify-peptides has no --on-error skip. Every database is part of one search " +
                "space: dropping one that failed to load would report peptides it contains as unique.");

        (List<DatabaseInput> inputs, List<string>? peptideLines) = InputsFrom(arguments, "peptides");
        List<string> peptides = (peptideLines ?? new List<string>()).Select(p => p.Trim()).Where(p => p.Length > 0).ToList();
        if (peptides.Count == 0)
            throw new Program.UsageException(
                "No peptides were supplied on stdin. Give one unmodified base sequence per line " +
                (arguments.Flag("paths-stdin") ? "after the database paths and a line holding only '--'." : "."));

        LoadedDatabase[] loaded = LoadDatabases(inputs, ThreadsFrom(arguments), OnError.Fail);
        List<Protein> proteins = loaded.SelectMany(d => d.Proteins!).ToList();
        foreach (LoadedDatabase database in loaded)
            database.RecordCount = database.Proteins!.Count(p => !p.IsDecoy);

        IReadOnlyList<PeptideUniqueness> results;
        try
        {
            results = PeptideUniquenessClassifier.Classify(peptides, proteins);
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            // mzLib's own validation, reclassified: a lower-case or modified sequence is the caller's
            // input, not a fault, and must exit 2 like every other malformed argument.
            throw new Program.UsageException(
                invalid.Message.Split(" (Parameter")[0] +
                " Strip modifications and pass the base sequence in upper case, e.g. 'PEPTIDEK'.");
        }

        var columns = NewColumns(ClassificationColumns);
        var sharingCounts = Enum.GetValues<PeptideSharing>().ToDictionary(s => s.ToString(), _ => 0);
        foreach (PeptideUniqueness result in results)
        {
            sharingCounts[result.Sharing.ToString()]++;
            AddRow(columns, result.Peptide, result.Sharing.ToString(), result.Accessions.Count,
                result.Accessions.ToList(), result.SharedGeneKeys.ToList());
        }

        return new
        {
            file_count = loaded.Length,
            peptide_count = results.Count,
            target_protein_count = proteins.Count(p => !p.IsDecoy),
            decoy_proteins_ignored = proteins.Count(p => p.IsDecoy),
            i_and_l_equivalent = true,
            sharing_counts = sharingCounts,
            files = loaded.Select(d => d.ToWire(extra: null)).ToList(),
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    // ---- the shared loader --------------------------------------------------------------------

    /// <summary>One database to load: its position in the input, its path, and whether it is a contaminant.</summary>
    internal sealed record DatabaseInput(int Index, string Path, bool Contaminant);

    internal enum OnError { Fail, Skip }

    /// <summary>A database after loading: its proteins, or the failure that took their place.</summary>
    internal sealed class LoadedDatabase
    {
        public required DatabaseInput Input { get; init; }
        public required string AbsolutePath { get; init; }
        public string? FileType { get; set; }
        public string? Reader { get; set; }
        public List<Protein>? Proteins { get; set; }
        public List<string> Caveats { get; } = new();
        public List<string> AbsentFields { get; } = new();
        public string? Sha256 { get; set; }
        public object? Error { get; set; }

        /// <summary>Rows this database contributed to the verb's main table.</summary>
        public int RecordCount { get; set; }

        /// <summary>The per-file block (bridge BULK.md §3).</summary>
        public object ToWire(string? extra) => new
        {
            source_index = Input.Index,
            path = AbsolutePath,
            file_type = FileType,
            reader = Reader,
            contaminant = Input.Contaminant,
            protein_count = Proteins?.Count ?? 0,
            decoy_count = Proteins?.Count(p => p.IsDecoy) ?? 0,
            record_count = RecordCount,
            search_database_sha256 = extra,
            caveats = Caveats,
            absent_fields = AbsentFields,
            error = Error,
        };
    }

    /// <summary>
    /// Loads every database, <paramref name="threads"/> at a time, and returns them in input order.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Output order never depends on <paramref name="threads"/>: each load writes its own slot, and a
    /// failure under <see cref="OnError.Fail"/> is rethrown for the LOWEST failing index, so the error a
    /// caller sees is the same at any thread count.
    /// </para>
    /// <para>
    /// mzLib's own <c>maxThreads</c> is pinned to 1. With no decoys generated it has nothing to do, and a
    /// value of -1 per load would multiply against this loop's degree — the unowned oversubscription
    /// bridge PARALLELISM.md warns about.
    /// </para>
    /// </remarks>
    internal static LoadedDatabase[] LoadDatabases(
        IReadOnlyList<DatabaseInput> inputs, int threads, OnError onError, bool hashDatabase = false)
    {
        var loaded = new LoadedDatabase[inputs.Count];
        var failures = new ConcurrentDictionary<int, Exception>();
        var options = new ParallelOptions { MaxDegreeOfParallelism = threads == -1 ? Environment.ProcessorCount : threads };

        Parallel.For(0, inputs.Count, options, i =>
        {
            var database = new LoadedDatabase { Input = inputs[i], AbsolutePath = Path.GetFullPath(inputs[i].Path) };
            loaded[i] = database;
            try
            {
                Load(database);
                if (hashDatabase)
                    database.Sha256 = DecompressedSha256(database.AbsolutePath);
            }
            catch (Exception exception)
            {
                database.Proteins = null;
                failures[i] = exception;
            }
        });

        if (failures.IsEmpty)
            return loaded;

        if (onError == OnError.Fail)
        {
            Exception first = failures[failures.Keys.Min()];
            System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(first).Throw();
        }

        foreach ((int index, Exception exception) in failures)
        {
            Exception cause = Program.Unwrap(exception);
            loaded[index].Error = new
            {
                kind = cause is Program.UsageException ? "usage" : "correctness",
                type = cause is Program.UsageException ? "usage" : Program.ClassifyError(exception),
                message = cause.Message,
            };
        }
        return loaded;
    }

    /// <summary>Loads one database with mzLib's loader, choosing XML or FASTA from its name.</summary>
    private static void Load(LoadedDatabase database)
    {
        string path = database.AbsolutePath;
        DatabaseFormat format = FormatOf(path);
        if (!File.Exists(path))
            throw new Program.UsageException($"Protein database not found: '{database.Input.Path}'.");

        bool contaminant = database.Input.Contaminant;

        if (format == DatabaseFormat.UniProtXml)
        {
            database.FileType = "UniProtXml";
            database.Reader = "ProteinDbLoader.LoadProteinXML";
            // maxHeterozygousVariants: 0 applies no sequence variant, so each entry is one row under its
            // own accession. The default (4) would add "P12345_S70N"-style variant proteins, which a
            // database read is not asking for.
            database.Proteins = ProteinDbLoader.LoadProteinXML(
                path, generateTargets: true, DecoyType.None, Array.Empty<Modification>(), contaminant,
                modTypesToExclude: null, out _, maxThreads: 1, maxHeterozygousVariants: 0);
        }
        else
        {
            database.FileType = "Fasta";
            database.Reader = "ProteinDbLoader.LoadProteinFasta";
            database.Proteins = ProteinDbLoader.LoadProteinFasta(
                path, generateTargets: true, DecoyType.None, contaminant, out List<string> errors, maxThreads: 1);
            database.Caveats.AddRange(errors);
            // A FASTA header has no dbReference lines, so it can say nothing about GO or Ensembl. An empty
            // go_terms table for a FASTA is "this format cannot say", never "no annotations exist".
            database.AbsentFields.AddRange(new[] { "go_terms", "ensembl_genes", "ensembl_gene_ids" });
            database.Caveats.Add(
                "FASTA carries no GO terms or Ensembl links, so this file contributes no go_terms or " +
                "ensembl_genes rows and every ensembl_gene_ids is empty. Read the UniProt XML for those.");
        }

        int variantApplied = database.Proteins.Count(p => p.AppliedSequenceVariations.Count > 0);
        if (variantApplied > 0)
            database.Caveats.Add(
                $"{variantApplied} protein(s) carry sequence variants mzLib applied from the genotype this " +
                "database records, so their accession is the entry's plus a variant suffix (e.g. " +
                "'P38936_C117Y') and their length and mass are the variant's, not the reference entry's.");

        if (database.Proteins.Count == 0)
            database.Caveats.Add("No proteins were read from this database.");
    }

    internal enum DatabaseFormat { UniProtXml, Fasta }

    /// <summary>The format a database path names, by extension, with an optional trailing <c>.gz</c>.</summary>
    internal static DatabaseFormat FormatOf(string path)
    {
        string name = Path.GetFileName(path).ToLowerInvariant();
        if (name.EndsWith(".gz", StringComparison.Ordinal))
            name = name[..^3];

        if (name.EndsWith(".xml", StringComparison.Ordinal))
            return DatabaseFormat.UniProtXml;
        if (new[] { ".fasta", ".fa", ".faa", ".fas" }.Any(e => name.EndsWith(e, StringComparison.Ordinal)))
            return DatabaseFormat.Fasta;

        throw new Program.UsageException(
            $"'{Path.GetFileName(path)}' is not a protein database this verb reads. Give a UniProt XML " +
            "(.xml or .xml.gz) or a FASTA (.fasta, .fa, .faa or .fas, optionally .gz).");
    }

    /// <summary>
    /// Lower-case hex sha256 of the database's DECOMPRESSED bytes — what
    /// <see cref="EnsemblGeneResolver.Resolve"/> documents as the key of a resolution, so the same
    /// database gzipped or not keys the same rows.
    /// </summary>
    internal static string DecompressedSha256(string path)
    {
        using FileStream file = File.OpenRead(path);
        using Stream content = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
            ? new GZipStream(file, CompressionMode.Decompress)
            : file;
        return Convert.ToHexString(SHA256.HashData(content)).ToLowerInvariant();
    }

    // ---- argument shaping ---------------------------------------------------------------------

    /// <summary>
    /// The databases to load, and the second stdin list when this verb takes one.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <c>--path FILE</c> names one database; <c>--contaminant</c> marks it as one. <c>--paths-stdin</c>
    /// reads one per line, each optionally followed by a tab and <c>contaminant</c> or <c>target</c>.
    /// Giving both is a usage error, as are a repeated path and a blank one (bridge BULK.md §1).
    /// </para>
    /// <para>
    /// When the verb also reads a list of its own from stdin (accessions, peptides) and the paths are on
    /// stdin too, the two lists are separated by a line holding only <c>--</c>. One stream, two sections,
    /// and no temporary file: argv cannot carry thousands of accessions, and a binding must not have to
    /// invent a side channel.
    /// </para>
    /// </remarks>
    internal static (List<DatabaseInput> Inputs, List<string>? Second) InputsFrom(
        Program.Arguments arguments, string? secondSection)
    {
        bool pathsFromStdin = arguments.Flag("paths-stdin");
        RequireValueIfProvided(arguments, "path");
        string? single = arguments.Optional("path");

        if (pathsFromStdin && single is not null)
            throw new Program.UsageException("Give --path or --paths-stdin, not both.");
        if (!pathsFromStdin && single is null)
            throw new Program.UsageException("Missing required option --path (or --paths-stdin for many databases).");
        if (pathsFromStdin && arguments.Flag("contaminant"))
            throw new Program.UsageException(
                "--contaminant applies to --path. With --paths-stdin, mark a database by following its path " +
                "with a tab and 'contaminant'.");

        // "peptides" is always on stdin for classify-peptides; "accessions-stdin" only when asked for.
        bool secondOnStdin = secondSection == "peptides" || (secondSection is not null && arguments.Flag(secondSection));
        List<string> stdin = pathsFromStdin || secondOnStdin ? Program.ReadStdinLines() : new List<string>();

        List<string> pathLines;
        List<string>? second = null;
        if (pathsFromStdin && secondOnStdin)
        {
            int separator = stdin.FindIndex(l => l.Trim() == StdinSectionSeparator);
            if (separator < 0)
                throw new Program.UsageException(
                    $"With --paths-stdin, stdin holds the database paths, then a line holding only " +
                    $"'{StdinSectionSeparator}', then the {(secondSection == "peptides" ? "peptides" : "accessions")}. " +
                    "No such line was found.");
            pathLines = stdin.Take(separator).ToList();
            second = stdin.Skip(separator + 1).ToList();
        }
        else if (pathsFromStdin)
        {
            pathLines = stdin;
        }
        else
        {
            pathLines = new List<string>();
            second = secondOnStdin ? stdin : null;
        }

        var inputs = new List<DatabaseInput>();
        if (single is not null)
        {
            inputs.Add(new DatabaseInput(0, single, arguments.Flag("contaminant")));
        }
        else
        {
            var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (string line in pathLines)
            {
                string[] fields = line.Split('\t');
                string path = fields[0].Trim();
                if (path.Length == 0)
                    throw new Program.UsageException("A stdin line has an empty database path.");
                string role = fields.Length > 1 ? fields[1].Trim().ToLowerInvariant() : "target";
                if (role is not ("target" or "contaminant" or ""))
                    throw new Program.UsageException(
                        $"Database '{path}' is marked '{fields[1].Trim()}'; the only roles are 'target' and 'contaminant'.");
                if (!seen.Add(Path.GetFullPath(path)))
                    throw new Program.UsageException($"Database '{path}' is listed twice.");
                inputs.Add(new DatabaseInput(inputs.Count, path, role == "contaminant"));
            }
            if (inputs.Count == 0)
                throw new Program.UsageException("--paths-stdin was given but no database paths were supplied.");
        }

        return (inputs, second);
    }

    /// <summary><c>--threads N</c>: databases loaded at once. Default 1; -1 means every core.</summary>
    internal static int ThreadsFrom(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "threads");
        int threads = arguments.OptionalInt("threads", 1);
        if (threads == 0 || threads < -1)
            throw new Program.UsageException($"Option --threads must be 1 or more, or -1 for every core; got {threads}.");
        return threads;
    }

    internal static OnError OnErrorFrom(Program.Arguments arguments, bool allowSkip)
    {
        RequireValueIfProvided(arguments, "on-error");
        return arguments.Optional("on-error") switch
        {
            null or "fail" => OnError.Fail,
            "skip" when allowSkip => OnError.Skip,
            string other => throw new Program.UsageException($"Option --on-error must be 'fail' or 'skip'; got '{other}'."),
        };
    }

    /// <summary><c>--tables</c>: a comma list of <see cref="TableNames"/>. Default: <c>proteins</c> alone.</summary>
    internal static HashSet<string> TablesFrom(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "tables");
        string raw = arguments.Optional("tables") ?? "proteins";
        var tables = raw.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .ToHashSet(StringComparer.Ordinal);
        List<string> unknown = tables.Where(t => !TableNames.Contains(t)).ToList();
        if (unknown.Count > 0 || tables.Count == 0)
            throw new Program.UsageException(
                $"Option --tables takes a comma list of {string.Join(", ", TableNames)}; got '{raw}'.");
        return tables;
    }

    private static void RequireValueIfProvided(Program.Arguments arguments, string name)
    {
        if (arguments.WasProvided(name) && string.IsNullOrWhiteSpace(arguments.Optional(name)))
            throw new Program.UsageException($"Option --{name} was given without a value.");
    }

    // ---- value shaping ------------------------------------------------------------------------

    /// <summary>
    /// The monoisotopic mass of the unmodified sequence, in daltons (mzLib's residue masses plus one
    /// water), or null when a letter has no defined mass.
    /// </summary>
    internal static double? MonoisotopicMassOf(string sequence)
    {
        try
        {
            return new Proteomics.AminoAcidPolymer.Peptide(sequence).MonoisotopicMass;
        }
        catch (MzLibException)
        {
            return null;
        }
    }

    private static string? PrimaryGeneName(Protein protein) =>
        NullIfEmpty(protein.GeneNames?.FirstOrDefault(n => n.Item1 == "primary")?.Item2);

    private static List<string> GeneNames(Protein protein) =>
        (protein.GeneNames ?? new List<Tuple<string, string>>())
            .Select(n => n.Item2).Where(n => !string.IsNullOrEmpty(n)).ToList();

    private static string? NullIfEmpty(string? value) => string.IsNullOrEmpty(value) ? null : value;

    private static List<string> Sorted(IEnumerable<string> values) => values.OrderBy(v => v, StringComparer.Ordinal).ToList();

    private static Dictionary<string, List<object?>> NewColumns(IEnumerable<string> names) =>
        names.ToDictionary(n => n, _ => new List<object?>());

    /// <summary>Appends one row, in the columns' declared order.</summary>
    private static void AddRow(Dictionary<string, List<object?>> columns, params object?[] values)
    {
        int i = 0;
        foreach (List<object?> column in columns.Values)
        {
            if (i >= values.Length)
                break;
            column.Add(values[i++]);
        }
    }

    private static object SubTable(Dictionary<string, List<object?>> columns, int rowCount) => new
    {
        row_count = rowCount,
        column_names = columns.Keys.ToList(),
        columns,
    };
}
