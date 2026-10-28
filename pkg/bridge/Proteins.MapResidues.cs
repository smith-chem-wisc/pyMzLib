using System.Globalization;
using UsefulProteomicsDatabases.Ensembl;

namespace MzLibBridge;

/// <summary>
/// <c>genes map-residues</c>: which residue of one protein sits in the same column of Ensembl Compara's
/// gene-tree alignment as a residue of another, answered by mzLib's
/// <see cref="ComparaGeneTreeAlignment.MapResidue"/> (mzLib #1465).
/// </summary>
/// <remarks>
/// <para>
/// <b>Projection only.</b> mzLib reads the gene-tree content (<see cref="ComparaGeneTreeContent.Load"/>),
/// keeps the canonical proteins it names from the alignment (<see cref="ComparaGeneTreeAlignment.Load"/>),
/// and maps each query. The bridge parses its own stdin, checks that the files exist, and shapes the
/// answer. The outcome is mzLib's <see cref="ComparaColumnOutcome"/> name, verbatim: mzLib has no
/// wire-name helper for it, as it has <see cref="OrthologySnapshot.StatusName"/> for orthology statuses,
/// and the bridge does not invent one.
/// </para>
/// <para>
/// <b>The gene sets only restrict.</b> The species lines are passed to
/// <see cref="ComparaGeneTreeContent.Load"/> as <c>restrictTo</c>, and nothing else. With none, mzLib's
/// default applies: every gene in the gene-tree file is kept, and so every canonical protein in the
/// alignment (millions, for the vertebrate collection), which a caveat says.
/// </para>
/// <para>
/// <b>A position outside the source protein stops the call.</b> <see cref="ComparaGeneTreeAlignment.MapResidue"/>
/// throws <see cref="ArgumentOutOfRangeException"/> for it, where every other case is an outcome. The
/// bridge does not check lengths itself (that would re-derive mzLib's rule); it reports mzLib's refusal
/// as a usage error naming the query.
/// </para>
/// <para>
/// <b>No <c>--threads</c>.</b> Neither load takes a degree of parallelism, and the queries are a loop over
/// a held alignment (bridge PARALLELISM.md).
/// </para>
/// </remarks>
internal static partial class Proteins
{
    /// <summary>The mapping table's columns, in wire order.</summary>
    internal static readonly string[] ResidueMappingColumns =
    {
        "source_protein_id", "position", "target_protein_id", "outcome", "column", "target_position", "target_residue",
    };

    /// <summary>The <c>--proteins</c> table's columns, in wire order: <see cref="ComparaAlignedProtein"/>.</summary>
    internal static readonly string[] AlignedProteinColumns = { "protein_id", "alignment_index", "tree_id", "sequence" };

    /// <summary>One stdin query: the source protein, a 1-based position in it, and the target protein.</summary>
    private sealed record ResidueQuery(string Source, int Position, string Target);

    /// <summary>
    /// <c>genes map-residues --alignment FILE --gene-trees FILE [--proteins]</c>, with any species lines,
    /// a <c>--</c> line and the queries on stdin.
    /// </summary>
    public static object MapResidues(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "alignment");
        RequireValueIfProvided(arguments, "gene-trees");
        // A flag. A value after it ("--proteins true") would land among the named options, and the flag
        // would read as absent.
        if (arguments.Optional("proteins") is { } stray)
            throw new Program.UsageException($"--proteins takes no value (got '{stray}'); give --proteins alone.");

        string alignmentPath = arguments.Required("alignment");
        string geneTreesPath = arguments.Required("gene-trees");
        bool withProteins = arguments.Flag("proteins");

        (List<GeneSetInput> geneSetInputs, List<ResidueQuery> queries) = MapResidueInputsFrom(Program.ReadStdinLines());

        foreach (GeneSetInput input in geneSetInputs)
            RequireFile(input.Path, input.Kind == GtfKind ? $"GTF of {input.Species}" : $"Gene-set table of {input.Species}");
        RequireFile(geneTreesPath, "Gene-trees file");
        RequireFile(alignmentPath, "Compara gene-tree alignment");

        var geneSets = new SortedDictionary<string, EnsemblGeneSet>(StringComparer.Ordinal);
        foreach (GeneSetInput input in geneSetInputs)
            geneSets[input.Species] = input.Kind == GtfKind
                ? EnsemblGeneSet.LoadGtf(input.Path)
                : EnsemblGeneSetReader.Load(input.Path);

        ComparaGeneTreeContent geneTrees = ComparaGeneTreeContent.Load(geneTreesPath, geneSets.Count > 0 ? geneSets.Values : null);
        ComparaGeneTreeAlignment alignment = ComparaGeneTreeAlignment.Load(alignmentPath, geneTrees);

        var columns = NewColumns(ResidueMappingColumns);
        var counts = Enum.GetValues<ComparaColumnOutcome>().ToDictionary(o => o.ToString(), _ => 0);
        for (int i = 0; i < queries.Count; i++)
        {
            ResidueQuery query = queries[i];
            ComparaColumnMapping mapping;
            try
            {
                mapping = alignment.MapResidue(query.Source, query.Position, query.Target);
            }
            catch (ArgumentOutOfRangeException refused)
            {
                // mzLib's one refusal, for a position outside the source protein. The caller's input, so
                // a usage error, in mzLib's words; and it stops the call, as mzLib's throw does.
                throw new Program.UsageException(
                    $"Query {i + 1} ({query.Source}, position {query.Position}, {query.Target}): " +
                    $"{refused.Message.Split(" (Parameter")[0]} mzLib's MapResidue refuses a position outside the " +
                    "source protein, so no query of this call was answered.");
            }

            string outcome = mapping.Outcome.ToString();
            counts[outcome]++;
            AddRow(columns, query.Source, query.Position, query.Target, outcome, mapping.Column,
                mapping.TargetPosition, mapping.TargetResidue?.ToString());
        }

        return new
        {
            alignment = new
            {
                source_file_name = alignment.SourceFileName,
                sha256 = alignment.SourceSha256,
                release = alignment.Release,
                collection = alignment.Collection,
                alignment_count = alignment.AlignmentCount,
                protein_count = alignment.Count,
                restricted_to_gene_tree_content = alignment.RestrictedToGeneTreeContent,
            },
            gene_trees = new
            {
                source_file_name = geneTrees.SourceFileName,
                sha256 = geneTrees.SourceSha256,
                release = geneTrees.Release,
                collection = geneTrees.Collection,
                member_count = geneTrees.Count,
                restricted_to_gene_sets = geneTrees.RestrictedToGeneSets,
            },
            gene_sets = geneSets.ToDictionary(p => p.Key, p => GeneSetWire(p.Value)),
            row_count = queries.Count,
            outcome_counts = counts,
            column_names = columns.Keys.ToList(),
            columns,
            proteins_included = withProteins,
            proteins = withProteins ? AlignedProteinsWire(alignment) : null,
            caveats = MapResidueCaveats(alignment, geneTrees, counts),
        };
    }

    /// <summary>
    /// The two stdin sections: species lines (<c>species TAB gtf|gene-set TAB path</c>, none or more), a line
    /// holding only <c>--</c>, then one query per line (<c>source TAB position TAB target</c>).
    /// </summary>
    private static (List<GeneSetInput> GeneSets, List<ResidueQuery> Queries) MapResidueInputsFrom(List<string> stdin)
    {
        int separator = stdin.FindIndex(l => l.Trim() == StdinSectionSeparator);
        if (separator < 0)
            throw new Program.UsageException(
                "stdin holds the species whose proteins to keep, one line each ('homo_sapiens<TAB>gtf<TAB>" +
                "Homo_sapiens.GRCh38.116.gtf.gz', or gene-set for a compact gene table; none keeps every species), " +
                $"then a line holding only '{StdinSectionSeparator}', then one query per line " +
                "('ENSP00000269305<TAB>15<TAB>ENSMUSP00000104298'). No such line was found.");

        List<GeneSetInput> geneSets = SpeciesLinesFrom(stdin.Take(separator));

        var queries = new List<ResidueQuery>();
        foreach (string line in stdin.Skip(separator + 1).Where(l => l.Trim().Length > 0))
        {
            string[] fields = line.Split('\t').Select(f => f.Trim()).ToArray();
            if (fields.Length != 3 || fields.Any(f => f.Length == 0))
                throw new Program.UsageException(
                    $"Query '{line}' is not three tab-separated fields: the source protein id, a 1-based position " +
                    "in it, and the target protein id, e.g. 'ENSP00000269305<TAB>15<TAB>ENSMUSP00000104298'.");
            if (!int.TryParse(fields[1], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int position))
                throw new Program.UsageException(
                    $"Query '{line}': the position '{fields[1]}' is not an integer. It is a 1-based residue number.");
            queries.Add(new ResidueQuery(fields[0], position, fields[2]));
        }
        if (queries.Count == 0)
            throw new Program.UsageException(
                $"No queries after the '{StdinSectionSeparator}' line. Give one per line: the source protein id, a " +
                "1-based position in it, and the target protein id, tab-separated.");
        return (geneSets, queries);
    }

    /// <summary><see cref="ComparaGeneTreeAlignment.Proteins"/>, in mzLib's order (ordinal by protein id).</summary>
    private static object AlignedProteinsWire(ComparaGeneTreeAlignment alignment)
    {
        var columns = NewColumns(AlignedProteinColumns);
        foreach (ComparaAlignedProtein protein in alignment.Proteins)
            AddRow(columns, protein.ProteinId, protein.AlignmentIndex, protein.TreeId, protein.Sequence);
        return new
        {
            row_count = alignment.Count,
            column_names = columns.Keys.ToList(),
            columns,
        };
    }

    /// <summary>This call's traps, each one a fact of these inputs that the rows alone would hide.</summary>
    private static List<string> MapResidueCaveats(ComparaGeneTreeAlignment alignment, ComparaGeneTreeContent geneTrees,
        Dictionary<string, int> counts)
    {
        var caveats = new List<string>();

        if (geneTrees.RestrictedToGeneSets is null)
            caveats.Add(
                $"No gene sets were given, so mzLib kept every gene in {geneTrees.SourceFileName} ({geneTrees.Count}) and " +
                $"every canonical protein of theirs in the alignment ({alignment.Count}). Give the species' gene sets on " +
                "stdin to keep only theirs: the answers for their proteins are the same, and the memory far smaller.");

        if (alignment.Release is null || geneTrees.Release is null)
            caveats.Add(
                $"{(alignment.Release is null ? alignment.SourceFileName : geneTrees.SourceFileName)} carries no release in " +
                "its file name, so it cannot be told whether the alignment and the gene trees are from one release. mzLib " +
                "does not compare them. Keep Ensembl's file names (Compara.<release>.protein_default.aa.fasta.gz, " +
                "<collection>.GeneTree_content.default.e<release>.txt.gz).");
        else if (alignment.Release != geneTrees.Release)
            caveats.Add(
                $"The alignment is named for release {alignment.Release} and the gene trees for release {geneTrees.Release}. " +
                "mzLib does not hold one to the other, and a gene's canonical protein and tree can change between " +
                "releases, so a protein may be missed or placed in the wrong tree. Use the two files of one release.");

        int missing = counts[nameof(ComparaColumnOutcome.SourceNotInAlignment)] + counts[nameof(ComparaColumnOutcome.TargetNotInAlignment)];
        if (missing > 0)
            caveats.Add(
                $"{missing} quer{(missing == 1 ? "y names" : "ies name")} a protein the alignment did not keep. Compara " +
                "aligns one protein per gene, its canonical one (the gene-tree row flagged Y), and mzLib keeps only the " +
                "canonical proteins of the genes it kept: another isoform's id, or a gene outside the gene sets, is not there. " +
                "--proteins lists the proteins kept.");

        return caveats;
    }
}
