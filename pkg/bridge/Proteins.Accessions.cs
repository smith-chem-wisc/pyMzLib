using MzLibUtil;
using Omics.BioPolymer;

namespace MzLibBridge;

/// <summary>
/// <c>proteins parse-accession</c>: read accessions back into their parent entry, isoform or
/// version, and applied sequence variants, with mzLib's <see cref="VariantApplication.ParseAccession"/>
/// (mzLib #1382) and, through it, <see cref="ProteinAccession.Parse"/>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Projection only.</b> The grammars (UniProt's published accession grammar, RefSeq's protein
/// prefixes) and the variant-suffix rule are mzLib's. The bridge hands each stdin line to mzLib exactly
/// as read, with no trimming, case folding or splitting, because mzLib parses and never repairs: a
/// stray space or a decoy prefix must come back unrecognized rather than be quietly fixed here.
/// mzLib #1382 exists so that no downstream tool re-implements this, and the bridge is one of them.
/// </para>
/// <para>
/// <b>Namespace is lower case on the wire</b> (<c>uniprot</c>, <c>refseq</c>, <c>unrecognized</c>),
/// as <c>genes resolve</c> already writes the same <see cref="AccessionNamespace"/>: that is mzLib's
/// own text form of it (GeneResolutionTsv's namespace column), so the two verbs join on one vocabulary.
/// </para>
/// </remarks>
internal static partial class Proteins
{
    internal static readonly string[] AccessionColumns =
    {
        "input", "entry_accession", "isoform", "version", "namespace", "parent_accession",
        "applied_variants", "has_applied_variants",
    };

    /// <summary>
    /// <c>proteins parse-accession</c>, with the accessions on stdin, one per line.
    /// </summary>
    public static object ParseAccession(Program.Arguments arguments)
    {
        List<string> inputs = Program.ReadStdinLines();
        if (inputs.Count == 0)
            throw new Program.UsageException(
                "No accessions on stdin: send one per line, e.g. 'P04637-2' or 'NP_000537.3_R72P'.");

        var columns = NewColumns(AccessionColumns);
        var namespaceCounts = Enum.GetValues<AccessionNamespace>().ToDictionary(NamespaceName, _ => 0);
        int withVariants = 0;

        foreach (string input in inputs)
        {
            ProteoformAccession parsed = VariantApplication.ParseAccession(input);
            ProteinAccession entry = parsed.Entry;
            namespaceCounts[NamespaceName(entry.Namespace)]++;
            if (parsed.HasAppliedVariants)
                withVariants++;

            AddRow(columns, parsed.Verbatim, entry.EntryAccession, entry.Isoform, entry.Version,
                NamespaceName(entry.Namespace), entry.Verbatim, parsed.AppliedVariants, parsed.HasAppliedVariants);
        }

        return new
        {
            record_count = inputs.Count,
            namespace_counts = namespaceCounts,
            with_applied_variants_count = withVariants,
            column_names = columns.Keys.ToList(),
            columns,
            caveats = AccessionCaveats(),
        };
    }

    /// <summary>mzLib's text form of a namespace, as GeneResolutionTsv writes it.</summary>
    private static string NamespaceName(AccessionNamespace value) => value.ToString().ToLowerInvariant();

    private static List<string> AccessionCaveats() => new()
    {
        "Parses, never repairs (mzLib's rule): an accession outside the UniProt and RefSeq protein grammars is " +
        "unrecognized and kept verbatim, entry_accession equal to input. That includes lower case, surrounding " +
        "spaces, a FASTA header (sp|P04637|P53_HUMAN), a decoy, contaminant or entrapment prefix (DECOY_P04637), " +
        "and ProteinDbLoader's collision counter (P04637_2, which names a different entry whose accession " +
        "collided).",

        "entry_accession names the entry, not a sequence: P04637 and P04637-2 share it but are different " +
        "sequences, and the isoform UniProt displays for an entry can carry any number. It answers 'same entry', " +
        "never 'same sequence'.",

        "A protein-group accession joined with '|' (P04637|P02768) is unrecognized as a whole: mzLib does not " +
        "split it. Split protein groups first.",

        "Only RefSeq protein accessions (NP_, XP_, YP_, WP_, AP_) are recognized; a RefSeq transcript (NM_) or an " +
        "Ensembl id (ENSP...) is unrecognized.",

        "A variant suffix is read only after an entry mzLib recognizes, in the form mzLib's GetAccession writes: " +
        "'_', the original residues, a 1-based position, the variant residues (S70N; a deletion A80; a stop gain " +
        "Q5*), repeated for several (S70N_A80T). applied_variants is that suffix as written: it is not checked " +
        "against any sequence.",

        "An isoform or version number too large for a 32-bit integer (over 2147483647) reads as null while the " +
        "accession still counts as UniProt or RefSeq: mzLib keeps the entry and drops the number.",

        "Each stdin line is one accession, given to mzLib exactly as read; blank lines are skipped. One row per " +
        "line, in input order, duplicates included.",
    };
}
