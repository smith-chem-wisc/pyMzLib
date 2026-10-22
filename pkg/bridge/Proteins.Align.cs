using Omics.SequenceAlignment;

namespace MzLibBridge;

/// <summary>
/// <c>proteins align</c>: global pairwise alignment of protein sequences with mzLib's
/// <see cref="PairwiseAligner"/> (mzLib #1428), and the residue-to-residue map
/// <see cref="PairwiseAlignment.TryMapSourceToTarget"/> and
/// <see cref="PairwiseAlignment.TryMapTargetToSource"/> read off it.
/// </summary>
/// <remarks>
/// <para>
/// <b>Projection only.</b> The aligner, its scoring and its tie-breaking are mzLib's. The bridge
/// builds one <see cref="PairwiseAligner"/> per call from the options the caller gave, and from
/// nothing else: an option left out is left to mzLib's own default, read back from a default-built
/// aligner rather than restated here, and the envelope echoes the values the aligner actually used.
/// The substitution matrix is always mzLib's default (<see cref="SubstitutionMatrix.Blosum62"/>);
/// the bridge passes <c>null</c> and lets mzLib choose.
/// </para>
/// <para>
/// <b>Pairs travel on stdin, always.</b> One pair per line, source and target separated by a tab.
/// Not on argv: two titin isoforms are about 70,000 characters, past argv's 32 KB ceiling on Windows,
/// so a one-pair form on argv would fail for exactly the proteins that most need aligning. Many pairs
/// are one call (bridge PARALLELISM.md: one process, one JIT warm-up). They are aligned one after
/// another: mzLib's aligner takes no degree of parallelism, so there is no <c>--threads</c>.
/// </para>
/// <para>
/// <b>Every residue is mapped, both ways.</b> The residue map has one row per residue of each
/// sequence: every source position through <see cref="PairwiseAlignment.TryMapSourceToTarget"/>,
/// then every target position through <see cref="PairwiseAlignment.TryMapTargetToSource"/>. A residue
/// facing a gap maps to null, never to the nearest residue, which is mzLib's rule (the Try methods
/// return false for it).
/// </para>
/// </remarks>
internal static partial class Proteins
{
    /// <summary>The pairs table: one row per input pair, in input order.</summary>
    internal static readonly string[] AlignmentColumns =
    {
        "pair_index", "source_length", "target_length", "score", "identical", "substituted",
        "source_residues_facing_gap", "target_residues_facing_gap", "aligned_source", "aligned_target",
    };

    /// <summary>The residue map: one row per residue of either sequence of each pair.</summary>
    internal static readonly string[] ResidueMapColumns =
    {
        "pair_index", "sequence", "position", "residue", "mapped_position", "mapped_residue",
    };

    /// <summary>
    /// <c>proteins align [--gap-open N] [--gap-extend N] [--end-gaps free|charged] [--max-cells N]
    /// [--no-residue-map]</c>, with the pairs on stdin, one per line: <c>source TAB target</c>.
    /// </summary>
    public static object Align(Program.Arguments arguments)
    {
        PairwiseAligner aligner = AlignerFrom(arguments);
        bool withMap = !arguments.Flag("no-residue-map");

        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No sequence pairs on stdin: send one pair per line, the source and target sequences " +
                "separated by a tab, e.g. 'MKTAYIAKQRQ<TAB>MKTAYIAKISF'.");

        var pairs = new List<(string Source, string Target)>(lines.Count);
        for (int i = 0; i < lines.Count; i++)
        {
            string[] fields = lines[i].Split('\t');
            if (fields.Length != 2)
                throw new Program.UsageException(
                    $"Pair {i} has {fields.Length} tab-separated field(s); each stdin line is one pair, " +
                    "the source and target sequences separated by exactly one tab.");
            pairs.Add((fields[0], fields[1]));
        }

        var columns = NewColumns(AlignmentColumns);
        var residueMap = withMap ? NewColumns(ResidueMapColumns) : null;
        int mapRows = 0;

        for (int i = 0; i < pairs.Count; i++)
        {
            (string source, string target) = pairs[i];
            PairwiseAlignment alignment;
            try
            {
                alignment = aligner.Align(source, target);
            }
            catch (ArgumentException refused) when (refused is not ArgumentNullException)
            {
                // mzLib's own refusals: a gap character in a sequence, or a table larger than
                // MaxCells. Both are the caller's input, so they exit 2, naming the pair.
                throw new Program.UsageException(
                    $"Pair {i} was refused by mzLib's PairwiseAligner: {refused.Message}");
            }

            AddRow(columns, i, source.Length, target.Length, alignment.Score, alignment.Identical,
                alignment.Substituted, alignment.SourceResiduesFacingGap, alignment.TargetResiduesFacingGap,
                alignment.AlignedSource, alignment.AlignedTarget);

            if (residueMap is null)
                continue;
            for (int position = 1; position <= source.Length; position++)
            {
                bool paired = alignment.TryMapSourceToTarget(position, out int mapped);
                AddRow(residueMap, i, "source", position, source[position - 1].ToString(),
                    paired ? mapped : null, paired ? target[mapped - 1].ToString() : null);
            }
            for (int position = 1; position <= target.Length; position++)
            {
                bool paired = alignment.TryMapTargetToSource(position, out int mapped);
                AddRow(residueMap, i, "target", position, target[position - 1].ToString(),
                    paired ? mapped : null, paired ? source[mapped - 1].ToString() : null);
            }
            mapRows += source.Length + target.Length;
        }

        return new
        {
            aligner_id = aligner.Id,
            matrix = aligner.Matrix.Name,
            gap_open = aligner.GapOpen,
            gap_extend = aligner.GapExtend,
            end_gaps = aligner.FreeEndGaps ? "free" : "charged",
            max_cells = aligner.MaxCells,
            record_count = pairs.Count,
            column_names = columns.Keys.ToList(),
            columns,
            residue_map = residueMap is null ? null : SubTable(residueMap, mapRows),
            caveats = AlignCaveats(),
        };
    }

    /// <summary>
    /// The aligner the options describe. An option not given is mzLib's default: the value a
    /// default-built <see cref="PairwiseAligner"/> carries, so no default is restated here.
    /// </summary>
    private static PairwiseAligner AlignerFrom(Program.Arguments arguments)
    {
        var defaults = new PairwiseAligner();
        bool any = false;

        int gapOpen = defaults.GapOpen;
        if (arguments.WasProvided("gap-open"))
        {
            gapOpen = arguments.OptionalInt("gap-open", defaults.GapOpen);
            any = true;
        }

        int gapExtend = defaults.GapExtend;
        if (arguments.WasProvided("gap-extend"))
        {
            gapExtend = arguments.OptionalInt("gap-extend", defaults.GapExtend);
            any = true;
        }

        bool freeEndGaps = defaults.FreeEndGaps;
        if (arguments.WasProvided("end-gaps"))
        {
            string? raw = arguments.Optional("end-gaps");
            if (string.IsNullOrWhiteSpace(raw))
                throw new Program.UsageException("Option --end-gaps was given without a value.");
            freeEndGaps = raw.ToLowerInvariant() switch
            {
                "free" => true,
                "charged" => false,
                _ => throw new Program.UsageException(
                    $"Option --end-gaps must be 'free' or 'charged' (mzLib's PairwiseAligner.FreeEndGaps " +
                    $"true or false); got '{raw}'."),
            };
            any = true;
        }

        long maxCells = defaults.MaxCells;
        if (arguments.WasProvided("max-cells"))
        {
            maxCells = arguments.OptionalLong("max-cells", defaults.MaxCells);
            any = true;
        }

        if (!any)
            return defaults;

        try
        {
            // matrix: null is mzLib's own choice of matrix (BLOSUM62), as for the default aligner.
            return new PairwiseAligner(matrix: null, gapOpen, gapExtend, freeEndGaps, maxCells);
        }
        catch (ArgumentOutOfRangeException refused)
        {
            // mzLib's constructor owns the ranges (gap open not negative, gap extend and max cells
            // above zero, max cells at most Array.MaxLength). Its message names its own parameter, so
            // the wire option is named beside it.
            string option = refused.ParamName switch
            {
                "gapOpen" => "--gap-open",
                "gapExtend" => "--gap-extend",
                "maxCells" => "--max-cells",
                _ => "an option",
            };
            throw new Program.UsageException(
                $"{option} was refused by mzLib's PairwiseAligner: {refused.Message.Split(" (Parameter")[0]}");
        }
    }

    private static List<string> AlignCaveats() => new()
    {
        "score is a raw substitution-matrix score in BLOSUM62's half-bit units, summed over the alignment: " +
        "it grows with length and is not a bit score, an E-value or a percent identity. Compare scores only " +
        "between alignments with the same aligner_id.",

        "Global alignment: every residue of both sequences is in the alignment, paired or facing a gap. With " +
        "end gaps free (mzLib's default) an overhang at either end costs nothing, so a short sequence is " +
        "placed on the best-scoring stretch of a long one even when it has no real match there.",

        "A residue facing a gap has no counterpart: its mapped_position is null, never the nearest residue.",

        "One optimal alignment is returned. Where several score the same, mzLib's traceback prefers a residue " +
        "pair, then a gap in the target, then a gap in the source, so the same input always gives the same " +
        "alignment; a residue beside an indel may sit on the other side of the gap in another aligner's equally " +
        "good answer.",

        "mzLib refuses only the gap characters '-' and '.'. Any other character BLOSUM62 does not name " +
        "(selenocysteine U, pyrrolysine O, J, a digit, a space, a bracket) is scored as X, so pass bare " +
        "one-letter sequences. Letters are compared case-insensitively; aligned_source, aligned_target and the " +
        "residue columns keep the case given.",

        "Time and memory grow with the product of the two lengths: mzLib builds a (source length + 1) x (target " +
        "length + 1) table of one byte per cell and refuses a pair larger than max_cells. Two identical " +
        "sequences free of X and * are aligned residue for residue without the table, so max_cells does not " +
        "limit them.",

        "Each stdin line is one pair, aligned in input order; blank lines are skipped. mzLib's aligner takes no " +
        "thread count, so pairs are aligned one after another.",
    };
}
