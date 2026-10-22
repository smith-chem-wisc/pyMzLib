using System.Globalization;
using FlashLFQ;
using Quantification.Differential;
using Readers;

namespace MzLibBridge;

/// <summary>
/// <c>diff observations</c>: mzLib's observation table (mzLib #1455), every peptide in every
/// biological sample as log2 intensity or no value with the reason, built from a stored FlashLFQ or
/// MetaMorpheus peptide table by <see cref="FlashLfqObservations.FromPeptideTable"/>, with mzLib's
/// shared-peptide normalization measured (and, if asked, applied) on it.
/// </summary>
/// <remarks>
/// <para>
/// <b>Runs come on stdin,</b> one per line: <c>file&lt;TAB&gt;sample[&lt;TAB&gt;fraction[&lt;TAB&gt;technical
/// replicate]]</c>, each an <see cref="ObservationRun"/>. A fraction or technical replicate left out
/// is not sent, so the record's own default (0) applies. Grouping comes only from these lines; mzLib
/// never parses a file name for it.
/// </para>
/// <para>
/// <b>Normalization is always measured, and applied only with <c>--normalize</c>.</b> mzLib measures
/// the shifts "whatever the setting, because the metadata always reports them"
/// (<see cref="SharedPeptideNormalization"/>), so they are on every answer. <c>--samples</c> narrows
/// the table to one stratum's samples first (<see cref="ObservationTable.ForSamples"/>), because the
/// shifts are measured over every sample of the table they are given.
/// </para>
/// <para>
/// <b>The responses table</b> (<c>--responses-out</c>) is the table's own peptide-by-sample shape,
/// <see cref="ObservationTable.Log2Intensity"/> under each sample, keyed by full sequence, which
/// mzLib keeps unique. It is not a roll-up and chooses no rows: mzLib does not own that step yet
/// (QuantProject M5b), so the bridge does not take it.
/// </para>
/// </remarks>
internal static partial class Differential
{
    /// <summary>Row columns of the <c>diff observations</c> table: one row per (peptide, sample).</summary>
    internal static readonly string[] ObservationColumns =
    {
        "full_sequence", "base_sequence", "protein_groups", "unique_protein_group", "sample", "log2_intensity", "state",
    };

    /// <summary>The wire's basis names: mzLib's own (<see cref="QuantBasisNames.ToMachineName"/>).</summary>
    private static readonly QuantBasis[] Bases = { QuantBasis.MsmsOnly, QuantBasis.MbrKept };

    /// <summary>Every <see cref="ObservationState"/>, in mzLib's declaration order (most to least informative).</summary>
    internal static readonly ObservationState[] States =
    {
        ObservationState.Quantified, ObservationState.MbrTransferred, ObservationState.AmbiguousPeak,
        ObservationState.IdentifiedNotQuantified, ObservationState.NotDetected,
    };

    /// <summary>
    /// mzLib's <see cref="ObservationState"/> as the wire writes it. mzLib gives the states no machine
    /// names of their own, so they are named here, every one explicitly: a state mzLib adds fails a
    /// test on the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string StateName(ObservationState state) => state switch
    {
        ObservationState.Quantified => "quantified",
        ObservationState.MbrTransferred => "mbr_transferred",
        ObservationState.AmbiguousPeak => "ambiguous_peak",
        ObservationState.IdentifiedNotQuantified => "identified_not_quantified",
        ObservationState.NotDetected => "not_detected",
        _ => throw new InvalidOperationException(
            $"ObservationState.{state} has no wire name. Name it in Differential.StateName."),
    };

    /// <summary>
    /// <c>diff observations --peptides FILE --basis msms_only|mbr_kept [--protein-groups FILE]
    /// [--samples FILE] [--normalize] [--min-peptides N] [--out FILE.tsv] [--responses-out FILE.tsv]</c>,
    /// with the runs on stdin.
    /// </summary>
    public static object Observations(Program.Arguments arguments)
    {
        string peptidesPath = arguments.Required("peptides");
        string basisText = arguments.Required("basis");
        QuantBasis[] named = Bases.Where(b => b.ToMachineName() == basisText).ToArray();
        if (named.Length == 0)
            throw new Program.UsageException(
                $"Option --basis must be one of {string.Join(", ", Bases.Select(b => b.ToMachineName()))}; got '{basisText}'. " +
                "mzLib has no default: whether match-between-runs values enter is a scientific choice.");
        QuantBasis basis = named[0];

        foreach (string name in new[] { "protein-groups", "samples", "out", "responses-out" })
            RequireValueIfProvided(arguments, name);
        string? groupsPath = arguments.Optional("protein-groups");
        string? samplesPath = arguments.Optional("samples");
        if (samplesPath is not null && !File.Exists(samplesPath))
            throw new Program.UsageException($"The --samples file was not found: '{samplesPath}'.");
        string? outPath = arguments.Optional("out");
        string? responsesPath = arguments.Optional("responses-out");
        Program.RequireTsvOutput(outPath);
        Program.RequireTsvOutput(responsesPath, "responses-out");
        bool normalize = arguments.Flag("normalize");
        // Absent, nothing is sent and mzLib's own floor applies (DefaultMinimumPeptides).
        int? minimumPeptides = arguments.WasProvided("min-peptides") ? arguments.OptionalInt("min-peptides", 0) : null;
        if (minimumPeptides is < 1)
            throw new Program.UsageException($"Option --min-peptides must be 1 or more; got {minimumPeptides}.");

        List<ObservationRun> runs = ReadRuns(Program.ReadStdinLines());

        var peptideFile = Reading.OpenAs<QuantifiedPeptideFile>(
            peptidesPath, "diff observations", nameof(SupportedFileType.FlashLFQQuantifiedPeptide),
            "a FlashLFQ QuantifiedPeptides.tsv, or a MetaMorpheus AllQuantifiedPeptides.tsv");
        IReadOnlyList<ProteinGroupInfo> groups = groupsPath is null
            ? Array.Empty<ProteinGroupInfo>()
            : FlashLfqObservations.ReadProteinGroups(Reading.OpenAs<ProteinGroupFromTsvFile>(
                groupsPath, "diff observations", nameof(SupportedFileType.MetaMorpheusQuantifiedProteinGroups),
                "a MetaMorpheus AllQuantifiedProteinGroups.tsv"));
        List<string>? keep = samplesPath is null ? null : FirstColumn(samplesPath, "samples");

        ObservationTable table;
        NormalizationResult normalization;
        try
        {
            table = FlashLfqObservations.FromPeptideTable(peptideFile, runs, groups, basis);
            if (keep is not null)
                table = table.ForSamples(keep);
            normalization = minimumPeptides is int floor
                ? SharedPeptideNormalization.Measure(table, floor)
                : SharedPeptideNormalization.Measure(table);
            if (normalize)
                table = SharedPeptideNormalization.Apply(table, normalization);
        }
        catch (ArgumentException refusal) when (refusal is not ArgumentNullException)
        {
            throw Refused(refusal);
        }

        var stateCounts = States.ToDictionary(StateName, _ => 0, StringComparer.Ordinal);
        var rows = new Sdrf.ColumnTable(ObservationColumns);
        object? written = null;
        int rowCount = 0;
        Reading.TsvSink? sink = outPath is null ? null : new Reading.TsvSink(outPath, ObservationColumns);
        try
        {
            foreach (Observation cell in table.Observations())
            {
                stateCounts[StateName(cell.State)]++;
                object?[] row =
                {
                    cell.Peptide.FullSequence, cell.Peptide.BaseSequence, cell.Peptide.ProteinGroups.ToList(),
                    cell.Peptide.UniqueProteinGroup, cell.SampleId, Finite(cell.Log2Intensity), StateName(cell.State),
                };
                if (sink is null)
                    rows.Add(row);
                else
                    sink.Write(row);
                rowCount++;
            }
            written = sink?.Complete();
        }
        finally
        {
            sink?.Dispose();
        }

        object? responsesWritten = null;
        if (responsesPath is not null)
        {
            WritePlain(responsesPath, new[] { "full_sequence" }.Concat(table.Samples),
                Enumerable.Range(0, table.Peptides.Count).Select(p =>
                    new[] { table.Peptides[p].FullSequence }.Concat(Enumerable.Range(0, table.Samples.Count)
                        .Select(s => table.Log2Intensity(p, s)).Select(v => double.IsNaN(v) ? "" : Number(v)))));
            responsesWritten = new
            {
                path = Path.GetFullPath(responsesPath),
                format = "tsv",
                row_count = table.Peptides.Count,
                sample_count = table.Samples.Count,
            };
        }

        return new
        {
            peptides_file = Path.GetFullPath(peptidesPath),
            protein_groups_file = groupsPath is null ? null : Path.GetFullPath(groupsPath),
            samples_file = samplesPath is null ? null : Path.GetFullPath(samplesPath),
            basis = basis.ToMachineName(),
            normalize,
            min_peptides = minimumPeptides,
            run_count = table.Runs.Count,
            runs = table.Runs.Select(r => new
            {
                file_name = r.FileName,
                sample_id = r.SampleId,
                fraction = r.Fraction,
                technical_replicate = r.TechnicalReplicate,
            }).ToList(),
            sample_count = table.Samples.Count,
            samples = table.Samples.ToList(),
            samples_without_values = table.SamplesWithoutValues.ToList(),
            peptide_count = table.Peptides.Count,
            protein_group_count = table.ProteinGroups.Count,
            protein_groups = table.ProteinGroups.Values.Select(g => new
            {
                name = g.Name,
                q_value = Finite(g.QValue),
                is_decoy = g.IsDecoy,
                is_contaminant = g.IsContaminant,
                is_entrapment = g.IsEntrapment,
                genes = g.Genes,
            }).ToList(),
            state_counts = stateCounts,
            normalization = new
            {
                setting = normalization.Setting,
                reference_set_size = normalization.ReferenceSetSize,
                minimum_reference_set_size = normalization.MinimumReferenceSetSize,
                shifts = table.Samples.ToDictionary(
                    s => s, s => Finite(normalization.PerSampleShift[s]), StringComparer.Ordinal),
                warnings = normalization.Warnings.ToList(),
            },
            warnings = table.Warnings.ToList(),
            written,
            responses_written = responsesWritten,
            row_count = rowCount,
            column_names = ObservationColumns,
            columns = rows.Columns(),
            caveats = ObservationCaveats(groupsPath is not null, normalize, normalization),
        };
    }

    private static List<ObservationRun> ReadRuns(List<string> lines)
    {
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No runs were given on stdin. Give one per line: file<TAB>sample[<TAB>fraction[<TAB>technical replicate]], " +
                "the file named as the peptide table's Intensity_ columns name it, without its extension.");
        var runs = new List<ObservationRun>(lines.Count);
        for (int i = 0; i < lines.Count; i++)
        {
            string[] f = lines[i].TrimEnd('\r').Split('\t').Select(x => x.Trim()).ToArray();
            string where = $"Line {i + 1} of stdin";
            if (f.Length is < 2 or > 4)
                throw new Program.UsageException(
                    $"{where} has {f.Length} tab-separated fields; a run is file, sample, then optionally fraction and technical replicate.");
            if (f[0].Length == 0 || f[1].Length == 0)
                throw new Program.UsageException($"{where} has an empty file or sample.");
            var run = new ObservationRun(f[0], f[1]);
            if (OptionalInteger(f, 2, where, "fraction") is int fraction)
                run = run with { Fraction = fraction };
            if (OptionalInteger(f, 3, where, "technical replicate") is int replicate)
                run = run with { TechnicalReplicate = replicate };
            runs.Add(run);
        }
        return runs;
    }

    private static int? OptionalInteger(string[] fields, int index, string where, string what)
    {
        if (fields.Length <= index || fields[index].Length == 0)
            return null;
        if (!int.TryParse(fields[index], NumberStyles.Integer, CultureInfo.InvariantCulture, out int value))
            throw new Program.UsageException($"{where}: the {what} must be an integer; got '{fields[index]}'.");
        return value;
    }

    private static List<string> ObservationCaveats(bool hasGroups, bool normalized, NormalizationResult normalization)
    {
        var caveats = new List<string>
        {
            "A missing value is null, never 0, and state says why. Fractions of one technical replicate are summed on the linear scale, then a sample's technical replicates are averaged on the log2 scale (mzLib's ObservationTable).",
            "Under msms_only a match-between-runs value enters no sum, but its cell still says mbr_transferred; under mbr_kept it enters, and the cell says mbr_transferred.",
            "The peptide table does not record the design FlashLFQ quantified with, so mzLib cannot check that these runs group the files as FlashLFQ did. FlashLFQ withholds every fraction of a sample whose strongest fraction is ambiguous; give the runs the design it quantified with, or a value it withheld can become a partial sum.",
            "FlashLFQ writes its peptide table in the culture of the machine it ran on and mzLib reads it as invariant (a decimal point), so a table written where the decimal separator is a comma does not read correctly.",
            "Every peptide is a row: decoy, contaminant and entrapment peptides, and peptides of groups at any q-value, included. mzLib does not yet choose the rows a differential test uses, or roll peptides up to proteins, so --responses-out is peptide by sample only. Filter rows with protein_groups and unique_protein_group.",
            "The normalization shifts are always measured (SharedPeptideNormalization.Measure) over every sample of the table, after --samples; they are subtracted only with --normalize. Measure one stratum at a time. A sample whose shift is null is not shifted.",
        };
        if (!hasGroups)
            caveats.Add(
                "No --protein-groups: no peptide can be told to be a decoy, contaminant or entrapment, so all of them were eligible for the normalization's reference set.");
        if (normalization.Setting == SharedPeptideNormalization.None)
            caveats.Add(
                normalized
                    ? "Normalization was asked for, but too few reference peptides were found, so no sample was shifted (setting none); normalization.warnings says how many."
                    : "Too few reference peptides were found to measure a shift (setting none); normalization.warnings says how many.");
        return caveats;
    }
}
