using Omics.BioPolymerGroup;
using Readers;

namespace MzLibBridge;

/// <summary>
/// The two quantification tables mzLib 1.0.592 reads (#1347) — MetaMorpheus protein groups and
/// FlashLFQ peptides — as <b>long</b> tables: one row per record per sample.
/// </summary>
/// <remarks>
/// <para>
/// <c>read-records</c> already reads both files, and drops exactly the part anyone opens them for.
/// Each record keeps its per-sample values in a dictionary keyed by the sample's column label
/// (<see cref="ProteinGroupFromTsv.SampleGroups"/>, <see cref="QuantifiedPeptideFromTsv.Samples"/>),
/// and a dictionary has no column shape, so <c>read-records</c> names it in <c>excluded_fields</c>
/// and moves on.
/// </para>
/// <para>
/// <b>Why long, and not one column per sample.</b> A wide table would need column names made up
/// from sample labels (<c>intensity_QE-002106_GM1_a-calib</c>) — a schema mzLib never published,
/// different in every experiment, and one a multi-file read could not give a single header. Long
/// format has a fixed column set: the record's identity, then <c>sample_label</c> and that sample's
/// values. It is the shape pandas, polars, R's tidyverse and every plotting library want anyway; a
/// wide matrix is one <c>pivot</c> away, and the guide shows it.
/// </para>
/// <para>
/// PTM site occupancy is a third level of nesting — sites within the group's members within each
/// sample — so it gets its own table, <c>read-occupancy</c>, rather than a column of encoded text.
/// </para>
/// </remarks>
internal static partial class Reading
{
    // ---------------------------------------------------------------------------------------
    // Sample columns, read from the header
    // ---------------------------------------------------------------------------------------

    // mzLib matches the per-sample columns by these prefixes (ProteinGroupFromTsvFile and
    // QuantifiedPeptideFile, both private). Mirrored ONLY to report which of them a file has at all
    // — the values themselves come from mzLib's dictionaries. The order matters exactly as it does
    // there: the first prefix a header starts with wins, and "IntensityOccupancy_" must be tried
    // before "Intensity_".
    private static readonly string[] ProteinGroupSamplePrefixes =
        ["IntensityOccupancy_", "CountOccupancy_", "SpectralCount_", "Intensity_"];

    private static readonly string[] PeptideSamplePrefixes =
        ["Intensity_", "Detection Type_", "RetentionTime (min)_"];

    /// <summary>The per-sample columns of a header: which prefixes it has, and its labels in order.</summary>
    private sealed record SampleHeader(HashSet<string> Prefixes, List<string> Labels)
    {
        public static SampleHeader Of(string path, string[] prefixes)
        {
            var present = new HashSet<string>(StringComparer.Ordinal);
            var labels = new List<string>();
            foreach (string cell in HeaderCells(path))
            {
                string? prefix = prefixes.FirstOrDefault(p => cell.StartsWith(p, StringComparison.Ordinal));
                if (prefix is null)
                    continue;
                present.Add(prefix);
                string label = cell[prefix.Length..];
                if (!labels.Contains(label))
                    labels.Add(label);
            }

            return new SampleHeader(present, labels);
        }
    }

    /// <summary>
    /// Opens a path as the given reader type, or says what it is instead.
    /// </summary>
    private static TFile OpenAs<TFile>(string path, string verb, string fileType, string fileName)
        where TFile : class, IResultFile
    {
        IResultFile resultFile = OpenAny(path);
        if (resultFile is TFile typed)
            return typed;

        throw new Program.UsageException(
            $"'{resultFile.FileType}' files cannot be read by {verb}, which reads {fileName} " +
            $"(file type {fileType}). Every file type can be read with read-records, which returns " +
            "that format's own fields.");
    }

    // ---------------------------------------------------------------------------------------
    // read-protein-groups
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// <c>readers read-protein-groups --path FILE [--limit N] [--offset N] [--out FILE]</c>, or the
    /// <c>--paths-stdin</c> form — a MetaMorpheus <c>AllQuantifiedProteinGroups.tsv</c> as one row
    /// per protein group per sample group.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each row carries the group's identity and the fields a caller filters on (its q-value and
    /// decoy/contaminant/target label), then one sample group's spectral count and intensity. The
    /// group's other fields — coverage, masses, member counts — are in <c>read-records</c>, joined on
    /// <c>protein_group_name</c>. <c>--offset</c> and <c>--limit</c> count GROUPS, and each group
    /// yields one row per sample group.
    /// </para>
    /// <para>
    /// The file is written unfiltered, and so is this table: decoys, contaminants and groups above 1%
    /// FDR are all rows. That is a caveat on every read, because counting rows is the obvious thing
    /// to do and gives the wrong answer.
    /// </para>
    /// </remarks>
    public static object ReadProteinGroups(Program.Arguments arguments) =>
        RunTable(arguments, ProteinGroupsVerb);

    private static TableVerb ProteinGroupsVerb => new(
        "read-protein-groups", ReadProteinGroupsFile, ProteinGroupColumns.Select(c => c.Name).ToList(),
        ["sample_labels"], ReportsRows: true);

    /// <summary>One protein group in one sample group; <c>Sample</c> is null when the file has no sample columns.</summary>
    private sealed record GroupSample(ProteinGroupFromTsv Group, SampleGroupMeasurement? Sample);

    private static IReadOnlyList<Column<GroupSample>> ProteinGroupColumns { get; } =
    [
        // ProteinGroupFromTsv's own names. protein_group_name is the members' accessions,
        // '|'-joined in accession order — the join key back to read-records.
        new("protein_group_name", r => r.Group.ProteinGroupName),
        new("gene", r => r.Group.Gene),
        new("organism", r => r.Group.Organism),
        new("decoy_contaminant_target", r => r.Group.DecoyContaminantTarget),
        new("is_decoy", r => r.Group.IsDecoy),
        new("is_contaminant", r => r.Group.IsContaminant),
        new("is_entrapment", r => r.Group.IsEntrapment),
        new("q_value", r => r.Group.QValue),
        // SampleGroupMeasurement's own names; the label is the column suffix, verbatim.
        new("sample_label", r => r.Sample?.Label),
        new("spectral_count", r => r.Sample?.SpectralCount),
        new("intensity", r => r.Sample?.Intensity),
    ];

    /// <summary>The group-level columns mzLib reads from an optional file column, by property.</summary>
    private static readonly (string Column, string Property)[] ProteinGroupOptionalColumns =
    [
        ("gene", nameof(ProteinGroupFromTsv.Gene)),
        ("organism", nameof(ProteinGroupFromTsv.Organism)),
    ];

    /// <summary>One file through <c>read-protein-groups</c>.</summary>
    private static FileTable ReadProteinGroupsFile(string path, Window window)
    {
        var file = OpenAs<ProteinGroupFromTsvFile>(
            path, "read-protein-groups", nameof(SupportedFileType.MetaMorpheusQuantifiedProteinGroups),
            "a MetaMorpheus AllQuantifiedProteinGroups.tsv");
        List<ProteinGroupFromTsv> all = file.Results;
        IReadOnlyList<ProteinGroupFromTsv> selected = window.Apply(all, out bool truncated);

        SampleHeader header = SampleHeader.Of(path, ProteinGroupSamplePrefixes);
        var rows = new List<GroupSample>();
        foreach (ProteinGroupFromTsv group in selected)
        {
            if (group.SampleGroups.Count == 0)
                rows.Add(new GroupSample(group, null));
            else
                rows.AddRange(group.SampleGroups.Values.Select(sample => new GroupSample(group, sample)));
        }

        HashSet<string> absentProperties = AbsentProperties(typeof(ProteinGroupFromTsv), path);
        var absent = ProteinGroupOptionalColumns
            .Where(c => absentProperties.Contains(c.Property)).Select(c => c.Column).ToList();
        if (header.Labels.Count == 0)
            absent.Add("sample_label");
        if (!header.Prefixes.Contains("SpectralCount_"))
            absent.Add("spectral_count");
        if (!header.Prefixes.Contains("Intensity_"))
            absent.Add("intensity");

        var columns = ProteinGroupColumns;
        List<string> names = columns.Select(c => c.Name).ToList();

        return new FileTable
        {
            Block = Block(
                path,
                fileType: file.FileType.ToString(),
                reader: file.GetType().Name,
                recordCount: all.Count,
                rowsNotRead: null,
                retentionTimeUnit: null,
                caveats: ProteinGroupCaveats,
                columnNames: names,
                absent: absent,
                failed: FailedColumns(columns, rows),
                excluded:
                [
                    new
                    {
                        field = "count_occupancy",
                        type = nameof(ModificationOccupancyCell),
                        reason = "modified sites per sample group, nested within the group's members; one row per site in read-occupancy",
                        verb = "readers read-occupancy",
                    },
                    new
                    {
                        field = "intensity_occupancy",
                        type = nameof(ModificationOccupancyCell),
                        reason = "modified sites per sample group, nested within the group's members; one row per site in read-occupancy",
                        verb = "readers read-occupancy",
                    },
                ],
                ("sample_labels", header.Labels)),
            ColumnNames = names,
            RowCount = rows.Count,
            RecordsReturned = selected.Count,
            Row = i => Cells(columns, rows[i]),
            Offset = window.Offset,
            Truncated = truncated,
        };
    }

    private static List<string> ProteinGroupCaveats =>
    [
        "The table is UNFILTERED, as MetaMorpheus writes it: decoys, contaminants and groups above 1% " +
        "FDR are all rows (ProteinGroupFromTsv.cs:15). Filter on q_value and " +
        "decoy_contaminant_target before counting or comparing groups.",
        "sample_label is the column label verbatim, e.g. QE-002106_GM1_a-calib. Condition, replicate " +
        "and channel cannot be recovered from it (ProteinGroupFromTsv.cs:103); map labels to your " +
        "design yourself, from an SDRF or your own sample sheet.",
        "intensity null means the cell was blank: the group was not quantified in that sample group, " +
        "which is not zero (ProteinGroupFromTsv.cs:114). An isobaric file has one counting label and " +
        "one intensity label per channel, so a row can carry only spectral_count or only intensity " +
        "(ProteinGroupFromTsv.cs:104).",
        "No member is a leading protein: protein_group_name lists the members sorted by accession, so " +
        "the first accession is only the one that sorts first (ProteinGroupFromTsv.cs:17).",
        "offset and limit count GROUPS; each group yields one row per sample group.",
    ];

    // ---------------------------------------------------------------------------------------
    // read-quantified-peptides
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// <c>readers read-quantified-peptides --path FILE [--limit N] [--offset N] [--out FILE]</c>, or
    /// the <c>--paths-stdin</c> form — a FlashLFQ <c>QuantifiedPeptides.tsv</c> (or MetaMorpheus's
    /// <c>AllQuantifiedPeptides.tsv</c>) as one row per peptide per sample.
    /// </summary>
    /// <remarks>
    /// <b>A zero intensity is not a measurement.</b> FlashLFQ writes a literal <c>0</c> for a peptide
    /// it did not quantify in a sample, and mzLib keeps what was written. <c>detection_type</c> is
    /// what tells "not detected" from a measured value, which is why it is on every row.
    /// </remarks>
    public static object ReadQuantifiedPeptides(Program.Arguments arguments) =>
        RunTable(arguments, QuantifiedPeptidesVerb);

    private static TableVerb QuantifiedPeptidesVerb => new(
        "read-quantified-peptides", ReadQuantifiedPeptidesFile, PeptideColumns.Select(c => c.Name).ToList(),
        ["sample_labels"], ReportsRows: true);

    /// <summary>One peptide in one sample; <c>Sample</c> is null when the file has no sample columns.</summary>
    private sealed record PeptideSample(QuantifiedPeptideFromTsv Peptide, QuantifiedPeptideSample? Sample);

    private static IReadOnlyList<Column<PeptideSample>> PeptideColumns { get; } =
    [
        // QuantifiedPeptideFromTsv's own names. sequence is the full (modified) sequence.
        new("sequence", r => r.Peptide.Sequence),
        new("base_sequence", r => r.Peptide.BaseSequence),
        new("peak_order", r => r.Peptide.PeakOrder),
        new("protein_groups", r => r.Peptide.ProteinGroups),
        new("gene_names", r => r.Peptide.GeneNames),
        new("organism", r => r.Peptide.Organism),
        // QuantifiedPeptideSample's own names.
        new("sample_label", r => r.Sample?.Label),
        new("intensity", r => r.Sample?.Intensity),
        new("detection_type", r => r.Sample?.DetectionType),
        new("retention_time", r => r.Sample?.RetentionTime),
    ];

    private static readonly (string Column, string Property)[] PeptideOptionalColumns =
    [
        ("peak_order", nameof(QuantifiedPeptideFromTsv.PeakOrder)),
        ("protein_groups", nameof(QuantifiedPeptideFromTsv.ProteinGroups)),
        ("gene_names", nameof(QuantifiedPeptideFromTsv.GeneNames)),
        ("organism", nameof(QuantifiedPeptideFromTsv.Organism)),
    ];

    /// <summary>One file through <c>read-quantified-peptides</c>.</summary>
    private static FileTable ReadQuantifiedPeptidesFile(string path, Window window)
    {
        var file = OpenAs<QuantifiedPeptideFile>(
            path, "read-quantified-peptides", nameof(SupportedFileType.FlashLFQQuantifiedPeptide),
            "a FlashLFQ QuantifiedPeptides.tsv or a MetaMorpheus AllQuantifiedPeptides.tsv");
        List<QuantifiedPeptideFromTsv> all = file.Results;
        IReadOnlyList<QuantifiedPeptideFromTsv> selected = window.Apply(all, out bool truncated);

        SampleHeader header = SampleHeader.Of(path, PeptideSamplePrefixes);
        var rows = new List<PeptideSample>();
        foreach (QuantifiedPeptideFromTsv peptide in selected)
        {
            if (peptide.Samples.Count == 0)
                rows.Add(new PeptideSample(peptide, null));
            else
                rows.AddRange(peptide.Samples.Values.Select(sample => new PeptideSample(peptide, sample)));
        }

        HashSet<string> absentProperties = AbsentProperties(typeof(QuantifiedPeptideFromTsv), path);
        var absent = PeptideOptionalColumns
            .Where(c => absentProperties.Contains(c.Property)).Select(c => c.Column).ToList();
        if (header.Labels.Count == 0)
            absent.Add("sample_label");
        if (!header.Prefixes.Contains("Intensity_"))
            absent.Add("intensity");
        if (!header.Prefixes.Contains("Detection Type_"))
            absent.Add("detection_type");
        // IsoTracker writes a per-sample retention time; plain FlashLFQ does not.
        if (!header.Prefixes.Contains("RetentionTime (min)_"))
            absent.Add("retention_time");

        var columns = PeptideColumns;
        List<string> names = columns.Select(c => c.Name).ToList();

        return new FileTable
        {
            Block = Block(
                path,
                fileType: file.FileType.ToString(),
                reader: file.GetType().Name,
                recordCount: all.Count,
                rowsNotRead: null,
                // The header says "(min)"; mzLib reads it as written.
                retentionTimeUnit: "minutes",
                caveats: PeptideCaveats,
                columnNames: names,
                absent: absent,
                failed: FailedColumns(columns, rows),
                excluded: [],
                ("sample_labels", header.Labels)),
            ColumnNames = names,
            RowCount = rows.Count,
            RecordsReturned = selected.Count,
            Row = i => Cells(columns, rows[i]),
            Offset = window.Offset,
            Truncated = truncated,
        };
    }

    private static List<string> PeptideCaveats =>
    [
        "intensity 0 is NOT a measured zero: FlashLFQ writes a literal 0 for a peptide it did not " +
        "quantify in a sample, and mzLib keeps what was written (QuantifiedPeptideFromTsv.cs:14). " +
        "detection_type (MSMS, MBR, NotDetected, ...) is what tells them apart; filter on it before " +
        "taking a mean or a log. null means the cell was blank.",
        "sample_label is the column label verbatim; condition and replicate cannot be recovered from " +
        "it. Map labels to your design yourself.",
        "retention_time is filled only by IsoTracker output (QuantifiedPeptideFromTsv.cs:61); plain " +
        "FlashLFQ tables have no per-sample retention time, and absent_fields says so.",
        "offset and limit count PEPTIDES; each peptide yields one row per sample.",
    ];

    // ---------------------------------------------------------------------------------------
    // read-occupancy
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// <c>readers read-occupancy --path FILE [--limit N] [--offset N] [--out FILE]</c>, or the
    /// <c>--paths-stdin</c> form — the PTM site occupancy of a MetaMorpheus
    /// <c>AllQuantifiedProteinGroups.tsv</c>, one row per group, sample group, basis and site.
    /// </summary>
    /// <remarks>
    /// <para>
    /// MetaMorpheus writes two occupancy cells per group per sample group, one computed from PSM
    /// COUNTS and one from INTENSITIES. mzLib parses each into sites grouped by entity
    /// (<see cref="ModificationOccupancyCell"/>, #1347). Every site of every cell is a row here, with
    /// <c>basis</c> saying which cell it came from; the two are not paired, because they need not name
    /// the same sites.
    /// </para>
    /// <para>
    /// A malformed cell is not guessed at: mzLib refuses it with a <see cref="FormatException"/>, the
    /// cell gives no rows, and the field is named in <c>failed_fields</c>. A cell the writer cut short
    /// keeps its complete sites, each marked <c>cell_is_truncated</c>; <c>truncated_cell_count</c>
    /// counts every cut cell, including those with no complete site left to show.
    /// </para>
    /// </remarks>
    public static object ReadOccupancy(Program.Arguments arguments) =>
        RunTable(arguments, OccupancyVerb);

    private static TableVerb OccupancyVerb => new(
        "read-occupancy", ReadOccupancyFile, OccupancyColumns.Select(c => c.Name).ToList(),
        ["sample_labels", "truncated_cell_count"], ReportsRows: true);

    /// <summary>One site of one occupancy cell.</summary>
    private sealed record SiteRow(
        ProteinGroupFromTsv Group, string SampleLabel, string Basis, int EntityIndex, OccupancySite Site, bool CellIsTruncated);

    private static IReadOnlyList<Column<SiteRow>> OccupancyColumns { get; } =
    [
        new("protein_group_name", r => r.Group.ProteinGroupName),
        new("sample_label", r => r.SampleLabel),
        // "count" (CountOccupancy_) or "intensity" (IntensityOccupancy_).
        new("basis", r => r.Basis),
        new("entity_index", r => r.EntityIndex),
        // OccupancySite's own names; ModificationIdWithMotif crosses as "modification".
        new("position", r => r.Site.Position),
        new("is_n_terminus", r => r.Site.IsNTerminus),
        new("modification", r => r.Site.ModificationIdWithMotif),
        new("fraction", r => r.Site.Fraction),
        new("numerator", r => r.Site.Numerator),
        new("denominator", r => r.Site.Denominator),
        new("cell_is_truncated", r => r.CellIsTruncated),
    ];

    /// <summary>One file through <c>read-occupancy</c>.</summary>
    private static FileTable ReadOccupancyFile(string path, Window window)
    {
        var file = OpenAs<ProteinGroupFromTsvFile>(
            path, "read-occupancy", nameof(SupportedFileType.MetaMorpheusQuantifiedProteinGroups),
            "a MetaMorpheus AllQuantifiedProteinGroups.tsv");
        List<ProteinGroupFromTsv> all = file.Results;
        IReadOnlyList<ProteinGroupFromTsv> selected = window.Apply(all, out bool truncated);

        var rows = new List<SiteRow>();
        var failures = new SortedSet<string>(StringComparer.Ordinal);
        int truncatedCells = 0;

        foreach (ProteinGroupFromTsv group in selected)
        {
            foreach (SampleGroupMeasurement sample in group.SampleGroups.Values)
            {
                AddSites(group, sample.Label, "count", "count_occupancy", sample.CountOccupancyText, () => sample.CountOccupancy);
                AddSites(group, sample.Label, "intensity", "intensity_occupancy", sample.IntensityOccupancyText, () => sample.IntensityOccupancy);
            }
        }

        void AddSites(ProteinGroupFromTsv group, string label, string basis, string field, string? text,
            Func<ModificationOccupancyCell> parse)
        {
            if (text is null)
                return;

            ModificationOccupancyCell cell;
            try
            {
                cell = parse();
            }
            catch (FormatException exception)
            {
                failures.Add($"{field}: {exception.GetType().Name}");
                return;
            }

            if (cell.IsTruncated)
                truncatedCells++;

            for (int entity = 0; entity < cell.Entities.Count; entity++)
            {
                foreach (OccupancySite site in cell.Entities[entity])
                    rows.Add(new SiteRow(group, label, basis, entity, site, cell.IsTruncated));
            }
        }

        SampleHeader header = SampleHeader.Of(path, ProteinGroupSamplePrefixes);
        var absent = new List<string>();
        if (!header.Prefixes.Contains("CountOccupancy_") && !header.Prefixes.Contains("IntensityOccupancy_"))
            absent.AddRange(OccupancyColumns.Select(c => c.Name).Skip(2));

        var columns = OccupancyColumns;
        List<string> names = columns.Select(c => c.Name).ToList();

        return new FileTable
        {
            Block = Block(
                path,
                fileType: file.FileType.ToString(),
                reader: file.GetType().Name,
                recordCount: all.Count,
                rowsNotRead: null,
                retentionTimeUnit: null,
                caveats: OccupancyCaveats,
                columnNames: names,
                absent: absent,
                failed: [.. failures],
                excluded: [],
                ("sample_labels", header.Labels),
                ("truncated_cell_count", truncatedCells)),
            ColumnNames = names,
            RowCount = rows.Count,
            RecordsReturned = selected.Count,
            Row = i => Cells(columns, rows[i]),
            Offset = window.Offset,
            Truncated = truncated,
        };
    }

    private static List<string> OccupancyCaveats =>
    [
        "Trust numerator/denominator on basis 'count' rows and fraction on basis 'intensity' rows. " +
        "MetaMorpheus prints a count cell's fraction to two decimals and an intensity cell's " +
        "numerator and denominator to four significant digits (ModificationOccupancyCell). Exact " +
        "intensities are in read-protein-groups.",
        "entity_index is the position among the cell's entities THAT HAVE SITES: an entity with none " +
        "is dropped from the cell, so it cannot be zipped with the group's accessions by position " +
        "(ModificationOccupancyCell).",
        "position counts 0 as the N-terminus and residues from 1; the C-terminus is length + 1 and " +
        "cannot be recognised without the sequence, so only is_n_terminus is given (OccupancySite).",
        "A cell the writer cut short, or replaced with 'Output too long for Excel', keeps its complete " +
        "sites with cell_is_truncated true; a cut cell with no complete site gives no rows. " +
        "truncated_cell_count counts both, so a short table is never mistaken for a complete one.",
        "offset and limit count GROUPS; each group yields one row per sample group, basis and site, " +
        "and a group with no modified sites yields none.",
    ];
}
