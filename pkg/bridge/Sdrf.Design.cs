using MassSpectrometry;
using Readers;

namespace MzLibBridge;

/// <summary>
/// <c>sdrf design</c>: a label-free experimental design read out of an SDRF by mzLib's
/// <see cref="SdrfLabelFreeDesign"/>, or the list of reasons it was refused.
/// </summary>
/// <remarks>
/// <para>
/// <b>A refusal is a result, not an error.</b> mzLib's reader never throws for a bad document: it
/// returns <see cref="SdrfLabelFreeDesign.IsValid"/> = false and every reason at once in
/// <see cref="SdrfLabelFreeDesign.Refusals"/>. So the verb succeeds with <c>is_valid</c> = false,
/// exactly as <c>sdrf validate</c> succeeds on an invalid document. Only a call the bridge cannot
/// make (no file, an <c>--out</c> that is not <c>.tsv</c>) is a usage error.
/// </para>
/// <para>
/// <b>Coordinates are mzLib's model, 0-based.</b> The rows are the <see cref="SpectraFileInfo"/>
/// list mzLib built, so they feed <c>quant flashlfq</c> and <c>quant median-polish</c> unchanged.
/// The file <c>--out</c> writes is <c>ExperimentalDesign.tsv</c>, which is 1-based; mzLib converts in
/// <see cref="SdrfLabelFreeDesign.WriteExperimentalDesignTsv"/> and nowhere else, and so does this
/// verb: it never adds 1 itself.
/// </para>
/// <para>
/// <b>One document per call, by design.</b> The condition columns and the searched files are facts
/// about one deposit, so a batch would have to share them or take a list per document. Neither is
/// specified yet (the spec's <c>open_questions</c>), so there is no <c>--paths-stdin</c>.
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    private static readonly string[] DesignColumns =
        { "full_path", "file_name", "condition", "biological_replicate", "technical_replicate", "fraction" };

    /// <summary>
    /// <c>sdrf design --path FILE [--condition-columns NAMES] [--searched-files-stdin] [--out FILE.tsv]</c>.
    /// <c>NAMES</c> is one or more exact <c>factor value[...]</c> column names separated by tabs (the
    /// one character an SDRF column name cannot hold). With <c>--searched-files-stdin</c>, stdin
    /// holds the files the search will read, one per line, as paths or bare names.
    /// </summary>
    public static object Design(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        // PYB-1: an --out that is not .tsv is refused before anything is read, so a typo never costs
        // a read and never leaves a file of the wrong kind behind.
        string? outputPath = arguments.Optional("out");
        if (arguments.WasProvided("out") && string.IsNullOrWhiteSpace(outputPath))
            throw new Program.UsageException("Option --out needs a file path ending in .tsv, e.g. 'ExperimentalDesign.tsv'.");
        if (outputPath is not null && !outputPath.EndsWith(".tsv", StringComparison.OrdinalIgnoreCase))
            throw new Program.UsageException(
                $"Option --out writes ExperimentalDesign.tsv's format, so the path must end in .tsv; got '{outputPath}'.");

        if (!File.Exists(path))
            throw new Program.UsageException($"SDRF file not found: '{path}'.");

        List<string>? conditionColumns = null;
        if (arguments.WasProvided("condition-columns"))
        {
            conditionColumns = (arguments.Optional("condition-columns") ?? string.Empty)
                .Split('\t')
                .Select(c => c.Trim())
                .ToList();
            if (conditionColumns.Count == 0 || conditionColumns.Any(c => c.Length == 0))
                throw new Program.UsageException(
                    "Option --condition-columns needs one or more factor value[...] column names, " +
                    "separated by tabs, with none blank.");
        }

        List<string>? searchedFiles = arguments.Flag("searched-files-stdin")
            ? Program.ReadStdinLines()
            : null;

        var options = new SdrfLabelFreeDesignOptions
        {
            ConditionColumns = conditionColumns,
            SearchedFiles = searchedFiles,
        };
        SdrfLabelFreeDesign design = SdrfLabelFreeDesign.Read(path, options);

        var table = new ColumnTable(DesignColumns);
        foreach (SpectraFileInfo file in design.Files)
            table.Add(file.FullFilePathWithExtension, file.FilenameWithoutExtension, file.Condition,
                file.BiologicalReplicate, file.TechnicalReplicate, file.Fraction);

        object? written = null;
        if (outputPath is not null && design.IsValid)
        {
            // mzLib's writer throws on a refused design; the check above means it is never asked to,
            // so a refusal leaves no file for MetaMorpheus to find.
            design.WriteExperimentalDesignTsv(outputPath);
            written = new { path = outputPath, file_count = design.Files.Count };
        }

        return new
        {
            path,
            is_valid = design.IsValid,
            file_key_column = design.FileKeyColumn,
            condition_columns = design.ConditionColumns.ToList(),
            condition_columns_declared = conditionColumns is not null,
            searched_files_given = searchedFiles is not null,
            file_count = design.Files.Count,
            refusals = design.Refusals.ToList(),
            notes = design.Notes.ToList(),
            report = design.Report(),
            column_names = table.Names,
            columns = table.Columns(),
            written,
            caveats = DesignCaveats(),
        };
    }

    private static List<string> DesignCaveats() => new()
    {
        "LABEL-FREE ONLY. mzLib's SdrfLabelFreeDesign builds MetaMorpheus's label-free " +
        "ExperimentalDesign.tsv; an isobaric (TMT, iTRAQ) SDRF needs a channel design this verb does " +
        "not produce.",

        "is_valid = false is an answer, not a failure: refusals lists EVERY reason at once and the " +
        "table is empty. mzLib refuses rather than repairs because MetaMorpheus skips quantification " +
        "with only a warning when its design file is invalid.",

        "biological_replicate, technical_replicate and fraction are 0-based, as mzLib's " +
        "SpectraFileInfo and the quant verbs take them. ExperimentalDesign.tsv (written by --out) is " +
        "1-based; mzLib converts when it writes.",

        "notes records every relabelling: biological replicates ranked within each condition (with " +
        "the mapping), and rows dropped because the search does not read their file. Read it; it is " +
        "the only record of a renumbering.",

        "With no condition columns declared, mzLib uses the document's only factor value column and " +
        "refuses a document that has several. A factor holding 'not available' or 'not applicable' " +
        "is refused, because it would pair samples nobody said were alike.",

        "MetaMorpheus finds the design only as a file named ExperimentalDesign.tsv beside the " +
        "spectra. --out writes any .tsv path; naming and placing it is the caller's job.",
    };
}
