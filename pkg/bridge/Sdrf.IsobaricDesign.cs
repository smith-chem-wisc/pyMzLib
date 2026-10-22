using Omics.Modifications;
using Readers;

namespace MzLibBridge;

/// <summary>
/// <c>sdrf isobaric-design</c>: an isobaric (TMT, TMTpro, iTRAQ, DiLeu) experimental design read out
/// of an SDRF by mzLib's <see cref="SdrfIsobaricDesign"/>, or every reason it was refused. The
/// channel-design counterpart of <c>sdrf design</c> (<c>Sdrf.Design.cs</c>).
/// </summary>
/// <remarks>
/// <para>
/// <b>A refusal is a result, not an error.</b> mzLib's reader never throws for a bad document or bad
/// options: a missing kit, no plex source, two plex sources, a bad regular expression and a channel
/// the kit does not have are all <see cref="SdrfIsobaricDesign.Refusals"/>. So the verb succeeds with
/// <c>is_valid</c> = false. Only a call the bridge cannot make is a usage error: no file, a kit name
/// mzLib cannot resolve (there is then no <see cref="IsobaricMassTag"/> to hand it), an option given
/// with no value.
/// </para>
/// <para>
/// <b>1-based, as written.</b> mzLib keeps an isobaric design's fraction, technical replicate and
/// biological replicate exactly as the SDRF writes them, 1-based, which is also what MetaMorpheus's
/// <c>TmtDesign.txt</c> holds. <c>sdrf design</c> is 0-based, because its rows are mzLib's
/// <see cref="MassSpectrometry.SpectraFileInfo"/>. That split is mzLib's, and this verb converts
/// nothing.
/// </para>
/// <para>
/// <b>The kit is resolved by <c>isobaric kits</c>' own code path</b> (<see cref="Isobaric.ResolveKit"/>),
/// so a name one verb accepts the other never refuses.
/// </para>
/// <para>
/// <b><c>--out</c> names a folder.</b> The file's name is MetaMorpheus's (<c>TmtDesign.txt</c>,
/// <c>GlobalVariables.TmtExperimentalDesignFileName</c>) and its format mzLib's
/// (<see cref="SdrfIsobaricDesign.WriteTmtDesign"/>), so neither is the bridge's to choose. mzLib
/// holds no constant for the name; <see cref="TmtDesignFileName"/> is the one copy here, as
/// <c>quant median-polish</c> holds <c>QuantifiedProteins.tsv</c>, until mzLib publishes one.
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    /// <summary>
    /// MetaMorpheus's name for the file <see cref="SdrfIsobaricDesign.WriteTmtDesign"/> writes
    /// (<c>GlobalVariables.TmtExperimentalDesignFileName</c>). MetaMorpheus reads it only from the
    /// folder of the first spectra file. Not an mzLib constant yet: a PUSH-TO-MAINLAND item.
    /// </summary>
    internal const string TmtDesignFileName = "TmtDesign.txt";

    private static readonly string[] IsobaricDesignColumns =
    {
        "full_path", "plex", "plex_id", "sample_name", "channel", "condition",
        "biological_replicate", "fraction", "technical_replicate",
    };

    /// <summary>
    /// <c>sdrf isobaric-design --path FILE [--kit NAME] [--plex-column NAME] [--plex-file-name-pattern
    /// REGEX] [--single-plex NAME] [--condition-columns NAMES] [--searched-files-stdin] [--out DIR]</c>.
    /// Every option is passed to <see cref="SdrfIsobaricDesignOptions"/> as given, or not at all.
    /// </summary>
    public static object IsobaricDesign(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        // --out names a folder, checked before anything is read, so a typo never costs a read.
        string? outputDirectory = arguments.Optional("out");
        if (arguments.WasProvided("out") && string.IsNullOrWhiteSpace(outputDirectory))
            throw new Program.UsageException(
                $"Option --out needs a folder; {TmtDesignFileName} is written inside it, e.g. '--out design'.");
        if (outputDirectory is not null && File.Exists(outputDirectory))
            throw new Program.UsageException(
                $"Option --out names a folder, and '{outputDirectory}' is a file. {TmtDesignFileName} " +
                "(MetaMorpheus's name) is written inside the folder you name.");

        if (!File.Exists(path))
            throw new Program.UsageException($"SDRF file not found: '{path}'.");

        // Absent, mzLib is given no tag and refuses ("No isobaric tag was given"): the kit is never
        // guessed from the labels, by mzLib or here.
        string? kit = arguments.Optional("kit");
        IsobaricMassTag? tag = arguments.WasProvided("kit") ? Isobaric.ResolveKit(kit) : null;

        string? plexColumn = OptionalText(arguments, "plex-column", "a column name, e.g. 'comment[plex]'");
        string? plexPattern = OptionalText(arguments, "plex-file-name-pattern", "a regular expression, e.g. 'TMT_?pool(?<plex>\\d+)'");
        string? singlePlex = OptionalText(arguments, "single-plex", "a plex name, e.g. 'plex1'");
        List<string>? conditionColumns = ConditionColumnsFrom(arguments);
        List<string>? searchedFiles = arguments.Flag("searched-files-stdin")
            ? Program.ReadStdinLines()
            : null;

        var options = new SdrfIsobaricDesignOptions
        {
            Tag = tag,
            PlexColumn = plexColumn,
            PlexFileNamePattern = plexPattern,
            SinglePlex = singlePlex,
            ConditionColumns = conditionColumns,
            SearchedFiles = searchedFiles,
        };
        SdrfIsobaricDesign design = SdrfIsobaricDesign.Read(path, options);

        IReadOnlyDictionary<string, int> plexIds = design.PlexIds;
        var table = new ColumnTable(IsobaricDesignColumns);
        foreach (SdrfIsobaricFile file in design.Files)
            foreach (SdrfIsobaricChannel channel in file.Channels)
                table.Add(file.FilePath, file.Plex, plexIds[file.Plex], channel.SampleName, channel.Label,
                    channel.Condition, channel.BiologicalReplicate, file.Fraction, file.TechnicalReplicate);

        object? written = null;
        if (outputDirectory is not null && design.IsValid)
        {
            // mzLib's writer throws on a refused design; the check above means it is never asked to,
            // so a refusal leaves no file for MetaMorpheus to find. The folder is created (PYB-5).
            Directory.CreateDirectory(outputDirectory);
            string outputPath = Path.Combine(outputDirectory, TmtDesignFileName);
            design.WriteTmtDesign(outputPath);
            written = new { path = outputPath, file_count = design.Files.Count, record_count = table.RowCount };
        }

        return new
        {
            path,
            is_valid = design.IsValid,
            kit,
            kit_type = design.Tag?.TagType.ToString(),
            kit_channels = design.Tag?.ChannelLabels.ToList(),
            plex_column = plexColumn,
            plex_file_name_pattern = plexPattern,
            single_plex = singlePlex,
            plex_source = design.PlexSource,
            plex_ids = plexIds.OrderBy(p => p.Value).ToDictionary(p => p.Key, p => p.Value),
            plex_count = plexIds.Count,
            file_key_column = design.FileKeyColumn,
            condition_columns = design.ConditionColumns.ToList(),
            condition_columns_declared = conditionColumns is not null,
            searched_files_given = searchedFiles is not null,
            file_count = design.Files.Count,
            record_count = table.RowCount,
            refusals = design.Refusals.ToList(),
            notes = design.Notes.ToList(),
            report = design.Report(),
            column_names = table.Names,
            columns = table.Columns(),
            written,
            caveats = IsobaricDesignCaveats(),
        };
    }

    /// <summary>
    /// A string option passed through as given: null when absent, a usage error when given with no
    /// value (a blank would otherwise read to mzLib as "not declared", which is not what was typed).
    /// </summary>
    private static string? OptionalText(Program.Arguments arguments, string name, string example)
    {
        if (!arguments.WasProvided(name))
            return null;
        string? value = arguments.Optional(name);
        if (string.IsNullOrWhiteSpace(value))
            throw new Program.UsageException($"Option --{name} needs {example}.");
        return value;
    }

    private static List<string> IsobaricDesignCaveats() => new()
    {
        "1-BASED, AS WRITTEN. fraction, technical_replicate and biological_replicate are the SDRF's " +
        "own numbers, as mzLib's SdrfIsobaricDesign keeps them and MetaMorpheus's TmtDesign.txt " +
        "holds them. sdrf design is 0-based (mzLib's SpectraFileInfo). The split is mzLib's; the " +
        "bridge converts neither table, so never mix the two without adding or removing 1.",

        "is_valid = false is an answer, not a failure: refusals lists EVERY reason at once and the " +
        "table is empty. A missing kit, no plex source, two plex sources and a bad pattern are " +
        "refusals too. mzLib refuses rather than repairs because MetaMorpheus quantifies nothing, " +
        "with only a warning, from a design it rejects.",

        "Nothing is guessed. The kit is only the one named, because one kit's labels can be a " +
        "subset of another's (TMT11 with 131C unannotated reads exactly like TMT10). The plex is " +
        "only the one source declared: comment[sample preparation batch] is often not a plex " +
        "(in PXD061609 it holds three fractionation arms of one plex).",

        "One row per file and ANNOTATED channel, in SDRF row order and then reporter m/z order. A " +
        "kit channel the plex does not annotate is absent; mzLib's quantification design " +
        "(ToExperimentalDesign) carries it with biological replicate 0 and no sample.",

        "A plex is one channel-to-sample map: every file of a plex carries all of its channels, so " +
        "fractions share their samples. notes names each file that has no row of its own for one " +
        "of its plex's channels.",

        "Sample type is not read (mzLib MAP-08): every annotated channel is a study sample, so a " +
        "reference or pool channel is a sample under the name the SDRF gives it.",

        "--out names a FOLDER. TmtDesign.txt (MetaMorpheus's name) is written inside it, 1-based, " +
        "replacing one already there (MetaMorpheus's own --sdrfDesign never replaces one). " +
        "MetaMorpheus reads it only from the folder of the first spectra file, and refuses a " +
        "folder that also holds ExperimentalDesign.tsv.",

        "TMT16 is mzLib's TMTpro 16-plex. MetaMorpheus searches it with its TMT18 label; a TMT16 " +
        "design is read as is by such a run, with channels 134C and 135N empty.",
    };
}
