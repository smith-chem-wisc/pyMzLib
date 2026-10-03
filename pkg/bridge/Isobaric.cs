using Omics.Modifications;

namespace MzLibBridge;

/// <summary>
/// <c>isobaric kits</c>: the isobaric labelling kits mzLib can name (TMT, TMTpro, iTRAQ, DiLeu), with
/// every channel's label and the theoretical m/z of its reporter ion, from
/// <see cref="IsobaricMassTag"/>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Nothing here is typed in.</b> mzLib holds only the channel labels; each m/z is a <c>DI HCD:</c>
/// diagnostic-ion line of the kit's <c>Multiplex Label</c> modification in mzLib's embedded
/// <c>TMT.txt</c>, plus one proton, sorted ascending and paired with the labels by position. The
/// bridge copies those arrays into a table and adds nothing.
/// </para>
/// <para>
/// <b>Its own module.</b> The kits are labelling chemistry, not a property of one peptidoform or one
/// file, so they sit neither under <c>peptidoform</c> nor under <c>readers</c>. Reading reporter
/// intensities out of spectra (<see cref="IsobaricMassTag.GetReporterIonIntensities"/>) would be the
/// module's second verb; it is not built.
/// </para>
/// </remarks>
internal static class Isobaric
{
    private static readonly string[] KitColumns =
        { "kit", "channel_index", "channel_label", "reporter_ion_mz", "mz_min", "mz_max" };

    /// <summary>
    /// <c>isobaric kits [--kit NAME]</c>: every kit, or the one <c>NAME</c> resolves to. <c>NAME</c>
    /// is matched by mzLib's <see cref="IsobaricMassTag.TryGetTagType"/>: the whole name,
    /// case-insensitively, optionally followed by <c>" on &lt;motif&gt;"</c>, never a substring.
    /// </summary>
    public static object Kits(Program.Arguments arguments)
    {
        string? requested = null;
        List<IsobaricMassTagType> types;
        if (arguments.WasProvided("kit"))
        {
            requested = arguments.Optional("kit");
            if (string.IsNullOrWhiteSpace(requested))
                throw new Program.UsageException("Option --kit needs a kit name, e.g. 'TMT10' or 'iTRAQ-4plex'.");
            if (!IsobaricMassTag.TryGetTagType(requested, out IsobaricMassTagType type))
                throw new Program.UsageException(
                    $"mzLib knows no isobaric kit named '{requested}'. Known kits: " +
                    $"{string.Join(", ", Enum.GetNames<IsobaricMassTagType>())}; MetaMorpheus modification " +
                    "names such as 'TMT6-plex' or 'iTRAQ-4plex on K' are accepted too. The name must match " +
                    "whole: 'TMT10plex' is not 'TMT10'.");
            types = new List<IsobaricMassTagType> { type };
        }
        else
        {
            types = Enum.GetValues<IsobaricMassTagType>().ToList();
        }

        var table = new Sdrf.ColumnTable(KitColumns);
        var kits = new List<object>(types.Count);
        foreach (IsobaricMassTagType type in types)
        {
            if (!IsobaricMassTag.TryGetIsobaricMassTag(type, out IsobaricMassTag? tag) || tag is null)
                throw new InvalidOperationException(
                    $"mzLib could not build the {type} kit from its embedded modifications (its Multiplex " +
                    "Label entry is missing, has no HCD diagnostic ions, or has the wrong ion count).");

            for (int i = 0; i < tag.ReporterIonMzs.Length; i++)
                table.Add(type.ToString(), i, tag.ChannelLabels[i], tag.ReporterIonMzs[i],
                    tag.ReporterIonMzRanges[i].Minimum, tag.ReporterIonMzRanges[i].Maximum);

            kits.Add(new { kit = type.ToString(), channel_count = tag.ReporterIonMzs.Length });
        }

        return new
        {
            kit = requested,
            kit_count = kits.Count,
            record_count = table.RowCount,
            absolute_tolerance = IsobaricMassTag.AbsoluteToleranceValue,
            kits,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = KitCaveats(),
        };
    }

    private static List<string> KitCaveats() => new()
    {
        "reporter_ion_mz is THEORETICAL, charge 1: a DI HCD diagnostic-ion mass from mzLib's embedded " +
        "TMT.txt plus one proton. It is not a calibrated or observed value.",

        "mz_min and mz_max are reporter_ion_mz -/+ absolute_tolerance (0.003 Da, the window mzLib " +
        "takes from doi:10.1021/acs.jproteome.1c00168). mzLib uses the same window to read reporter " +
        "intensities from a spectrum, keeping the most intense peak inside it.",

        "TMT16 and TMT18 are TMTpro. TMT16 is the lowest sixteen channels of the TMTpro 18-plex set " +
        "and has no modification entry of its own, so its m/z are TMT18's first sixteen.",

        "iTRAQ8 has no 120 channel: the phenylalanine immonium ion sits at 120.081, so the eighth " +
        "reagent is 121.",

        "channel_index is 0-based and ascends with reporter_ion_mz; channel_label at that index names " +
        "the channel whose reporter ion it is.",
    };
}
