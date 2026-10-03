using Omics.SequenceConversion;
using Readers.ProForma;

namespace MzLibBridge;

/// <summary>
/// <c>peptidoform convert</c>: rewrite full sequences from one notation to another with mzLib's
/// <see cref="SequenceConversionService"/>, for example a MetaMorpheus full sequence
/// <c>[UniProt:N-acetylserine on S]SEQK</c> to <c>[UNIMOD:1]SEQK</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Projection only.</b> Every conversion is <see cref="SequenceConversionService.Default"/>'s
/// <see cref="ISequenceConverter.Convert"/>, with the format names and the
/// <see cref="SequenceConversionHandlingMode"/> passed through as given. The bridge validates the
/// two names against the service's own registry, reads mzLib's <see cref="ConversionWarnings"/>
/// back out, and shapes one row per input. It resolves no modification itself and picks no
/// lookup: which modifications a target can name is the serializer's business, and the ProForma
/// target at the pin is known not to name UniProt-sourced ones (see the caveats).
/// </para>
/// <para>
/// <b>ProForma is registered first.</b> mzLib keeps ProForma in <c>Readers</c> and registers it at
/// runtime through <see cref="ProFormaSequenceConversion.RegisterWithDefault"/>, which is idempotent
/// and is what mzLib's own psmtsv reader calls. Without it the Default service would not list
/// ProForma at all.
/// </para>
/// </remarks>
internal static partial class Peptidoform
{
    private static readonly string[] ConvertColumns =
        { "input", "output", "status", "failure_reason", "incompatible_items", "warnings", "errors" };

    /// <summary>
    /// <c>peptidoform convert [--from mzLib] [--to Unimod] [--mode ReturnNull] [--threads 1]</c>,
    /// with the sequences on stdin, one per line.
    /// </summary>
    public static object Convert(Program.Arguments arguments)
    {
        ProFormaSequenceConversion.RegisterWithDefault();
        SequenceConversionService service = SequenceConversionService.Default;

        string source = FormatName(arguments, "from", "mzLib", service.AvailableSourceFormats, "source");
        string target = FormatName(arguments, "to", "Unimod", service.AvailableTargetFormats, "target");
        SequenceConversionHandlingMode mode = ModeFrom(arguments);
        int threads = Proteins.ThreadsFrom(arguments);

        List<string> inputs = Program.ReadStdinLines();
        if (inputs.Count == 0)
            throw new Program.UsageException(
                "No sequences on stdin: send one full sequence per line, e.g. '[UniProt:N-acetylserine on S]SEQK'.");

        // Both names were checked above, so the service hands back a converter: a registered one
        // (mzLib-Unimod, mzLib-ProForma, Modomics-mzLib, ...) or the parser x serializer pair it
        // composes itself, exactly as SequenceConversionService.Convert would.
        ISequenceConverter converter = service.GetConverter(source, target)
            ?? throw new InvalidOperationException($"mzLib registered no converter from {source} to {target}.");

        var rows = new Row[inputs.Count];
        var options = new ParallelOptions { MaxDegreeOfParallelism = threads };
        try
        {
            Parallel.For(0, inputs.Count, options, i => rows[i] = ConvertOne(converter, inputs[i], i, mode));
        }
        catch (AggregateException ex) when (ex.InnerExceptions.All(e => e is Program.UsageException))
        {
            // ThrowException: report the first failing input in input order, whatever --threads is.
            throw ex.InnerExceptions.Cast<Program.UsageException>()
                .OrderBy(e => (int)e.Data["index"]!).First();
        }

        var table = new Sdrf.ColumnTable(ConvertColumns);
        int converted = 0, warned = 0, failed = 0;
        foreach (Row row in rows)
        {
            table.Add(row.Input, row.Output, row.Status, row.FailureReason, row.IncompatibleItems, row.Warnings, row.Errors);
            switch (row.Status)
            {
                case StatusConverted: converted++; break;
                case StatusWarned: warned++; break;
                default: failed++; break;
            }
        }

        return new
        {
            source_format = source,
            target_format = target,
            mode = mode.ToString(),
            source_formats = Sorted(service.AvailableSourceFormats),
            target_formats = Sorted(service.AvailableTargetFormats),
            record_count = table.RowCount,
            converted_count = converted,
            warned_count = warned,
            failed_count = failed,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = ConvertCaveats(),
        };
    }

    internal const string StatusConverted = "converted";
    internal const string StatusWarned = "converted_with_warnings";
    internal const string StatusFailed = "failed";

    private sealed record Row(string Input, string? Output, string Status, string? FailureReason,
        List<string> IncompatibleItems, List<string> Warnings, List<string> Errors);

    /// <summary>One input through mzLib's converter, with mzLib's own account of what went wrong.</summary>
    private static Row ConvertOne(ISequenceConverter converter, string input, int index, SequenceConversionHandlingMode mode)
    {
        var warnings = new ConversionWarnings();
        string? output;
        ConversionFailureReason? reason;
        List<string> incompatible;
        bool clean;
        try
        {
            output = converter.Convert(input, warnings, mode);
            reason = warnings.FailureReason;
            incompatible = warnings.IncompatibleItems.ToList();
            clean = warnings.IsClean;
        }
        catch (SequenceConversionException ex) when (mode == SequenceConversionHandlingMode.ThrowException)
        {
            var usage = new Program.UsageException(
                $"mzLib could not convert '{input}' from {converter.SourceFormatName} to {converter.TargetFormatName} " +
                $"({ex.FailureReason}): {ex.Message} Use mode ReturnNull to get one row per sequence instead.");
            usage.Data["index"] = index;
            throw usage;
        }
        catch (SequenceConversionException ex)
        {
            // mzLib threw although the mode asked it not to. It is still mzLib's verdict on this
            // input, so it is reported on the row rather than failing every other sequence.
            output = null;
            reason = ex.FailureReason;
            incompatible = (ex.IncompatibleItems ?? Array.Empty<string>()).Concat(warnings.IncompatibleItems).Distinct().ToList();
            warnings.AddError(ex.Message);
            clean = false;
        }

        // mzLib's own verdict, read straight off ConversionWarnings: no output is a failure, and an
        // output with anything recorded against it (a removed or skipped modification, a skipped
        // character) is not clean.
        string status = output is null ? StatusFailed
            : clean ? StatusConverted
            : StatusWarned;

        return new Row(input, output, status, reason?.ToString(), incompatible,
            warnings.Warnings.ToList(), warnings.Errors.ToList());
    }

    private static string FormatName(Program.Arguments arguments, string option, string fallback,
        IReadOnlyCollection<string> registered, string role)
    {
        if (arguments.WasProvided(option) && string.IsNullOrWhiteSpace(arguments.Optional(option)))
            throw new Program.UsageException($"Option --{option} was given without a value.");
        string requested = arguments.Optional(option) ?? fallback;

        // The service matches names case-insensitively; the wire reports the name mzLib registered.
        string? match = registered.FirstOrDefault(n => string.Equals(n, requested, StringComparison.OrdinalIgnoreCase));
        return match ?? throw new Program.UsageException(
            $"mzLib has no {role} format named '{requested}'. Registered {role} formats: " +
            $"{string.Join(", ", Sorted(registered))}.");
    }

    private static SequenceConversionHandlingMode ModeFrom(Program.Arguments arguments)
    {
        if (arguments.WasProvided("mode") && string.IsNullOrWhiteSpace(arguments.Optional("mode")))
            throw new Program.UsageException("Option --mode was given without a value.");
        string? raw = arguments.Optional("mode");
        if (raw is null)
            return SequenceConversionHandlingMode.ReturnNull;
        if (Enum.TryParse(raw, ignoreCase: true, out SequenceConversionHandlingMode mode)
            && Enum.IsDefined(mode) && !int.TryParse(raw, out _))
            return mode;
        throw new Program.UsageException(
            $"Option --mode must be one of mzLib's SequenceConversionHandlingMode names: " +
            $"{string.Join(", ", Enum.GetNames<SequenceConversionHandlingMode>())}; got '{raw}'.");
    }

    private static List<string> Sorted(IEnumerable<string> names) =>
        names.OrderBy(n => n, StringComparer.Ordinal).ToList();

    private static List<string> ConvertCaveats() => new()
    {
        "ProForma does not resolve UniProt-sourced modifications at the bridge's mzLib pin: " +
        "'[UniProt:N-acetylserine on S]SEQK' comes out as '[UniProt:N-acetylserine on S]-SEQK', status " +
        "converted, not '[UNIMOD:1]-SEQK'. mzLib's ProForma serializer looks modifications up in the " +
        "MetaMorpheus list only (ProFormaSequenceSerializer.cs:34, MzLibModificationLookup.cs:54); the " +
        "Unimod target looks them up in every list mzLib loads, UniProt's included " +
        "(GlobalModificationLookup.cs:23). Convert to Unimod when the sequences carry UniProt modifications.",

        "ProForma writes a modification it cannot resolve under its mzLib name, and the row is still " +
        "converted: status converted does not mean every bracket in a ProForma output is a controlled " +
        "vocabulary term. Unimod either names a modification by its UNIMOD accession or reports it in " +
        "incompatible_items.",

        "What a modification mzLib cannot write in the target does to a row depends on mode. ReturnNull: " +
        "the row fails (output null). RemoveIncompatibleElements and UsePrimarySequence: the modification " +
        "is dropped from output, named in incompatible_items, and the row is converted_with_warnings. " +
        "ThrowException: the first such input fails the whole call as a usage error naming it.",

        "status is mzLib's own verdict: failed when mzLib returned no output, converted_with_warnings " +
        "when it returned one but recorded a warning, error or incompatible item (ConversionWarnings." +
        "IsClean false), converted otherwise. failure_reason is mzLib's ConversionFailureReason when it " +
        "recorded one; the Unimod serializer under ReturnNull fails a row without recording a reason, so " +
        "read incompatible_items.",

        "An ambiguous MetaMorpheus full sequence ('|'-joined candidates) is not refused: mzLib's parser " +
        "skips each '|' with a warning and joins the candidates into one sequence, status " +
        "converted_with_warnings. Split ambiguous full sequences before converting them.",

        "Each line of stdin is one input, given to mzLib exactly as read; blank lines are skipped. Rows " +
        "are in input order, duplicates included, at any --threads.",
    };
}
