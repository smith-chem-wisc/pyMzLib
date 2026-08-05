using System.Collections;
using System.Globalization;
using System.Reflection;
using Chemistry;
using MassSpectrometry;
using Readers;

namespace MzLibBridge;

/// <summary>
/// The <i>exhaustive</i> half of the readers surface: every file type mzLib recognises, read into a
/// table — plus the three typed views that only some of them have.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="Reading.ReadResults"/> reads the <c>quantifiable</c> view, which exactly three of
/// mzLib's twenty-nine file types implement. That is a real and useful view — it is what FlashLFQ
/// consumes — but it left twenty-six formats readable by mzLib and unreachable from a binding. The
/// verbs here close that gap, in two different ways, because the gap has two different shapes.
/// </para>
/// <para>
/// <b>The typed views</b> — <c>read-features</c> (<see cref="IMs1FeatureFile"/>) and
/// <c>read-matches</c> (<see cref="ISpectralMatch"/>) and <c>read-spectra</c>
/// (<see cref="MsDataFile"/>) — are cross-format like <c>read-results</c>: a fixed column set, the
/// same for every format that offers the view, safe to compare between files. They add nine of the
/// twenty-six.
/// </para>
/// <para>
/// <b><c>read-records</c> is the exhaustive one</b>, and it is a deliberately different animal. It
/// works on all twenty-nine, including the thirteen that belong to no cross-format interface at all
/// (TopPIC, Crux, MSFragger's peptide and protein tables, the FlashDeconv formats, …), by
/// projecting each format's <i>own</i> record type. So its columns are <b>not</b> uniform: reading a
/// TopPIC file gives you TopPIC's thirty-odd columns under TopPIC's own names, and reading a Crux
/// file gives you Crux's. That is the honest shape of the data — mzLib does not normalise these
/// formats onto a common record, and inventing a normalisation here would be the bridge answering a
/// question mzLib cannot. The column names are therefore reported in every response, and a caller
/// who wants comparable numbers across formats wants a typed view, not this verb.
/// </para>
/// <para>
/// What <c>read-records</c> buys, in exchange for giving up uniformity, is that <b>no field is
/// silently dropped</b>. Every property mzLib parsed is either a column or is named in
/// <c>excluded_fields</c> with the reason it could not cross the wire. A format-specific column that
/// simply vanished would be indistinguishable from one the file did not contain.
/// </para>
/// </remarks>
internal static partial class Reading
{
    /// <summary>
    /// <c>readers read-records --path FILE [--limit N] [--offset N] [--out FILE]</c> — any of the
    /// 29 file types, as a table of its own native fields.
    /// </summary>
    /// <remarks>
    /// The only verb here with no view requirement: if <c>readers identify</c> succeeds on a path,
    /// this reads it. See the class remarks for why the columns are per-format rather than uniform.
    /// </remarks>
    public static object ReadRecords(Program.Arguments arguments)
    {
        Window window = Window.From(arguments);
        string path = window.Path;

        IResultFile resultFile = OpenAny(path);

        // MsDataFileToResultFileAdapter.Results is a plain auto-property that stays null until
        // LoadResults() is called, unlike ResultFile<T>.Results which lazy-loads on first access.
        // FileReader.ReadResultFile never calls it, so reading .Results off a spectra file without
        // this line is a NullReferenceException rather than an empty list. Calling it on a
        // ResultFile<T> is a harmless re-parse of a file that is about to be parsed anyway.
        resultFile.LoadResults();

        IReadOnlyList<object> all = RecordsOf(resultFile);
        IReadOnlyList<object> selected = window.Apply(all, out bool truncated);

        RecordProjection projection = RecordProjection.For(ElementTypeOf(resultFile, all));

        object? written = null;
        if (window.OutputPath is not null)
            written = WriteProjectedTable(window.OutputPath, projection, selected);

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            reader = resultFile.GetType().Name,
            // The record type's own name, e.g. "ToppicPrsm". The columns belong to THIS type, and
            // saying which one they came from is what makes them cross-referenceable against the
            // mzLib source rather than a bare list of strings.
            record_type = projection.RecordTypeName,
            views = ViewsOf(resultFile.GetType()),
            record_count = all.Count,
            returned_count = written is null ? selected.Count : 0,
            offset = window.Offset,
            truncated,
            // Every property that exists on the record type and could NOT become a column, with the
            // reason. A dropped field must never be indistinguishable from an absent one.
            excluded_fields = projection.Excluded,
            // Properties that threw when read, with the exception type. mzLib has several computed
            // properties that assume a UniProt-style header and throw on anything else — Crux's
            // Accession is `ProteinId.Split('|')[1]`. Those become null cells rather than a failed
            // read of the whole file, but silence would misreport a parse failure as absent data.
            failed_fields = projection.FailedFieldsFor(selected),
            column_names = projection.ColumnNames,
            columns = written is null ? projection.BuildColumns(selected) : null,
            output = written,
        };
    }

    /// <summary>
    /// <c>readers read-features --path FILE [--limit N] [--offset N] [--out FILE]</c> — the
    /// deconvolved-MS1-feature view.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The <c>ms1_features</c> view: <see cref="IMs1FeatureFile.GetMs1Features"/> projected onto
    /// <see cref="ISingleChargeMs1Feature"/>. Two file types offer it — TopFD/FLASHDeconv
    /// <c>_ms1.feature</c> and Dinosaur <c>.feature.tsv</c>.
    /// </para>
    /// <para>
    /// <b>One row here is not one row in the file.</b> An <c>_ms1.feature</c> row is a
    /// deconvolved neutral mass spanning a charge range, and mzLib expands it into one
    /// single-charge feature per charge in <c>[ChargeStateMin, ChargeStateMax]</c>. So a file of a
    /// hundred features can read as a thousand rows, and <c>record_count</c> counts the expansion,
    /// not the lines. That is mzLib's model and the one the deconvolution consumers use; the raw
    /// per-file rows are available through <c>read-records</c>.
    /// </para>
    /// </remarks>
    public static object ReadFeatures(Program.Arguments arguments)
    {
        Window window = Window.From(arguments);
        string path = window.Path;

        IResultFile resultFile = OpenAny(path);
        if (resultFile is not IMs1FeatureFile featureFile)
            throw NoSuchView(resultFile, "ms1_features", "read-features");

        resultFile.LoadResults();
        List<Feature> all = FeaturesOf(featureFile, resultFile);
        IReadOnlyList<Feature> selected = window.Apply(all, out bool truncated);

        var columns = FeatureColumns;
        object? written = null;
        if (window.OutputPath is not null)
            written = WriteTable(window.OutputPath, columns, selected);

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            record_count = all.Count,
            returned_count = written is null ? selected.Count : 0,
            offset = window.Offset,
            truncated,
            retention_time_unit = FeatureRetentionTimeUnitOf(resultFile),
            caveats = FeatureCaveatsFor(resultFile, all),
            column_names = columns.Select(c => c.Name).ToList(),
            columns = written is null ? BuildColumns(columns, selected) : null,
            output = written,
        };
    }

    /// <summary>
    /// <c>readers read-matches --path FILE [--limit N] [--offset N] [--out FILE]</c> — the
    /// spectral-match view.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The <c>spectral_match</c> view: <see cref="ISpectralMatch"/>, which MsPathFinderT's three
    /// result types and Casanovo's <c>.mztab</c> implement on their <i>record</i> rather than on
    /// the file — which is why <c>readers identify</c> has to look at the record type to report it.
    /// </para>
    /// <para>
    /// Note that this is <b>not</b> the same interface as the identically-named
    /// <c>Omics.SpectralMatch.ISpectralMatch</c>, which carries a score and no modifications.
    /// mzLib has two unrelated types of that name and aliases them in its own source to tell them
    /// apart. This verb projects the <c>Readers</c> one, and like <c>read-results</c> it carries no
    /// confidence field at all — so nothing here is FDR-filtered.
    /// </para>
    /// </remarks>
    public static object ReadMatches(Program.Arguments arguments)
    {
        Window window = Window.From(arguments);
        string path = window.Path;

        IResultFile resultFile = OpenAny(path);
        resultFile.LoadResults();

        List<ISpectralMatch> all = RecordsOf(resultFile).OfType<ISpectralMatch>().ToList();
        if (all.Count == 0 && !HasSpectralMatchRecords(resultFile.GetType()))
            throw NoSuchView(resultFile, "spectral_match", "read-matches");

        IReadOnlyList<ISpectralMatch> selected = window.Apply(all, out bool truncated);

        var columns = MatchColumns;
        object? written = null;
        if (window.OutputPath is not null)
            written = WriteTable(window.OutputPath, columns, selected);

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            record_count = all.Count,
            returned_count = written is null ? selected.Count : 0,
            offset = window.Offset,
            truncated,
            caveats = MatchCaveatsFor(resultFile.FileType),
            column_names = columns.Select(c => c.Name).ToList(),
            columns = written is null ? BuildColumns(columns, selected) : null,
            output = written,
        };
    }

    /// <summary>
    /// <c>readers read-spectra --path FILE [--limit N] [--offset N] [--ms-order N] [--peaks]
    /// [--out FILE]</c> — the scan headers of a spectra file, and optionally its peaks.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The <c>spectra</c> view: the seven <see cref="MsDataFile"/>-backed types — mzML, MGF,
    /// ms1/ms2 msalign, Thermo <c>.raw</c>, Bruker <c>.d</c>, and timsTOF <c>.d</c>.
    /// </para>
    /// <para>
    /// <b>Peaks are opt-in</b>, and that default is the whole design of this verb. A scan header is
    /// tens of bytes; a scan's peak list is thousands, and a mid-size mzML holds tens of thousands
    /// of scans. Returning peaks by default would make the ordinary call — "what is in this
    /// file?" — serialise hundreds of megabytes through a JSON envelope. With <c>--peaks</c> the
    /// <c>mz</c> and <c>intensity</c> columns each carry one array per scan, so use it with
    /// <c>--limit</c> or a scan filter, or write to <c>--out</c>, where peaks are emitted as
    /// <c>;</c>-joined lists.
    /// </para>
    /// <para>
    /// <c>--ms-order</c> filters to MS1, MS2, … <i>before</i> the offset/limit window, so
    /// <c>--ms-order 2 --limit 10</c> means the first ten MS2 scans rather than the MS2 scans among
    /// the first ten of any order. <c>record_count</c> counts the scans that passed the filter, and
    /// <c>scan_count</c> reports the file's true total, so a filter can never look like an empty
    /// file.
    /// </para>
    /// <para>
    /// <b>Two of the seven need native code that is not in the wheel.</b> Bruker <c>.d</c> and
    /// timsTOF <c>.d</c> reach vendor DLLs (<c>baf2sql_c.dll</c>, <c>timsdata.dll</c>) through
    /// P/Invoke and are Windows-x64 only; Thermo <c>.raw</c> uses managed vendor assemblies and
    /// works anywhere. A vendor format that cannot load surfaces as its real failure rather than as
    /// an empty read.
    /// </para>
    /// </remarks>
    public static object ReadSpectra(Program.Arguments arguments)
    {
        Window window = Window.From(arguments);
        string path = window.Path;

        RequireValueIfProvided(arguments, "ms-order");
        int? msOrder = arguments.WasProvided("ms-order") ? arguments.OptionalInt("ms-order", 0) : null;
        if (msOrder is < 1)
            throw new Program.UsageException($"Option --ms-order must be 1 or greater; got {msOrder}.");

        bool includePeaks = arguments.Flag("peaks");

        MsDataFile dataFile = OpenSpectra(path);

        // A vendor reader can hold a native handle (timsTOF is the only MsDataFile that is
        // IDisposable), so the file is disposed even when projection throws.
        try
        {
            List<MsDataScan> allScans = dataFile.GetAllScansList();
            List<MsDataScan> filtered = msOrder is null
                ? allScans
                : allScans.Where(s => s.MsnOrder == msOrder).ToList();

            IReadOnlyList<MsDataScan> selected = window.Apply(filtered, out bool truncated);

            var columns = includePeaks ? ScanColumnsWithPeaks : ScanColumns;
            object? written = null;
            if (window.OutputPath is not null)
                written = WriteTable(window.OutputPath, columns, selected);

            return new
            {
                path = Path.GetFullPath(path),
                file_type = FileTypeOf(path),
                reader = dataFile.GetType().Name,
                // The file's true scan count, always — so a --ms-order filter that matches nothing
                // is visibly a filter that matched nothing, not an empty file.
                scan_count = allScans.Count,
                ms_order = msOrder,
                record_count = filtered.Count,
                returned_count = written is null ? selected.Count : 0,
                offset = window.Offset,
                truncated,
                peaks_included = includePeaks,
                // MsDataScan.RetentionTime is minutes for every MsDataFile reader in mzLib: the
                // readers convert at the boundary rather than passing the vendor's unit through,
                // which is exactly what the result-file readers do NOT do. Stated as a value so it
                // reads the same way as every other retention time this bridge emits.
                retention_time_unit = "minutes",
                caveats = SpectraCaveatsFor(path, includePeaks),
                column_names = columns.Select(c => c.Name).ToList(),
                columns = written is null ? BuildColumns(columns, selected) : null,
                output = written,
            };
        }
        finally
        {
            (dataFile as IDisposable)?.Dispose();
        }
    }

    // ---------------------------------------------------------------------------------------
    // Opening
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// One MS1 feature, and whether its intensity is a measurement or mzLib's zero.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="ISingleChargeMs1Feature.Intensity"/> is a non-nullable double, and
    /// <c>Ms1Feature.GetSingleChargeFeatures</c> fills it with <c>IntensityApex ?? 0</c>
    /// (Ms1Feature.cs:86). <c>Apex_intensity</c> is an <c>[Optional]</c> column, and the
    /// FLASHDeconv/OpenMS <c>_ms1.feature</c> schema <b>does not have it at all</b> — so for that
    /// whole family every feature's intensity is a fabricated zero, indistinguishable from a real
    /// measurement of nothing.
    /// </para>
    /// <para>
    /// Carrying the availability alongside the feature is what lets the column cross as null
    /// instead. This is the bridge-principle case in miniature: the value is genuinely absent, the
    /// core contract cannot say so because the interface types it non-nullable, and the honest
    /// projection is optionality rather than a plausible number. It is recorded on
    /// <c>bridge/UPSTREAM.md</c> as a candidate for a nullable <c>Intensity</c> upstream, which
    /// would let this carrier go away.
    /// </para>
    /// </remarks>
    private sealed record Feature(ISingleChargeMs1Feature Value, bool IntensityMeasured);

    /// <summary>
    /// The file's features, each tagged with whether mzLib had an intensity to report.
    /// </summary>
    /// <remarks>
    /// The per-charge expansion is contiguous — <c>Ms1Feature</c> yields one feature per charge in
    /// <c>[ChargeStateMin, ChargeStateMax]</c>, in order — so a record's availability applies to a
    /// known run of outputs. If that ever stops holding, the counts disagree and this falls back to
    /// reporting every intensity as measured: a wrong "measured" is the status quo, whereas a
    /// misaligned null would mark the wrong rows.
    /// </remarks>
    private static List<Feature> FeaturesOf(IMs1FeatureFile featureFile, IResultFile resultFile)
    {
        List<ISingleChargeMs1Feature> features = featureFile.GetMs1Features().ToList();

        if (resultFile is not Ms1FeatureFile ms1File)
            return features.Select(feature => new Feature(feature, true)).ToList();

        var measured = new List<bool>(features.Count);
        foreach (Ms1Feature record in ms1File.Results)
        {
            int charges = record.ChargeStateMax - record.ChargeStateMin + 1;
            for (int i = 0; i < charges; i++)
                measured.Add(record.IntensityApex is not null);
        }

        if (measured.Count != features.Count)
            return features.Select(feature => new Feature(feature, true)).ToList();

        return features.Select((feature, index) => new Feature(feature, measured[index])).ToList();
    }

    /// <summary>Opens a path as whatever mzLib says it is, or explains why it cannot.</summary>
    private static IResultFile OpenAny(string path)
    {
        if (!File.Exists(path) && !Directory.Exists(path))
            throw new Program.UsageException($"File not found: '{path}'.");

        try
        {
            return FileReader.ReadResultFile(path);
        }
        catch (MzLibUtil.MzLibException exception)
        {
            throw new Program.UsageException(
                $"{exception.Message}: '{path}'. The formats listing enumerates every file type mzLib recognises.");
        }
    }

    /// <summary>Opens a path as a spectra file, or says what it is instead.</summary>
    private static MsDataFile OpenSpectra(string path)
    {
        if (!File.Exists(path) && !Directory.Exists(path))
            throw new Program.UsageException($"File not found: '{path}'.");

        IResultFile resultFile = OpenAny(path);
        // The seven spectra types all resolve to the same adapter, which wraps the real MsDataFile.
        // Going through MsDataFileReader rather than unwrapping the adapter keeps this on mzLib's
        // own documented entry point.
        if (resultFile is not MsDataFile && !IsSpectraType(resultFile.FileType))
            throw NoSuchView(resultFile, "spectra", "read-spectra");

        return MsDataFileReader.GetDataFile(path).LoadAllStaticData();
    }

    /// <summary>The records mzLib parsed, whatever the reader's record type is.</summary>
    /// <remarks>
    /// <see cref="IResultFile"/> does not expose the results — only <c>ResultFile&lt;T&gt;</c> does,
    /// through a property whose type differs per reader. Every reader in mzLib is either a
    /// <c>ResultFile&lt;T&gt;</c> or the MsDataFile adapter, and both expose <c>Results</c> as a
    /// list, so it is read by name. Returning an empty list rather than throwing for a shape that
    /// is neither keeps a future reader from taking the verb down.
    /// </remarks>
    private static IReadOnlyList<object> RecordsOf(IResultFile resultFile)
    {
        PropertyInfo? results = resultFile.GetType().GetProperty("Results",
            BindingFlags.Public | BindingFlags.Instance | BindingFlags.FlattenHierarchy);

        if (results?.GetValue(resultFile) is IEnumerable enumerable and not string)
            return enumerable.Cast<object>().Where(item => item is not null).ToList();

        return [];
    }

    /// <summary>The element type to project: the declared record type, or the runtime one.</summary>
    /// <remarks>
    /// The declared <c>T</c> of <c>ResultFile&lt;T&gt;</c> is preferred over the first record's own
    /// type, so an empty file still reports its column names — a table with no rows must still say
    /// what its columns would have been — and so a file whose first record happens to be a subclass
    /// (a glyco PSM among ordinary ones) does not silently give every row the subclass's columns.
    /// </remarks>
    private static Type ElementTypeOf(IResultFile resultFile, IReadOnlyList<object> records) =>
        RecordTypeOf(resultFile.GetType())
        ?? (records.Count > 0 ? records[0].GetType() : typeof(MsDataScan));

    /// <summary>Whether this reader's records implement <see cref="ISpectralMatch"/>.</summary>
    private static bool HasSpectralMatchRecords(Type readerType)
    {
        Type? recordType = RecordTypeOf(readerType);
        return recordType is not null && typeof(ISpectralMatch).IsAssignableFrom(recordType);
    }

    /// <summary>The seven file types mzLib parses into an <see cref="MsDataFile"/>.</summary>
    private static bool IsSpectraType(SupportedFileType fileType) => fileType
        is SupportedFileType.MzML
        or SupportedFileType.Mgf
        or SupportedFileType.ThermoRaw
        or SupportedFileType.BrukerD
        or SupportedFileType.BrukerTimsTof
        or SupportedFileType.Ms1Align
        or SupportedFileType.Ms2Align;

    /// <summary>The file type of a path, for a verb that did not otherwise need to open it.</summary>
    private static string FileTypeOf(string path)
    {
        try
        {
            return path.ParseFileType().ToString();
        }
        catch (MzLibUtil.MzLibException)
        {
            return "unknown";
        }
    }

    /// <summary>Rejects a file that does not offer the view a verb needs, naming what it does have.</summary>
    private static Program.UsageException NoSuchView(IResultFile resultFile, string view, string verb)
    {
        List<string> views = ViewsOf(resultFile.GetType());
        string has = views.Count == 0
            ? "it has no cross-format view at all"
            : $"it offers the {string.Join(", ", views)} view";

        return new Program.UsageException(
            $"'{resultFile.FileType}' files do not offer the {view} view, so {verb} cannot read them - " +
            $"{has}. Every file type can be read with read-records, which returns that format's own fields.");
    }

    // ---------------------------------------------------------------------------------------
    // The offset/limit/out window, shared by all four verbs
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The <c>--path</c>/<c>--offset</c>/<c>--limit</c>/<c>--out</c> quartet, parsed and validated
    /// once.
    /// </summary>
    /// <remarks>
    /// Declared once rather than repeated per verb so the four verbs cannot drift on the rules that
    /// were argued out for <see cref="Reading.ReadResults"/> — an option supplied without a value is
    /// a usage error rather than a silent default, a negative offset or limit is rejected, and
    /// <c>--out</c> may not name the input file.
    /// </remarks>
    private sealed record Window(string Path, int Offset, int Limit, string? OutputPath)
    {
        public static Window From(Program.Arguments arguments)
        {
            string path = arguments.Required("path");

            RequireValueIfProvided(arguments, "out");
            RequireValueIfProvided(arguments, "limit");
            RequireValueIfProvided(arguments, "offset");

            int offset = arguments.OptionalInt("offset", 0);
            if (offset < 0)
                throw new Program.UsageException($"Option --offset must be zero or greater; got {offset}.");

            bool limited = arguments.WasProvided("limit");
            int limit = arguments.OptionalInt("limit", int.MaxValue);
            if (limited && limit < 0)
                throw new Program.UsageException($"Option --limit must be zero or greater; got {limit}.");

            string? outputPath = arguments.Optional("out");
            if (outputPath is not null &&
                string.Equals(
                    System.IO.Path.GetFullPath(outputPath),
                    System.IO.Path.GetFullPath(path),
                    StringComparison.OrdinalIgnoreCase))
            {
                throw new Program.UsageException(
                    $"Option --out must differ from --path: writing to '{path}' would overwrite the input file.");
            }

            return new Window(path, offset, limit, string.IsNullOrWhiteSpace(outputPath) ? null : outputPath);
        }

        /// <summary>The selected window, and whether anything was left outside it.</summary>
        public IReadOnlyList<T> Apply<T>(IReadOnlyList<T> all, out bool truncated)
        {
            int start = Math.Min(Offset, all.Count);
            int count = (int)Math.Min((long)Limit, all.Count - start);

            var selected = new List<T>(count);
            for (int i = 0; i < count; i++)
                selected.Add(all[start + i]);

            // "Were any records left behind", by either the limit or the offset — deliberately not
            // `offset + count < all.Count`, which reports a complete answer for an offset past the
            // end. Same rule as ReadResults.
            truncated = selected.Count < all.Count;
            return selected;
        }
    }

    // ---------------------------------------------------------------------------------------
    // Reflective projection of an arbitrary mzLib record type
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The columns an arbitrary mzLib record type projects onto, worked out once per type.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is the one place in the bridge that reads mzLib through reflection rather than through
    /// a named interface, and it is confined to <c>read-records</c> deliberately. The alternative
    /// was to hand-transcribe the fields of nine record types carrying two hundred properties
    /// between them — which would be a second, silently-drifting copy of mzLib's schema, exactly
    /// what <see cref="Reading.Formats"/> already refuses to be for the file-type table.
    /// </para>
    /// <para>
    /// The projection is deliberately narrow about what a column may be. A scalar crosses as
    /// itself; a list of scalars crosses <c>;</c>-joined, matching how the quantifiable view already
    /// renders protein groups; <b>anything else does not cross at all</b> and is named in
    /// <see cref="Excluded"/>. Flattening a nested object into invented column names would be the
    /// bridge inventing a schema mzLib never published.
    /// </para>
    /// </remarks>
    private sealed class RecordProjection
    {
        private static readonly Dictionary<Type, RecordProjection> Cache = [];

        private readonly List<(string Name, PropertyInfo Property)> _fields;
        private readonly HashSet<string> _failed = [];
        private readonly HashSet<string> _sentinelFields;

        public string RecordTypeName { get; }
        public IReadOnlyList<string> ColumnNames { get; }
        public IReadOnlyList<object> Excluded { get; }

        private RecordProjection(Type recordType)
        {
            RecordTypeName = recordType.Name;

            var fields = new List<(string, PropertyInfo)>();
            var excluded = new List<object>();

            foreach (PropertyInfo property in PropertiesOf(recordType))
            {
                string? why = WhyNotAColumn(property);
                if (why is null)
                    fields.Add((SnakeCase(property.Name), property));
                else
                    excluded.Add(new { field = SnakeCase(property.Name), type = FriendlyTypeName(property.PropertyType), reason = why });
            }

            _fields = fields;
            ColumnNames = fields.Select(f => f.Item1).ToList();
            Excluded = excluded;
            _sentinelFields = SentinelFieldsOf(recordType);
        }

        /// <summary>
        /// The fields on this record type where mzLib documents <c>-1</c> as "absent".
        /// </summary>
        /// <remarks>
        /// <para>
        /// The general rule here is that <c>-1</c> crosses through untouched, because in a format's
        /// own columns it is usually a real measurement — a mass difference, a delta, TopPIC's
        /// <c>feature_score</c>. That rule has exactly one documented exception, and it is not
        /// optional: <see cref="IQuantifiableRecord.RetentionTime"/> and
        /// <see cref="IQuantifiableRecord.MonoisotopicMass"/> are typed as non-nullable doubles, so
        /// mzLib assigns literal <c>-1</c> when the column is missing (SpectrumMatchFromTsv.cs:198,
        /// :89).
        /// </para>
        /// <para>
        /// Those two members are projected by <c>read-results</c> as null, and
        /// <c>read-records</c> reaches <b>the same properties on the same record types</b> — so
        /// without this, one verb would answer <c>null</c> and the other <c>-1</c> for the same
        /// column of the same file, and the <c>-1</c> would enter a mean. The exception is scoped
        /// to the interface that documents it rather than to a name, so a format-specific
        /// <c>retention_time</c> elsewhere is untouched.
        /// </para>
        /// </remarks>
        private static HashSet<string> SentinelFieldsOf(Type recordType)
        {
            if (!typeof(IQuantifiableRecord).IsAssignableFrom(recordType))
                return [];

            return
            [
                nameof(IQuantifiableRecord.RetentionTime),
                nameof(IQuantifiableRecord.MonoisotopicMass),
            ];
        }

        public static RecordProjection For(Type recordType)
        {
            lock (Cache)
            {
                if (!Cache.TryGetValue(recordType, out RecordProjection? projection))
                    Cache[recordType] = projection = new RecordProjection(recordType);
                return projection;
            }
        }

        /// <summary>
        /// The public instance properties, base-class-first and in declaration order.
        /// </summary>
        /// <remarks>
        /// <c>Type.GetProperties</c> does not promise an order, and for a derived record it
        /// interleaves inherited and declared members unpredictably — which would make the column
        /// order of the same file differ between runs of the same build. Walking the base chain
        /// explicitly makes the order stable and puts the shared fields first, which is also how a
        /// reader of the mzLib source encounters them. Indexers are skipped: they take arguments
        /// and are not fields.
        /// </remarks>
        private static IEnumerable<PropertyInfo> PropertiesOf(Type recordType)
        {
            var chain = new List<Type>();
            for (Type? type = recordType; type is not null && type != typeof(object); type = type.BaseType)
                chain.Insert(0, type);

            var seen = new HashSet<string>(StringComparer.Ordinal);
            foreach (Type type in chain)
            {
                PropertyInfo[] declared = type.GetProperties(
                    BindingFlags.Public | BindingFlags.Instance | BindingFlags.DeclaredOnly);

                foreach (PropertyInfo property in declared)
                {
                    if (property.GetIndexParameters().Length > 0)
                        continue;
                    if (property.GetMethod is null || !property.CanRead)
                        continue;
                    if (seen.Add(property.Name))
                        yield return property;
                }
            }
        }

        /// <summary>Why a property cannot be a column, or null if it can.</summary>
        private static string? WhyNotAColumn(PropertyInfo property)
        {
            Type type = Nullable.GetUnderlyingType(property.PropertyType) ?? property.PropertyType;

            if (IsScalar(type))
                return null;

            if (SequenceElementOf(type) is Type element && IsScalar(Nullable.GetUnderlyingType(element) ?? element))
                return null;

            if (typeof(IDictionary).IsAssignableFrom(type))
                return "a dictionary has no faithful column shape; read it through the typed view or mzLib directly";

            if (SequenceElementOf(type) is not null)
                return "a list of composite values has no faithful column shape";

            return "a composite value; flattening it here would invent column names mzLib does not publish";
        }

        /// <summary>A value that can be a single cell.</summary>
        private static bool IsScalar(Type type) =>
            type.IsPrimitive || type.IsEnum
            || type == typeof(string) || type == typeof(decimal)
            || type == typeof(DateTime) || type == typeof(DateTimeOffset)
            || type == typeof(TimeSpan) || type == typeof(Guid);

        /// <summary>The element type of a non-string sequence, or null if it is not one.</summary>
        private static Type? SequenceElementOf(Type type)
        {
            if (type == typeof(string) || !typeof(IEnumerable).IsAssignableFrom(type))
                return null;

            if (type.IsArray)
                return type.GetElementType();

            foreach (Type contract in type.GetInterfaces().Prepend(type))
            {
                if (contract.IsGenericType && contract.GetGenericTypeDefinition() == typeof(IEnumerable<>))
                    return contract.GetGenericArguments()[0];
            }

            return typeof(object);
        }

        /// <summary>One record's value for one column, or null if it could not be read.</summary>
        /// <remarks>
        /// Each property is read under its own guard because several of mzLib's computed properties
        /// throw on real files rather than returning a default — <c>CruxResult.Accession</c> and
        /// <c>MsPathFinderTResult.Accession</c> both index into
        /// <c>ProteinId.Split('|')</c> and raise <see cref="IndexOutOfRangeException"/> on a header
        /// that is not UniProt-shaped. Letting one such property abort the read would make an
        /// entire file unreadable because of one derived field the caller may not even want; the
        /// cell becomes null and the field is named in <c>failed_fields</c>.
        /// </remarks>
        private object? Read(PropertyInfo property, object record)
        {
            object? value;
            try
            {
                value = property.GetValue(record);
            }
            catch (Exception exception)
            {
                lock (_failed)
                    _failed.Add($"{SnakeCase(property.Name)}: {(exception.InnerException ?? exception).GetType().Name}");
                return null;
            }

            if (value is double number && _sentinelFields.Contains(property.Name))
                return NullIfSentinel(number);

            return Normalize(value);
        }

        /// <summary>A read value in its wire shape.</summary>
        /// <remarks>
        /// <para>
        /// Non-finite doubles cross as null: <c>NaN</c> is mzLib's "column absent" default for the
        /// confidence fields, and JSON cannot represent it anyway — serialising it would produce
        /// invalid JSON rather than a wrong number.
        /// </para>
        /// <para>
        /// <b>The -1 sentinel is deliberately NOT applied here</b>, unlike in the quantifiable view.
        /// There, -1 is a documented "absent" marker on two specific interface members. Here the
        /// columns are each format's own fields, where -1 is frequently a real measurement — a mass
        /// difference, a delta, a log ratio — and silently nulling it would destroy data. A value
        /// that means "absent" in one format's column is that format's business to document, not
        /// this projection's to guess.
        /// </para>
        /// </remarks>
        private static object? Normalize(object? value) => value switch
        {
            null => null,
            // Delegated so the reflective projection and the typed views cannot end up with two
            // different answers for the same non-finite double.
            double or float => WireValue(value),
            DateTime moment => moment.ToString("O", CultureInfo.InvariantCulture),
            DateTimeOffset moment => moment.ToString("O", CultureInfo.InvariantCulture),
            TimeSpan span => span.ToString(null, CultureInfo.InvariantCulture),
            Enum name => name.ToString(),
            string text => text,
            IEnumerable sequence => string.Join(";", sequence.Cast<object?>()
                .Select(item => Normalize(item)?.ToString() ?? string.Empty)),
            _ => value,
        };

        public Dictionary<string, List<object?>> BuildColumns(IReadOnlyList<object> records)
        {
            var built = new Dictionary<string, List<object?>>(_fields.Count);
            foreach ((string name, PropertyInfo property) in _fields)
            {
                var values = new List<object?>(records.Count);
                foreach (object record in records)
                    values.Add(Read(property, record));
                built[name] = values;
            }

            return built;
        }

        /// <summary>The cell values of one record, in column order.</summary>
        public IEnumerable<object?> Cells(object record) =>
            _fields.Select(field => Read(field.Property, record));

        /// <summary>
        /// The fields that threw while reading the given records, with the exception type.
        /// </summary>
        /// <remarks>
        /// Reads the records first so the answer describes this call rather than a previous one:
        /// the projection is cached per type across verbs within one process, and a stale failure
        /// list would attribute another file's bad rows to this one.
        /// </remarks>
        public IReadOnlyList<string> FailedFieldsFor(IReadOnlyList<object> records)
        {
            lock (_failed)
                _failed.Clear();

            foreach ((_, PropertyInfo property) in _fields)
                foreach (object record in records)
                    Read(property, record);

            lock (_failed)
                return _failed.Order(StringComparer.Ordinal).ToList();
        }

        /// <summary>A type name a non-.NET caller can read.</summary>
        private static string FriendlyTypeName(Type type)
        {
            Type? underlying = Nullable.GetUnderlyingType(type);
            if (underlying is not null)
                return FriendlyTypeName(underlying) + "?";

            if (type.IsArray)
                return FriendlyTypeName(type.GetElementType()!) + "[]";

            if (!type.IsGenericType)
                return type.Name;

            string name = type.Name[..type.Name.IndexOf('`')];
            return $"{name}<{string.Join(", ", type.GetGenericArguments().Select(FriendlyTypeName))}>";
        }
    }

    /// <summary>mzLib's PascalCase property names as the wire's snake_case column names.</summary>
    /// <remarks>
    /// The envelope is snake_case throughout, so a column called <c>PrecursorCharge</c> would be
    /// the only camel-cased thing a caller ever sees. Consecutive capitals are kept together, so
    /// <c>EValue</c> becomes <c>e_value</c> and <c>MIScore</c> becomes <c>mi_score</c> rather than
    /// <c>m_i_score</c>; a digit starts no new word, so <c>MS2RetentionTime</c> becomes
    /// <c>ms2_retention_time</c>; and a pluralising <c>s</c> belongs to the acronym it follows, so
    /// <c>FixedPTMs</c> becomes <c>fixed_ptms</c> rather than <c>fixed_pt_ms</c>.
    /// </remarks>
    private static string SnakeCase(string name)
    {
        var built = new System.Text.StringBuilder(name.Length + 8);

        for (int i = 0; i < name.Length; i++)
        {
            char current = name[i];
            if (char.IsUpper(current) && i > 0)
            {
                bool previousWasLowerOrDigit = char.IsLower(name[i - 1]) || char.IsDigit(name[i - 1]);
                if (previousWasLowerOrDigit || (char.IsUpper(name[i - 1]) && StartsNewWord(name, i)))
                    built.Append('_');
            }

            built.Append(char.ToLowerInvariant(current));
        }

        return built.ToString();
    }

    /// <summary>
    /// Whether the capital at <paramref name="index"/>, inside a run of capitals, begins a new word.
    /// </summary>
    /// <remarks>
    /// The ordinary rule is "a capital followed by a lowercase letter begins a word" — the <c>S</c>
    /// of <c>MIScore</c>. It has one exception, and mzLib hits it: a lone trailing <c>s</c> is
    /// pluralising the acronym rather than opening a word, so <c>FixedPTMs</c> must not split as
    /// <c>fixed_pt_ms</c>. A following <c>s</c> only counts as a new word when a word actually
    /// follows it.
    /// </remarks>
    private static bool StartsNewWord(string name, int index)
    {
        if (index + 1 >= name.Length || !char.IsLower(name[index + 1]))
            return false;

        bool pluralising = name[index + 1] == 's'
            && (index + 2 >= name.Length || char.IsUpper(name[index + 2]));

        return !pluralising;
    }

    /// <summary>Writes a reflectively-projected table, in the same shape as <see cref="WriteTable"/>.</summary>
    private static object WriteProjectedTable(
        string outputPath, RecordProjection projection, IReadOnlyList<object> records)
    {
        string? directory = Path.GetDirectoryName(Path.GetFullPath(outputPath));
        if (!string.IsNullOrEmpty(directory))
            Directory.CreateDirectory(directory);

        var configuration = new CsvHelper.Configuration.CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
        };

        using (var writer = new StreamWriter(File.Create(outputPath)))
        using (var csv = new CsvHelper.CsvWriter(writer, configuration))
        {
            foreach (string name in projection.ColumnNames)
                csv.WriteField(name);
            csv.NextRecord();

            foreach (object record in records)
            {
                foreach (object? cell in projection.Cells(record))
                    csv.WriteField(Render(cell));
                csv.NextRecord();
            }
        }

        return new
        {
            path = Path.GetFullPath(outputPath),
            format = "tsv",
            row_count = records.Count,
        };
    }
}
