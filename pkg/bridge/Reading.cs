using System.Globalization;
using MassSpectrometry;
using Readers;

namespace MzLibBridge;

/// <summary>
/// Reading proteomics result files: what a file <i>is</i>, and what can be done with it.
/// </summary>
/// <remarks>
/// <para>
/// mzLib's <c>Readers</c> recognises 29 file types written by a dozen different search and
/// deconvolution tools, and dispatches each to a parser it maintains. That dispatch is the whole
/// value here — the bridge adds no parsing of its own; it asks mzLib what a path is
/// (<see cref="SupportedFileTypeExtensions.ParseFileType"/>) and reports the answer.
/// </para>
/// <para>
/// The important honesty in this verb is <c>views</c>. It is tempting to describe mzLib as reading
/// 29 formats into one uniform shape; it does not. Only three types implement
/// <see cref="IQuantifiableResultFile"/> (MetaMorpheus <c>.psmtsv</c> and <c>.osmtsv</c>, and
/// MSFragger <c>psm.tsv</c>); deconvolution features have their own interface
/// (<see cref="IMs1FeatureFile"/>); spectra files are an <see cref="MsDataFile"/>; and several
/// formats belong to no common interface at all. Rather than discovering that as a cast failure
/// deep in a later call, a caller is told up front which views a given file actually supports.
/// </para>
/// </remarks>
internal static class Reading
{
    /// <summary>
    /// <c>readers formats</c> — every file type mzLib can recognise, with its extension and the
    /// views it supports.
    /// </summary>
    /// <remarks>
    /// Enumerated from <see cref="SupportedFileType"/> rather than transcribed, so the published
    /// table of supported formats cannot drift from what mzLib actually dispatches. This mirrors
    /// how mzLib's own test suite enumerates the enum instead of maintaining a list.
    /// </remarks>
    public static object Formats(Program.Arguments arguments)
    {
        var formats = new List<object>();

        foreach (SupportedFileType fileType in Enum.GetValues<SupportedFileType>())
        {
            // A member missing its extension or reader mapping is a broken mzLib build, not a
            // reason to fail the whole listing: report what is known and keep going, so this verb
            // stays usable for diagnosing exactly that.
            string? extension = Try(() => fileType.GetFileExtension());
            Type? readerType = Try(() => fileType.GetResultFileType());

            formats.Add(new
            {
                file_type = fileType.ToString(),
                extension,
                reader = readerType?.Name,
                views = readerType is null ? new List<string>() : ViewsOf(readerType),
            });
        }

        return new
        {
            format_count = formats.Count,
            formats,
        };
    }

    /// <summary>
    /// <c>readers identify --path FILE</c> — what kind of result file this is, and what can be
    /// done with it, without parsing its contents.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Deliberately cheap: <see cref="FileReader.ReadResultFile"/> is lazy — it resolves the type
    /// and sets the path but does not call <c>LoadResults()</c> — so identifying a million-row file
    /// costs no more than identifying an empty one.
    /// </para>
    /// <para>
    /// "Without parsing" is not quite "without reading". mzLib's own type detection opens the file
    /// for three families: a bare <c>.tsv</c> is disambiguated by its first line, a <c>.mztab</c>
    /// by its first five, and a <c>.d</c> by which analysis file the directory contains. So this
    /// verb can fail on an unreadable file, which is why the path is checked first.
    /// </para>
    /// <para>
    /// <b>No <c>software</c> field, deliberately.</b> <see cref="IResultFile.Software"/> looks like
    /// the answer to "which tool wrote this" and is not. Every reader carries its software constant
    /// on its <i>path-taking</i> constructor (<c>MsFraggerPsmFile(string) : base(filePath,
    /// Software.MsFragger)</c>), but <see cref="FileReader"/> builds readers through the
    /// parameterless one, which sets <see cref="Software.Unspecified"/> — so the property is
    /// <c>Unspecified</c> for everything mzLib's own factory returns. Nor is it reliable through
    /// the other constructor: <c>MsFraggerPeptideFile</c> passes <c>Unspecified</c> outright and the
    /// psmtsv readers route to a base overload that never sets it. Reconstructing the value here —
    /// by probing the other constructor, or by mapping file type to tool — would mean the bridge
    /// answering a question mzLib cannot, which is precisely the drift the translation-layer rule
    /// exists to prevent. <c>file_type</c> and <c>reader</c> already name the tool
    /// (<c>MsFraggerPsm</c>, <c>ToppicPrsm</c>), so nothing is lost but a false certainty.
    /// </para>
    /// </remarks>
    public static object Identify(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        // Checked here rather than left to mzLib because FileReader throws a bare
        // FileNotFoundException carrying neither a message nor a file name, which would reach the
        // caller as an empty error. A missing path is the caller's mistake, so it is a usage
        // failure (exit 2), not a runtime one.
        // Directory.Exists as well as File.Exists: Bruker's .d "file" is a directory.
        if (!File.Exists(path) && !Directory.Exists(path))
            throw new Program.UsageException($"File not found: '{path}'.");

        IResultFile resultFile;
        try
        {
            resultFile = FileReader.ReadResultFile(path);
        }
        catch (MzLibUtil.MzLibException exception)
        {
            // mzLib recognises the file's extension or it does not; there is no sentinel and no
            // "unknown" enum member to return. An unsupported file is a bad argument rather than a
            // fault, so it surfaces as a usage error pointing at the verb that lists what IS
            // supported — not as a BridgeError, which would read as "pyMzLib is broken".
            throw new Program.UsageException(
                $"{exception.Message}: '{path}'. The formats listing enumerates every file type mzLib recognises.");
        }

        Type readerType = resultFile.GetType();

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            extension = Try(() => resultFile.FileType.GetFileExtension()),
            reader = readerType.Name,
            views = ViewsOf(readerType),
        };
    }

    /// <summary>
    /// <c>readers read-results --path FILE [--limit N] [--offset N] [--out FILE]</c> — the
    /// cross-format record view of a result file.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Only the three file types offering the <c>quantifiable</c> view can be read this way; a file
    /// without it is rejected with a message naming the views it does have, rather than a cast
    /// failure. Use <c>readers identify</c> first to find out.
    /// </para>
    /// <para>
    /// Records come back <b>columnar</b> — one array per field rather than one object per record.
    /// That is three to five times smaller, but the reason is ergonomic: pyMzLib takes no
    /// third-party dependency, so it can never hand back a DataFrame, and a map of arrays is the
    /// one shape both pandas and polars ingest in a single call. The same shape is equally natural
    /// in any other language, so it costs the language-neutral contract nothing.
    /// </para>
    /// <para>
    /// <b>There is no default row cap.</b> A result file can carry a million rows, and truncating
    /// by default would mean the common call returns a table that looks complete and is not — the
    /// failure mode that quietly produces a wrong answer. <c>--limit</c> is available and reports
    /// <c>truncated</c> when it bites; for genuinely large files <c>--out</c> writes the table to
    /// disk and returns only a summary, which is the intended path rather than an escape hatch.
    /// </para>
    /// </remarks>
    public static object ReadResults(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        // An option written without a value lands in the flag set, not the named set, so Optional()
        // returns null and the option is silently discarded: '--out' with no path would skip the
        // write and serialise the whole table inline — the exact large-payload case --out exists to
        // avoid — and '--limit' would degenerate to no limit while the caller believes they asked
        // for a subset. This is the rule the PRIDE download verb already states: an option that was
        // ASKED FOR but degenerates to nothing must fail, never silently widen.
        RequireValueIfProvided(arguments, "out");
        RequireValueIfProvided(arguments, "limit");
        RequireValueIfProvided(arguments, "offset");

        string? outputPath = arguments.Optional("out");

        int offset = arguments.OptionalInt("offset", 0);
        if (offset < 0)
            throw new Program.UsageException($"Option --offset must be zero or greater; got {offset}.");

        bool limited = arguments.WasProvided("limit");
        int limit = arguments.OptionalInt("limit", int.MaxValue);
        if (limited && limit < 0)
            throw new Program.UsageException($"Option --limit must be zero or greater; got {limit}.");

        if (!File.Exists(path) && !Directory.Exists(path))
            throw new Program.UsageException($"File not found: '{path}'.");

        IQuantifiableResultFile resultFile = OpenQuantifiable(path);

        // Touching the results materialises the whole file: every reader's LoadResults ends in
        // ToList(), and GetQuantifiableResults() is `=> Results` on all of them. It LOOKS lazy and
        // is not, so --offset is a window over an already-parsed list, never a cursor that saves
        // work. Paging a large file re-reads and re-parses it once per page; --out exists so that
        // is never the right thing to do.
        List<IQuantifiableRecord> all = resultFile.GetQuantifiableResults().ToList();

        // GetRange over Skip().Take().ToList(): the records are already a materialised List, so the
        // LINQ form allocates a second full copy of what can be a million-row table purely to take
        // a window of it. GetRange copies only the window.
        int start = Math.Min(offset, all.Count);
        int count = (int)Math.Min((long)limit, all.Count - start);
        List<IQuantifiableRecord> selected = all.GetRange(start, count);

        // "Were any records left behind", by either the limit or the offset. Deliberately not
        // `offset + selected.Count < all.Count`, which reads plausibly and is wrong: an offset past
        // the end makes the sum exceed the total and reports a complete answer for an empty one.
        bool truncated = selected.Count < all.Count;

        var columns = Column.QuantifiableView;
        object? written = null;
        if (!string.IsNullOrWhiteSpace(outputPath))
            written = WriteTable(outputPath, columns, selected);

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            record_count = all.Count,
            returned_count = written is null ? selected.Count : 0,
            offset,
            // True whenever records were left behind, whether by --limit or by --offset. A short
            // answer and a complete one must never look alike.
            truncated,
            // mzLib's psmtsv reader catches a malformed line, adds it to a warnings list, and the
            // ResultFile wrapper discards that list — so a half-corrupt file reads "successfully"
            // with silently fewer rows. mzLib exposes no way to ask, so the rows are counted here
            // and the difference reported. Null when the count is not meaningful for this input.
            rows_not_read = UnreadRowCount(path, all.Count, resultFile.FileType),
            // What the uniform view cannot be trusted to mean for THIS format. See CaveatsFor.
            caveats = CaveatsFor(resultFile.FileType),
            // The unit of retention_time, as a value rather than as prose. The caveats say the same
            // thing in English, but a caller converting between formats should not have to grep a
            // sentence for the word "SECONDS" - which is exactly what a reader did before this
            // existed. "unknown" when mzLib gives no basis to claim one, never a guess.
            retention_time_unit = RetentionTimeUnitOf(resultFile.FileType),
            column_names = columns.Select(c => c.Name).ToList(),
            // Omitted entirely when writing to disk: materialising both would defeat the point.
            columns = written is null ? BuildColumns(columns, selected) : null,
            output = written,
        };
    }

    /// <summary>
    /// Rejects an option that was supplied with no value, rather than treating it as absent.
    /// </summary>
    private static void RequireValueIfProvided(Program.Arguments arguments, string name)
    {
        if (arguments.WasProvided(name) && string.IsNullOrWhiteSpace(arguments.Optional(name)))
            throw new Program.UsageException(
                $"Option --{name} was given but has no value; omit it to use the default.");
    }

    /// <summary>
    /// Opens a file for the quantifiable view, or explains precisely why it cannot be.
    /// </summary>
    /// <remarks>
    /// mzLib throws the same <c>MzLibException</c> for "this extension is unknown" and for "this
    /// type exists but implements the wrong interface", distinguishable only by the message text.
    /// The second case is the interesting one and deserves a real answer, so it is re-derived
    /// through <see cref="FileReader.ReadResultFile"/> and reported with the views the file
    /// actually has.
    /// </remarks>
    private static IQuantifiableResultFile OpenQuantifiable(string path)
    {
        try
        {
            return FileReader.ReadQuantifiableResultFile(path);
        }
        catch (FileNotFoundException)
        {
            // ReadQuantifiableResultFile checks File.Exists only, so a Bruker .d directory lands
            // here despite existing. Either way the caller gets a path they can act on.
            throw new Program.UsageException(
                $"File not found, or not a readable result file: '{path}'.");
        }
        catch (MzLibUtil.MzLibException)
        {
            string detail;
            try
            {
                IResultFile any = FileReader.ReadResultFile(path);
                List<string> views = ViewsOf(any.GetType());
                detail = views.Count == 0
                    ? $"'{any.FileType}' files have no cross-format record view at all — mzLib parses " +
                      "them into a format-specific shape only."
                    : $"'{any.FileType}' files offer the {string.Join(", ", views)} view, not quantifiable.";
            }
            catch (MzLibUtil.MzLibException)
            {
                // Only an mzLib dispatch failure means "unrecognised". Catching everything here
                // would report an I/O error, or a bug in the view projection, as a bad file type
                // and send the caller looking in the wrong place.
                detail = "mzLib does not recognise this file type.";
            }
            catch (FileNotFoundException)
            {
                detail = "mzLib does not recognise this file type.";
            }

            throw new Program.UsageException(
                $"Cannot read '{path}' into the uniform record view. {detail} " +
                "Identifying a file reports the views it supports; the formats listing names the " +
                "three types that offer this one.");
        }
    }

    /// <summary>
    /// One field of the uniform view: its wire name and how to read it off a record.
    /// </summary>
    /// <remarks>
    /// Declared once so the JSON columns and the written table cannot disagree about which fields
    /// exist, what they are called, or what order they are in.
    /// </remarks>
    private sealed record Column(string Name, Func<IQuantifiableRecord, object?> Read)
    {
        /// <summary>
        /// The <see cref="IQuantifiableRecord"/> fields, under mzLib's own names.
        /// </summary>
        /// <remarks>
        /// <c>retention_time</c> and <c>monoisotopic_mass</c> use <c>-1</c> as a "not present"
        /// sentinel in mzLib because the interface types them as non-nullable doubles. Passing that
        /// through would put a real-looking -1 into someone's arithmetic, so it crosses as null.
        /// The protein tuple list is flattened into three parallel <c>;</c>-joined fields, matching
        /// how the FlashLFQ tranche already renders protein groups.
        /// </remarks>
        public static IReadOnlyList<Column> QuantifiableView { get; } = new[]
        {
            new Column("file_name", r => r.FileName),
            new Column("base_sequence", r => r.BaseSequence),
            new Column("full_sequence", r => r.FullSequence),
            new Column("retention_time", r => NullIfSentinel(r.RetentionTime)),
            new Column("charge_state", r => r.ChargeState),
            new Column("monoisotopic_mass", r => NullIfSentinel(r.MonoisotopicMass)),
            // Null, not false, where the format cannot report decoys. mzLib hardcodes false for
            // MSFragger because psm.tsv carries no target/decoy column at all - so the value means
            // "unknown", and a boolean column that silently means "unknown" for part of a table is
            // the single most dangerous thing this view could hand back. Both bake-off rounds named
            // it. Reported upstream; until the contract can express it, the wire refuses to.
            new Column("is_decoy", r => DecoysAreKnown(r) ? r.IsDecoy : (bool?)null),
            // Named exactly as mzLib names the tuple fields —
            // List<(string proteinAccessions, string geneName, string organism)> — including the
            // singulars, even though each carries a ';'-joined list here. A caller reading the
            // mzLib source must not have to hold a translation table.
            new Column("protein_accessions", r => Join(r, p => p.proteinAccessions)),
            new Column("gene_name", r => Join(r, p => p.geneName)),
            new Column("organism", r => Join(r, p => p.organism)),
        };

        /// <summary>Whether this record's format can report decoy status at all.</summary>
        /// <remarks>
        /// <c>MsFraggerPsm.IsDecoy</c> is <c>=> false</c> with a comment saying decoy reading is
        /// unsupported, and the file genuinely has no target/decoy column - FragPipe strips decoys
        /// before writing it. So false is not an answer, and passing it through would let a caller
        /// group by a column that is fabricated for one format and real for another.
        /// </remarks>
        private static bool DecoysAreKnown(IQuantifiableRecord record) =>
            // An ALLOWLIST, not a denylist. Naming the formats that cannot report decoys would make
            // "decoys are known" the default for any type mzLib adds later — so a new reader that
            // also hardcodes false would start emitting fabricated booleans without anyone editing
            // this file. Defaulting to "unknown" is wrong in the harmless direction.
            record is SpectrumMatchFromTsv or LightWeightSpectralMatch;

        private static string Join(
            IQuantifiableRecord record, Func<(string proteinAccessions, string geneName, string organism), string> part)
            => string.Join(";", (record.ProteinGroupInfos ?? []).Select(part));

        /// <summary>mzLib's -1 "absent" sentinel, and any non-finite value, as null.</summary>
        private static double? NullIfSentinel(double value) =>
            double.IsFinite(value) && value != -1 ? value : null;
    }

    /// <summary>Builds the columnar payload: one array per field.</summary>
    private static Dictionary<string, List<object?>> BuildColumns(
        IReadOnlyList<Column> columns, List<IQuantifiableRecord> records)
    {
        var built = new Dictionary<string, List<object?>>(columns.Count);
        foreach (Column column in columns)
        {
            var values = new List<object?>(records.Count);
            foreach (IQuantifiableRecord record in records)
                values.Add(column.Read(record));
            built[column.Name] = values;
        }

        return built;
    }

    /// <summary>
    /// Writes the selected records as a tab-separated table and reports where they went.
    /// </summary>
    /// <remarks>
    /// Tab-separated, not comma-separated, because these fields contain commas: MSFragger's mapped
    /// proteins are a comma-separated list inside a single field, and joined accessions look the
    /// same. Tabs do not occur in them — if one did, mzLib's own tab-splitting readers would
    /// already be broken — so the delimiter is safe by inheritance rather than by hope. It is also
    /// what mzLib writes everywhere (every Delimiter in Readers is a tab bar one visualization
    /// format) and what the FlashLFQ verb already emits. Fields are still quoted when they would
    /// otherwise contain a delimiter or newline, so the file is lossless rather than lucky.
    /// </remarks>
    private static object WriteTable(
        string outputPath, IReadOnlyList<Column> columns, List<IQuantifiableRecord> records)
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
            foreach (Column column in columns)
                csv.WriteField(column.Name);
            csv.NextRecord();

            foreach (IQuantifiableRecord record in records)
            {
                foreach (Column column in columns)
                    csv.WriteField(Render(column.Read(record)));
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

    /// <summary>A value as written text: invariant, and empty for absent.</summary>
    private static string Render(object? value) => value switch
    {
        null => string.Empty,
        double number => number.ToString(CultureInfo.InvariantCulture),
        int number => number.ToString(CultureInfo.InvariantCulture),
        bool flag => flag ? "true" : "false",
        _ => value.ToString() ?? string.Empty,
    };

    /// <summary>
    /// How many data rows the file appears to hold that did not become records.
    /// </summary>
    /// <remarks>
    /// mzLib's psmtsv reader collects per-line parse failures into a warnings list that the
    /// <c>ResultFile</c> wrapper then discards, and the CsvHelper-backed readers configured with
    /// <c>BadDataFound = null</c> drop malformed rows without a sound. Either way the caller sees a
    /// successful read with fewer rows than the file contains, which is the kind of silent loss
    /// this library exists to refuse. Since mzLib offers no way to ask, the lines are counted.
    /// <para>
    /// Deliberately conservative: null rather than a guess whenever the count could mislead — for a
    /// directory-shaped input, an unreadable file, or a negative difference (which would mean the
    /// format emits more records than lines, as an expand-per-charge reader would).
    /// </para>
    /// </remarks>
    private static int? UnreadRowCount(string path, int recordCount, SupportedFileType fileType)
    {
        // Only where "one non-blank line after the header is one record" actually holds. It does
        // for the psmtsv family and MSFragger, and it is worth stating rather than relying on: a
        // TopPIC prsm file has a 29-line uncommented preamble AND continuation rows that repeat a
        // match's protein fields with every other column blank, so eight lines are four records.
        // If mzLib ever gives such a format the quantifiable view, a line count would silently
        // report half the file as unread. Null means "no basis to say", not "nothing missing".
        bool oneLinePerRecord = fileType is SupportedFileType.psmtsv
            or SupportedFileType.osmtsv
            or SupportedFileType.MsFraggerPsm;
        if (!oneLinePerRecord)
            return null;

        if (!File.Exists(path))
            return null;

        try
        {
            int dataLines = File.ReadLines(path).Count(line => !string.IsNullOrWhiteSpace(line)) - 1;
            int difference = dataLines - recordCount;
            return difference > 0 ? difference : difference == 0 ? 0 : null;
        }
        catch (IOException)
        {
            return null;
        }
    }

    /// <summary>
    /// The unit mzLib's <c>RetentionTime</c> carries for a given format.
    /// </summary>
    /// <remarks>
    /// mzLib's result-file readers do not normalise retention time, so the unit is whatever the
    /// tool wrote and differs per format. Reported as a value so a caller can convert
    /// programmatically; the alternative is every consumer hard-coding this table, which is what
    /// happened when only the prose caveat existed.
    /// <para>
    /// Verified per format against the committed fixtures, and pinned by a test. If mzLib starts
    /// normalising (see the upstream retention-time issue), these values change with it - which is
    /// the point of deriving them from the pinned mzLib rather than documenting them once.
    /// </para>
    /// </remarks>
    private static string RetentionTimeUnitOf(SupportedFileType fileType) => fileType switch
    {
        // MetaMorpheus writes minutes; verified against BottomUpExample.psmtsv, where scan 13955
        // reads 97.42 and scan 27567 reads 174.96 - 0.0057 per scan, a normal Orbitrap duty cycle.
        SupportedFileType.psmtsv or SupportedFileType.osmtsv => "minutes",
        // MSFragger writes seconds; the fixture advances ~1.1 per scan on an LTQ Orbitrap Velos.
        SupportedFileType.MsFraggerPsm => "seconds",
        _ => "unknown",
    };

    /// <summary>
    /// What the uniform view cannot be trusted to mean for a given format.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is a deliberate and narrow exception to the rule that the bridge holds no domain
    /// knowledge. mzLib's readers pass each tool's columns through without normalising them, so the
    /// same field means different things per format — and mzLib exposes no programmatic way to ask
    /// which. The alternatives were to say nothing (leaving callers to compare seconds against
    /// minutes) or to convert (inventing numbers mzLib never produced). Reporting the discrepancy
    /// is the only option that neither hides nor fabricates.
    /// </para>
    /// <para>
    /// Every entry is verified against mzLib at the pinned commit and pinned by a test, so a fix
    /// upstream shows up as a failure here rather than as a caveat that quietly became a lie.
    /// </para>
    /// </remarks>
    private static List<string> CaveatsFor(SupportedFileType fileType) => fileType switch
    {
        SupportedFileType.MsFraggerPsm =>
        [
            "retention_time is in SECONDS for this format, not minutes: MSFragger's Retention " +
            "column is passed through unconverted (MsFraggerPsm.cs:48). Do not compare it with a " +
            "MetaMorpheus file's retention_time, and do not quantify this file with FlashLFQ, " +
            "which reads the value as minutes.",
            "is_decoy is always false: mzLib does not read MSFragger decoys (MsFraggerPsm.cs:217). " +
            "False means 'unknown', not 'target'.",
            // Corrected after the readers bake-off: an earlier version of this caveat claimed the
            // psmtsv formats report the OBSERVED mass, which is false - they report the file's
            // "Peptide Monoisotopic Mass" column, which is theoretical, exactly as MSFragger does.
            // The caveat manufactured a cross-format discrepancy that does not exist and sent a
            // reader chasing it. Both formats agree; what is worth saying is only that neither is
            // the observed precursor mass.
            "monoisotopic_mass is the THEORETICAL peptide mass (MsFraggerPsm.cs:220, " +
            "CalculatedPeptideMass), not the observed precursor mass. The psmtsv formats report " +
            "the theoretical mass here too, so the two are consistent - but neither is what the " +
            "instrument measured.",
            "file_name is the full 'Spectrum File' path including its .pep.xml extension, whereas " +
            "the psmtsv formats report a bare base name. The field is not a join key across formats.",
        ],
        SupportedFileType.psmtsv or SupportedFileType.osmtsv =>
        [
            "full_sequence and monoisotopic_mass keep only the FIRST candidate of an ambiguous " +
            "identification; mzLib splits the '|'-separated list and discards the rest " +
            "(SpectrumMatchFromTsv.cs:89).",
            "monoisotopic_mass is the file's 'Peptide Monoisotopic Mass' - the THEORETICAL mass of " +
            "the identified peptide, not the observed precursor mass, which the file carries " +
            "separately as 'Precursor Mass'.",
            "There is no q-value, PEP or score in this view, so nothing here is FDR-filtered. " +
            "IQuantifiableRecord carries only the fields FlashLFQ needs; confidence columns present " +
            "in the file are not exposed. Filter before you report.",
            "A malformed row is dropped silently: mzLib collects a warning per unreadable line and " +
            "the reader discards the list (SpectrumMatchTsvReader.cs:71, PsmFromTsvFile.cs:17). " +
            "Check rows_not_read.",
        ],
        _ => [],
    };

    /// <summary>
    /// The common interfaces a reader implements, named as capabilities rather than as .NET types.
    /// </summary>
    /// <remarks>
    /// An empty list is a real and common answer - TopPIC, MsPathFinderT, Crux and Casanovo readers
    /// each parse into their own record type and share no cross-format interface. Saying so is the
    /// point of the verb.
    /// </remarks>
    private static List<string> ViewsOf(Type readerType)
    {
        var views = new List<string>();

        if (typeof(IQuantifiableResultFile).IsAssignableFrom(readerType))
            views.Add("quantifiable");
        if (typeof(IMs1FeatureFile).IsAssignableFrom(readerType))
            views.Add("ms1_features");
        if (typeof(MsDataFile).IsAssignableFrom(readerType))
            views.Add("spectra");

        // ISpectralMatch is implemented by the RECORD, not the file, so it cannot be found on the
        // reader itself - it is the type argument of the ResultFile<T> the reader derives from.
        Type? recordType = RecordTypeOf(readerType);
        if (recordType is not null && typeof(ISpectralMatch).IsAssignableFrom(recordType))
            views.Add("spectral_match");

        return views;
    }

    /// <summary>The <c>T</c> of the <c>ResultFile&lt;T&gt;</c> a reader derives from, if it does.</summary>
    private static Type? RecordTypeOf(Type readerType)
    {
        for (Type? type = readerType; type is not null; type = type.BaseType)
        {
            if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(ResultFile<>))
                return type.GetGenericArguments()[0];
        }

        return null;
    }

    /// <summary>
    /// Evaluates a lookup that throws rather than returning null for an unmapped value.
    /// </summary>
    /// <remarks>
    /// Used only by <see cref="Formats"/>, where one unmapped enum member must not take the whole
    /// listing down with it. Not a general-purpose swallow: every other path here lets failures
    /// propagate so the caller gets the real exception type.
    /// </remarks>
    private static T? Try<T>(Func<T> lookup) where T : class
    {
        try
        {
            return lookup();
        }
        catch (MzLibUtil.MzLibException)
        {
            return null;
        }
    }
}
