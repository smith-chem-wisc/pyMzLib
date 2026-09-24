using System.Globalization;
using MassSpectrometry;
using Readers;
using Readers.ExternalResults.IndividualResultRecords;
using Readers.ExternalResults.ResultFiles;

namespace MzLibBridge;

/// <summary>
/// The column sets and per-format caveats of the three typed views added alongside
/// <c>quantifiable</c>: <c>ms1_features</c>, <c>spectral_match</c>, and <c>spectra</c>.
/// </summary>
/// <remarks>
/// Each view's columns are mzLib's own interface members under mzLib's own names, exactly as the
/// quantifiable view is. The caveats follow the same rule the quantifiable ones do: they state only
/// what is verified against mzLib at the pinned commit and pinned by a test, so an upstream fix
/// surfaces here as a failing test rather than as a caveat that quietly became a lie.
/// </remarks>
internal static partial class Reading
{
    // ---------------------------------------------------------------------------------------
    // ms1_features
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The <see cref="ISingleChargeMs1Feature"/> fields, under mzLib's own names.
    /// </summary>
    private static IReadOnlyList<Column<Feature>> FeatureColumns { get; } = new[]
    {
        new Column<Feature>("mz", f => f.Value.Mz),
        new Column<Feature>("charge", f => f.Value.Charge),
        new Column<Feature>("retention_time_start", f => f.Value.RetentionTimeStart),
        new Column<Feature>("retention_time_end", f => f.Value.RetentionTimeEnd),
        // Null where mzLib had no apex intensity to report. The interface types Intensity as a
        // non-nullable double and Ms1Feature fills it with `IntensityApex ?? 0`, so for a file
        // whose schema lacks the optional Apex_intensity column — every FLASHDeconv/OpenMS
        // _ms1.feature — a plain projection hands back a whole column of fabricated zeros that
        // look exactly like measurements of nothing. Optionality is the honest projection.
        new Column<Feature>("intensity", f => f.IntensityMeasured ? f.Value.Intensity : (double?)null),
        // Genuinely nullable on the interface, and null for a whole format rather than for odd
        // rows: mzLib's _ms1.feature expansion never sets it. Crossing as null is the faithful
        // projection; a zero would read as "no isotopes were found".
        new Column<Feature>("number_of_isotopes", f => f.Value.NumberOfIsotopes),
    };

    /// <summary>The unit the feature view's retention times carry, for a given file.</summary>
    /// <remarks>
    /// Deliberately <c>unknown</c> for <c>_ms1.feature</c>. TopFD wrote seconds through v1.6.2 and
    /// minutes from v1.7.0 — <b>within the same file type</b> — and mzLib does not normalise
    /// either. mzLib's own deconvolution parameters guess with a heuristic (if the largest end time
    /// exceeds 500, divide everything by 60), which is a guess this bridge will not launder into a
    /// stated fact. Dinosaur writes minutes and is reported as such.
    /// </remarks>
    private static string FeatureRetentionTimeUnitOf(IResultFile resultFile) => resultFile.FileType switch
    {
        SupportedFileType.Tsv_Dinosaur => "minutes",
        _ => "unknown",
    };

    /// <summary>What the feature view cannot be trusted to mean, per format.</summary>
    /// <remarks>
    /// Per-format, not shared. The two implementers of this view differ in a way that matters:
    /// <c>Ms1FeatureFile</c> <b>expands</b> each row across its charge range, while
    /// <c>DinosaurTsvFile.GetMs1Features</c> is <c>=&gt; Results</c> and returns the file's rows
    /// one for one. A single "rows are expanded" caveat would therefore be false for Dinosaur —
    /// the same class of manufactured discrepancy the readers bake-off already caught once in the
    /// quantifiable caveats, and worth not repeating.
    /// </remarks>
    private static List<string> FeatureCaveatsFor(IResultFile resultFile, IReadOnlyList<Feature> features)
    {
        List<string> caveats = BaseFeatureCaveatsFor(resultFile.FileType);

        // Said only when it is true of THIS file, because it depends on the schema the writer used
        // rather than on the file type: TopFD writes Apex_intensity and FLASHDeconv does not, and
        // both are SupportedFileType.Ms1Feature.
        int unmeasured = features.Count(feature => !feature.IntensityMeasured);
        if (unmeasured == features.Count && features.Count > 0)
        {
            caveats.Add(
                "intensity is NULL for every row of this file. mzLib takes the per-charge intensity " +
                "from the optional Apex_intensity column (Ms1Feature.cs:86) and this file's schema " +
                "does not have it — the FLASHDeconv/OpenMS _ms1.feature layout omits it entirely. " +
                "mzLib substitutes zero, which is indistinguishable from a real measurement of " +
                "nothing, so the value crosses as null instead. read-records has the file's own " +
                "summed Intensity column.");
        }
        else if (unmeasured > 0)
        {
            caveats.Add(
                $"intensity is null for {unmeasured} of {features.Count} rows: those features carry " +
                "no Apex_intensity value, and mzLib substitutes zero (Ms1Feature.cs:86). Null here " +
                "means 'not reported', not 'no signal'.");
        }

        return caveats;
    }

    private static List<string> BaseFeatureCaveatsFor(SupportedFileType fileType) => fileType switch
    {
        SupportedFileType.Ms1Feature =>
        [
            "One row here is one CHARGE STATE of one feature, not one row of the file. mzLib expands " +
            "each deconvolved feature across [ChargeStateMin, ChargeStateMax] (Ms1Feature.cs:84), so " +
            "record_count exceeds the file's line count, and a charge the tool never observed appears " +
            "if the writer recorded a gapped charge range. Use read-records for the file's own rows.",
            "intensity is the per-charge APEX intensity (Ms1Feature.cs:86), not the summed intensity " +
            "over the feature. The file's own Intensity column is a different number, and " +
            "read-records has it.",
            "retention_time_start/_end are in UNKNOWN units for this format. TopFD wrote seconds " +
            "through v1.6.2 and minutes from v1.7.0 without changing the file type, and mzLib " +
            "normalises neither - its deconvolution parameters instead GUESS, dividing by 60 when " +
            "the largest end time exceeds 500. Check the values against your gradient length before " +
            "comparing them with anything.",
            "number_of_isotopes is null for every row of this format: the single-charge expansion " +
            "mzLib builds never sets it (Ms1Feature.cs:91). Null means 'not reported', not " +
            "'no isotopes found'.",
        ],
        SupportedFileType.Tsv_Dinosaur =>
        [
            // Deliberately NOT the expansion caveat: DinosaurTsvFile.GetMs1Features is `=> Results`
            // (DinosaurTsvFile.cs:15), so one row here is exactly one row of the file.
            "intensity is Dinosaur's intensityApex column, not intensitySum (DinosaurFeature.cs:18). " +
            "Both are in the file and read-records has both; this view can only carry one, and mzLib " +
            "chose the apex.",
            "mz is the feature's monoisotopic m/z. Dinosaur also reports mostAbundantMz, which for a " +
            "peptide above roughly 1.8 kDa is a different isotope; read-records has it.",
        ],
        _ => [],
    };

    // ---------------------------------------------------------------------------------------
    // spectral_match
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The <see cref="ISpectralMatch"/> fields, under mzLib's own names.
    /// </summary>
    /// <remarks>
    /// This is <c>Readers.ISpectralMatch</c> — modifications and no score — not the unrelated
    /// <c>Omics.SpectralMatch.ISpectralMatch</c>, which has a score and no modifications. mzLib
    /// carries both names and aliases them in its own source to tell them apart.
    /// </remarks>
    private static IReadOnlyList<Column<ISpectralMatch>> MatchColumns { get; } = new[]
    {
        new Column<ISpectralMatch>("file_name_without_extension", m => m.FileNameWithoutExtension),
        new Column<ISpectralMatch>("one_based_scan_number", m => m.OneBasedScanNumber),
        new Column<ISpectralMatch>("base_sequence", m => m.BaseSequence),
        new Column<ISpectralMatch>("full_sequence", m => m.FullSequence),
        new Column<ISpectralMatch>("accession", m => m.Accession),
        // Null, not false, where the format cannot report decoys — the same rule the quantifiable
        // view applies to MSFragger, and for the same reason. See DecoysAreReported.
        new Column<ISpectralMatch>("is_decoy", m => DecoysAreReported(m) ? m.IsDecoy : (bool?)null),
        // The modification dictionary flattened to "position:name" pairs, ';'-joined, in position
        // order. Position 1 is the N-terminus, following mzLib's own one-is-N-terminus convention —
        // renumbering it here would silently disagree with every mzLib document about the same file.
        new Column<ISpectralMatch>("modifications", FormatModifications),
        new Column<ISpectralMatch>("modification_count", m => m.AllModsOneIsNterminus?.Count ?? 0),
        // The three confidence fields a format may record, under mzLib's names. None is on
        // Readers.ISpectralMatch, so each is read from the record types that carry it, and a
        // format with no source for one names it in absent_fields (MatchAbsentFields), where the
        // engine nulls it in every row. mzIdentML's q-value is the PSM-level MS:1002354 or an
        // engine-specific child of it, null on an item that reports none (mzLib #1306).
        new Column<ISpectralMatch>("q_value", m => m switch
        {
            MzIdentMLRecord record => record.QValue,
            MsPathFinderTResult result => result.QValue,
            _ => null,
        }),
        new Column<ISpectralMatch>("rank", m => (m as MzIdentMLRecord)?.Rank),
        new Column<ISpectralMatch>("pass_threshold", m => (m as MzIdentMLRecord)?.PassThreshold),
    };

    /// <summary>
    /// The spectral-match fields this file's format has no source for: null in every row, and
    /// named, so a missing value is never mistaken for a measured one.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <c>is_decoy</c> is absent for every format but MsPathFinderT (see
    /// <see cref="DecoysAreReported"/>): Casanovo writes no label, and mzIdentML's <c>isDecoy</c>
    /// defaults to false when a writer omits it, so mzLib cannot tell a stated false from an unstated
    /// one (MzIdentMLResultFile.cs:173).
    /// </para>
    /// <para>
    /// <c>q_value</c> is absent for Casanovo, and for an MsPathFinderT file without a
    /// <c>QValue</c> column — the <c>_IcTarget.tsv</c> and <c>_IcDecoy.tsv</c> files, which
    /// MSPathFinder writes before target-decoy analysis. mzLib types that property as a
    /// non-nullable double, so it reads 0 there: a perfect q-value for every match, which is the
    /// most dangerous number this view could hand back. <c>rank</c> and <c>pass_threshold</c> exist
    /// only in mzIdentML.
    /// </para>
    /// </remarks>
    private static List<string> MatchAbsentFields(IResultFile resultFile, string path)
    {
        Type? recordType = RecordTypeOf(resultFile.GetType());
        bool msPathFinder = recordType is not null && typeof(MsPathFinderTResult).IsAssignableFrom(recordType);
        bool mzid = resultFile is MzIdentMLResultFile;

        var absent = new List<string>();
        if (!msPathFinder)
            absent.Add("is_decoy");

        bool qValueRead = mzid
            || (msPathFinder && !AbsentProperties(recordType, path).Contains(nameof(MsPathFinderTResult.QValue)));
        if (!qValueRead)
            absent.Add("q_value");

        if (!mzid)
            absent.AddRange(["rank", "pass_threshold"]);

        return absent;
    }

    /// <summary>Whether this record's format can report decoy status at all.</summary>
    /// <remarks>
    /// <c>CasanovoMzTabRecord.IsDecoy</c> is an auto-property initialised to <c>false</c> that the
    /// reader never assigns (CasanovoMzTabRecord.cs:84) — Casanovo is de novo and emits no
    /// target/decoy label at all — so <c>false</c> there means "mzLib cannot tell", exactly as it
    /// does for MSFragger in the quantifiable view. An ALLOWLIST for the same reason: a reader
    /// added later that also leaves the field at its default must not start emitting fabricated
    /// booleans without anyone editing this file.
    /// </remarks>
    private static bool DecoysAreReported(ISpectralMatch match) => match is MsPathFinderTResult;

    /// <summary>
    /// The <c>ms1_features</c> fields this file has no source for.
    /// </summary>
    /// <remarks>
    /// <c>intensity</c> is absent when the file's header has no <c>Apex_intensity</c> column — every
    /// FLASHDeconv/OpenMS <c>_ms1.feature</c> — because mzLib then substitutes zero for every
    /// feature (Ms1Feature.cs:86). <c>number_of_isotopes</c> is absent for the whole
    /// <c>_ms1.feature</c> format: the single-charge expansion never sets it (Ms1Feature.cs:91).
    /// </remarks>
    private static List<string> FeatureAbsentFields(IResultFile resultFile, string path)
    {
        var absent = new List<string>();
        if (resultFile is Ms1FeatureFile)
        {
            if (AbsentProperties(typeof(Ms1Feature), path).Contains(nameof(Ms1Feature.IntensityApex)))
                absent.Add("intensity");
            absent.Add("number_of_isotopes");
        }

        return absent;
    }

    /// <summary>
    /// What a spectra file says about the run: instrument, serial number, and when acquisition
    /// started (mzLib <see cref="SourceFile"/>, #1349).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <c>instrument_model</c> is the model's NAME and <c>instrument_model_accession</c> its PSI-MS
    /// accession, kept apart because the two readers that fill them do so differently: mzML carries
    /// both, Thermo <c>.raw</c> only the name. mzLib's own rule is to match on the accession, never on
    /// the name, so the accession is null rather than guessed when the file has none.
    /// </para>
    /// <para>
    /// <c>acquisition_start_time</c> is ISO-8601. It ends in <c>Z</c> only when the source fixed the
    /// instant (an mzML <c>startTimeStamp</c> with an offset); otherwise it is the acquisition
    /// computer's wall-clock time with no offset, because the instant is unknown and none is
    /// invented. <c>acquisition_start_time_is_utc</c> says which, and is always false for Thermo
    /// <c>.raw</c>. A <c>.raw</c> and ProteoWizard's mzML of it can therefore differ by the site's UTC
    /// offset: ProteoWizard assumes the converting machine's time zone and writes <c>Z</c>.
    /// </para>
    /// <para>
    /// Null as a whole only when the reader built no source description at all; a blank field is
    /// null individually, meaning the file does not record it.
    /// </para>
    /// </remarks>
    private static object? SourceOf(MsDataFile file)
    {
        SourceFile? source = file.SourceFile;
        if (source is null)
            return null;

        DateTime? started = source.AcquisitionStartTime;
        return new
        {
            instrument_model = NullIfBlank(source.InstrumentModel?.Name),
            instrument_model_accession = NullIfBlank(source.InstrumentModel?.Accession),
            instrument_serial_number = NullIfBlank(source.InstrumentSerialNumber),
            acquisition_start_time = started is { } moment ? IsoTimestamp(moment) : null,
            acquisition_start_time_is_utc = started?.Kind == DateTimeKind.Utc,
        };
    }

    private static string? NullIfBlank(string? value) => string.IsNullOrWhiteSpace(value) ? null : value.Trim();

    /// <summary>
    /// ISO-8601, to the second or the microsecond, with <c>Z</c> only for a UTC instant.
    /// </summary>
    /// <remarks>
    /// Never seven fractional digits (.NET's round-trip "O"), which Python before 3.11 cannot parse,
    /// and never a fraction of zeros that no instrument recorded.
    /// </remarks>
    private static string IsoTimestamp(DateTime moment)
    {
        string format = moment.Ticks % TimeSpan.TicksPerSecond == 0
            ? "yyyy-MM-dd'T'HH:mm:ss"
            : "yyyy-MM-dd'T'HH:mm:ss.ffffff";
        return moment.ToString(format, CultureInfo.InvariantCulture) + (moment.Kind == DateTimeKind.Utc ? "Z" : "");
    }

    /// <summary>
    /// The verb that carries a field <c>read-records</c> cannot, or null when none does.
    /// </summary>
    /// <remarks>
    /// So an <c>excluded_fields</c> entry is a pointer, not a dead end: the per-sample tables of the
    /// two mzLib 1.0.592 quantification readers (#1347) and mzIdentML's engine scores (#1306) each
    /// have a typed verb that projects them in long form.
    /// </remarks>
    private static string? VerbCarrying(Type recordType, string propertyName) => (recordType.Name, propertyName) switch
    {
        (nameof(ProteinGroupFromTsv), nameof(ProteinGroupFromTsv.SampleGroups)) => "readers read-protein-groups",
        (nameof(QuantifiedPeptideFromTsv), nameof(QuantifiedPeptideFromTsv.Samples)) => "readers read-quantified-peptides",
        (nameof(MzIdentMLRecord), nameof(MzIdentMLRecord.Scores)) => "readers read-matches",
        _ => null,
    };

    private static string FormatModifications(ISpectralMatch match)
    {
        Dictionary<int, Omics.Modifications.Modification>? mods = match.AllModsOneIsNterminus;
        if (mods is null || mods.Count == 0)
            return string.Empty;

        return string.Join(";", mods
            .OrderBy(pair => pair.Key)
            .Select(pair => $"{pair.Key}:{pair.Value?.IdWithMotif ?? pair.Value?.OriginalId ?? "unknown"}"));
    }

    /// <summary>What the spectral-match view cannot be trusted to mean, per format.</summary>
    private static List<string> MatchCaveatsFor(SupportedFileType fileType)
    {
        var caveats = new List<string>
        {
            "NOTHING here is FDR-filtered. q_value is the only confidence field this view carries, " +
            "and only mzIdentML and an MsPathFinderT file with a QValue column fill it (absent_fields " +
            "says when it is empty). Every one of these formats records scores this view does not " +
            "expose; read-records has them. Filter before you report.",
        };

        switch (fileType)
        {
            case SupportedFileType.MsPathFinderTTargets:
            case SupportedFileType.MsPathFinderTDecoys:
            case SupportedFileType.MsPathFinderTAllResults:
                caveats.Add(
                    "is_decoy is derived from the protein NAME, not from a column: mzLib reports a " +
                    "decoy when ProteinName starts with 'XXX' (MsPathFinderTResult.cs:92). A database " +
                    "whose decoys carry a different prefix reads entirely as targets.");
                caveats.Add(
                    "accession is ProteinName split on '|' taking the second field " +
                    "(MsPathFinderTResult.cs:89), which assumes a UniProt-style header. A FASTA with " +
                    "plain headers makes that property throw, and the cell arrives null with the " +
                    "field named in failed_fields rather than the whole file failing to read.");
                caveats.Add(
                    "modifications are resolved by a FUZZY match against mzLib's modification " +
                    "dictionary (ModificationConverter.GetClosestMod), not by an exact identifier, so " +
                    "a name here can differ from the one the search engine wrote.");
                break;

            case SupportedFileType.CasanovoMzTab:
                caveats.Add(
                    "is_decoy is null for this format. Casanovo is de novo and writes no target/decoy " +
                    "label; mzLib's record leaves the field at its default false and never assigns it " +
                    "(CasanovoMzTabRecord.cs:84), so false would mean 'unknown', not 'target'.");
                caveats.Add(
                    "one_based_scan_number is the mzTab spectrum INDEX plus one, not necessarily the " +
                    "instrument's scan number (CasanovoMzTabFile.cs:116). When Casanovo was run on an " +
                    "MGF the two are unrelated, so do not join this against a raw file on scan number.");
                caveats.Add(
                    "full_sequence and modifications are resolved by matching Casanovo's mass shifts " +
                    "against mzLib's modification dictionary (CasanovoMzTabFile.cs:124), not read " +
                    "from named annotations - Casanovo writes none. An empty value therefore means " +
                    "the peptide is unmodified, but a populated one is mzLib's interpretation of a " +
                    "mass, not the search engine's own call.");
                break;

            case SupportedFileType.MzIdentML:
            case SupportedFileType.MzIdentMLGz:
                caveats.Add(
                    "is_decoy is null for this format. mzIdentML's isDecoy attribute is optional and " +
                    "defaults to false, and mzLib reports a decoy only when every peptide evidence says " +
                    "so (MzIdentMLResultFile.cs:173), so false cannot be told apart from 'not stated'. " +
                    "read-records carries mzLib's boolean for a caller who knows the writer sets it.");
                caveats.Add(
                    "Every SpectrumIdentificationItem is a row, not only the matches the submitter " +
                    "accepted: lower-ranked candidates and items that fail the threshold are here too " +
                    "(MzIdentMLResultFile.cs:177). Filter on rank == 1 and pass_threshold before " +
                    "counting identifications.");
                caveats.Add(
                    "one_based_scan_number is parsed from the nativeID (MzIdentMLResultFile.cs:159). " +
                    "'scan=N' gives N, but 'index=N', which peak-list input carries, is a zero-based " +
                    "position in the file and gives N + 1, not an instrument scan number. -1 means the " +
                    "nativeID had neither.");
                caveats.Add(
                    "Items mzLib cannot represent as one linear match are skipped, not failed: " +
                    "crosslinks, modifications without a resolvable UNIMOD accession, substitutions, and " +
                    "two modifications on one residue (MzIdentMLResultFile.cs:123). They are not rows; " +
                    "skipped_count and skipped name each one and why, so record_count plus skipped_count " +
                    "is the number of items in the file.");
                caveats.Add(
                    "The engine's own scores (for example MS-GF:SpecEValue) have no common name across " +
                    "search engines (MzIdentMLResultFile.cs:179). Pass scores=true for them as long " +
                    "rows, one per match and score. q_value is null on an item that reports none.");
                caveats.Add(
                    "accession joins every protein the item's peptide evidence names with '|' " +
                    "(MzIdentMLRecord.cs:45), the same character a UniProt header uses inside one " +
                    "accession, so the cell cannot be split back into proteins reliably.");
                break;
        }

        return caveats;
    }

    // ---------------------------------------------------------------------------------------
    // spectra
    // ---------------------------------------------------------------------------------------

    /// <summary>The scan-header fields of an <see cref="MsDataScan"/>.</summary>
    /// <remarks>
    /// Header only: the peak arrays are excluded here and added by
    /// <see cref="ScanColumnsWithPeaks"/> under <c>--peaks</c>. The precursor fields are all
    /// genuinely nullable on <see cref="MsDataScan"/> and are left null for MS1 scans rather than
    /// zero-filled — a precursor m/z of 0 is a number someone will plot.
    /// </remarks>
    private static IReadOnlyList<Column<MsDataScan>> ScanColumns { get; } = BuildScanColumns();

    /// <summary>The scan headers plus the peak arrays.</summary>
    private static IReadOnlyList<Column<MsDataScan>> ScanColumnsWithPeaks { get; } =
    [
        .. BuildScanColumns(),
        // One array per scan, so the column is a list of lists in the JSON envelope and a
        // ';'-joined list per cell in a written table.
        new Column<MsDataScan>("mz", s => s.MassSpectrum?.XArray),
        new Column<MsDataScan>("intensity", s => s.MassSpectrum?.YArray),
    ];

    private static List<Column<MsDataScan>> BuildScanColumns() =>
    [
        new Column<MsDataScan>("one_based_scan_number", s => s.OneBasedScanNumber),
        new Column<MsDataScan>("ms_order", s => s.MsnOrder),
        // Minutes for every MsDataFile reader in mzLib — the spectra readers convert at the
        // boundary, unlike the result-file readers, which pass the tool's unit through unchanged.
        new Column<MsDataScan>("retention_time", s => s.RetentionTime),
        new Column<MsDataScan>("polarity", s => s.Polarity.ToString()),
        new Column<MsDataScan>("mz_analyzer", s => s.MzAnalyzer.ToString()),
        new Column<MsDataScan>("is_centroid", s => s.IsCentroid),
        new Column<MsDataScan>("total_ion_current", s => s.TotalIonCurrent),
        new Column<MsDataScan>("injection_time", s => s.InjectionTime),
        new Column<MsDataScan>("peak_count", s => s.MassSpectrum?.Size ?? 0),
        new Column<MsDataScan>("scan_window_lower_mz", s => s.ScanWindowRange?.Minimum),
        new Column<MsDataScan>("scan_window_upper_mz", s => s.ScanWindowRange?.Maximum),
        new Column<MsDataScan>("scan_filter", s => s.ScanFilter),
        new Column<MsDataScan>("native_id", s => s.NativeId),
        new Column<MsDataScan>("scan_description", s => s.ScanDescription),
        // Precursor fields — null on an MS1 scan, and null rather than zero wherever the file did
        // not record them.
        new Column<MsDataScan>("one_based_precursor_scan_number", s => s.OneBasedPrecursorScanNumber),
        new Column<MsDataScan>("isolation_mz", s => s.IsolationMz),
        new Column<MsDataScan>("isolation_width", s => s.IsolationWidth),
        new Column<MsDataScan>("selected_ion_mz", s => s.SelectedIonMZ),
        new Column<MsDataScan>("selected_ion_intensity", s => s.SelectedIonIntensity),
        new Column<MsDataScan>("selected_ion_charge_state_guess", s => s.SelectedIonChargeStateGuess),
        new Column<MsDataScan>("selected_ion_monoisotopic_guess_mz", s => s.SelectedIonMonoisotopicGuessMz),
        new Column<MsDataScan>("dissociation_type", s => s.DissociationType?.ToString()),
        new Column<MsDataScan>("hcd_energy", s => s.HcdEnergy),
        // FAIMS. Null on every instrument without it, which is most of them.
        new Column<MsDataScan>("compensation_voltage", s => s.CompensationVoltage),
    ];

    /// <summary>What the spectra view cannot be trusted to mean, for a given file.</summary>
    private static List<string> SpectraCaveatsFor(string path, bool includePeaks)
    {
        var caveats = new List<string>();

        SupportedFileType fileType;
        try
        {
            fileType = path.ParseFileType();
        }
        catch (MzLibUtil.MzLibException)
        {
            return caveats;
        }

        if (!includePeaks)
        {
            caveats.Add(
                "Peaks are not included. This is scan HEADERS only; peak_count reports how many " +
                "peaks each scan has but not what they are. Pass peaks=true to include the mz and " +
                "intensity arrays, and expect the payload to grow by roughly the size of the file.");
        }

        switch (fileType)
        {
            case SupportedFileType.Mgf:
                caveats.Add(
                    "MGF's 'scan numbers' come from the SCANS line when the writer supplied one, and are " +
                    "otherwise assigned in file order, not by an instrument, so they need not be " +
                    "contiguous or even unique. one_based_precursor_scan_number is null unless the file " +
                    "carries PRECURSORSCAN, an mzLib extension header (Mgf.cs:297): standard MGF does not " +
                    "record which survey scan a fragment scan came from, so expect null on files mzLib " +
                    "did not write. ms_order is no longer always 2. mzLib takes it from the MSLEVEL line " +
                    "when the writer supplied one, and otherwise reads a block with a precursor as MS2 and " +
                    "a block without one as MS1 (Mgf.cs:350). Files written before MSLEVEL existed all " +
                    "carry PEPMASS, so they still read as MS2 throughout.");
                caveats.Add(
                    "scan_window_lower_mz/_upper_mz are DERIVED, not recorded: MGF has no scan-window " +
                    "field, so mzLib reports the first and last observed peak (Mgf.cs:346). They are " +
                    "the fragment m/z range actually seen, which is narrower than the instrument's " +
                    "window and depends on the peak-picking threshold.");
                break;

            case SupportedFileType.Ms1Align:
            case SupportedFileType.Ms2Align:
                caveats.Add(
                    "msalign holds DECONVOLVED masses, not raw m/z. The mz column carries neutral " +
                    "monoisotopic masses that a deconvolution step already produced, so it is not " +
                    "comparable with the mz of an mzML or raw file and must not be re-deconvolved.");
                caveats.Add(
                    "scan_window_lower_mz/_upper_mz are in m/z while the mz column is in neutral " +
                    "MASS, so the two are not on the same axis and filtering peaks to the window " +
                    "would discard most of the spectrum. mzLib synthesises the window by converting " +
                    "each mass back to m/z at its reported charge (MsAlign.cs:526); msalign records " +
                    "no window of its own.");
                break;

            case SupportedFileType.BrukerD:
            case SupportedFileType.BrukerTimsTof:
                caveats.Add(
                    "Bruker formats are read through vendor native libraries (baf2sql / timsdata) " +
                    "and are Windows-x64 only. On any other platform this verb fails rather than " +
                    "returning an empty file.");
                break;

            case SupportedFileType.ThermoRaw:
                caveats.Add(
                    "Thermo .raw is read through Thermo's RawFileReader, whose licence terms you " +
                    "accept by using it. The reader is managed rather than native, so it works on " +
                    "any platform.");
                break;
        }

        if (fileType is SupportedFileType.BrukerTimsTof)
        {
            caveats.Add(
                "timsTOF data is ion-mobility-resolved and this view flattens it: each frame's " +
                "mobility dimension is collapsed into scans, so a 1/K0 value is not reported. Use " +
                "mzLib directly if you need the mobility axis.");
        }

        return caveats;
    }
}
