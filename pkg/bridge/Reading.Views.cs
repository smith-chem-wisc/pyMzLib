using MassSpectrometry;
using Readers;

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
    private static IReadOnlyList<Column<ISingleChargeMs1Feature>> FeatureColumns { get; } = new[]
    {
        new Column<ISingleChargeMs1Feature>("mz", f => f.Mz),
        new Column<ISingleChargeMs1Feature>("charge", f => f.Charge),
        new Column<ISingleChargeMs1Feature>("retention_time_start", f => f.RetentionTimeStart),
        new Column<ISingleChargeMs1Feature>("retention_time_end", f => f.RetentionTimeEnd),
        // Passed through exactly as mzLib reports it, INCLUDING the fabricated zero it substitutes
        // when a file has no Apex_intensity column. That zero is disclosed in the caveats rather
        // than repaired here: repairing it would make the wire disagree with mzLib about a number,
        // which is a different kind of change from adding a verb and belongs in its own PR — and
        // ultimately upstream, where the interface should be nullable (bridge/UPSTREAM.md U1).
        new Column<ISingleChargeMs1Feature>("intensity", f => f.Intensity),
        // Genuinely nullable on the interface, and null for a whole format rather than for odd
        // rows: mzLib's _ms1.feature expansion never sets it. Crossing as null is the faithful
        // projection — no repair involved, because the interface already says it may be absent.
        new Column<ISingleChargeMs1Feature>("number_of_isotopes", f => f.NumberOfIsotopes),
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
    private static List<string> FeatureCaveatsFor(IResultFile resultFile)
    {
        List<string> caveats = BaseFeatureCaveatsFor(resultFile.FileType);

        // Said only when it is true of THIS file, because it depends on the schema the writer used
        // rather than on the file type: TopFD writes Apex_intensity and FLASHDeconv does not, and
        // both are SupportedFileType.Ms1Feature. Tested against the file's own records, which is
        // the same check mzLib itself uses to tell the two writers apart (Ms1FeatureFile.cs:50).
        if (resultFile is Ms1FeatureFile file &&
            file.Results.Count > 0 &&
            file.Results.All(record => record.IntensityApex is null))
        {
            caveats.Add(
                "intensity is 0 for EVERY row of this file, and that zero is FABRICATED. mzLib takes " +
                "the per-charge intensity from the optional Apex_intensity column " +
                "(Ms1Feature.cs:86); this file's schema does not have it - the FLASHDeconv/OpenMS " +
                "_ms1.feature layout omits it entirely - so mzLib substitutes zero, which is " +
                "indistinguishable from a real measurement of nothing. Do not read these as " +
                "intensities, and do not sum or rank them. read-records has the file's own summed " +
                "Intensity column, which is real.");
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
    };

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
            "There is no score, E-value or q-value in this view, so NOTHING here is FDR-filtered. " +
            "Readers.ISpectralMatch carries identity fields only; every one of these formats records " +
            "confidence in columns this view does not expose. read-records has them. Filter before " +
            "you report.",
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
                    "MGF carries no MS1 scans, and its 'scan numbers' come from the TITLE line " +
                    "rather than from an instrument, so they need not be contiguous or even unique. " +
                    "one_based_precursor_scan_number is always null: the format does not record " +
                    "which survey scan a fragment scan came from.");
                caveats.Add(
                    "scan_window_lower_mz/_upper_mz are DERIVED, not recorded: MGF has no scan-window " +
                    "field, so mzLib reports the first and last observed peak (Mgf.cs:221). They are " +
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
