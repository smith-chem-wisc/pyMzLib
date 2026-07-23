using FlashLFQ;
using MassSpectrometry;
using Readers;

namespace MzLibBridge;

/// <summary>
/// The label-free quantification workflow: read a PSM result file, quantify its peptides across a
/// set of mzML runs with FlashLFQ, and report per-file peptide and protein intensities.
/// </summary>
/// <remarks>
/// <para>
/// This is ONE verb answering the question a quant workflow actually asks — "given these
/// identifications and these runs, how much of each peptide and protein is in each run?" — rather
/// than exposing FlashLFQ's object graph. The heavy lifting is entirely mzLib's: the result file is
/// read by <see cref="FileReader.ReadQuantifiableResultFile"/> (Readers), turned into FlashLFQ
/// <see cref="Identification"/>s by <see cref="MzLibExtensions.MakeIdentifications"/> (FlashLFQ), and
/// quantified by <see cref="FlashLfqEngine"/>. The bridge writes none of that — its only job is to
/// build the <see cref="SpectraFileInfo"/> list (the experimental-design surface), call the engine,
/// and flatten the results onto the wire. MetaMorpheus is not involved: mzLib does it all alone.
/// </para>
/// <para>
/// mzML-only for now. The Thermo <c>.raw</c> and Bruker readers that <c>Readers</c> can dispatch to
/// are deliberately not exercised, and a non-mzML path is rejected up front rather than failing deep
/// inside the indexing engine.
/// </para>
/// </remarks>
internal static class Quantification
{
    /// <summary>
    /// <c>quant flashlfq --psms PATH [--normalize] [--ppm 10] [--isotope-ppm 5] [--integrate]
    /// [--mbr] [--mbr-ppm 10] [--mbr-q 0.05] [--shared-peptides] [--bayesian] [--use-pep-q]
    /// [--threads N] [--out DIR]</c>
    /// </summary>
    /// <remarks>
    /// The spectra files come on <b>stdin</b>, one per line, tab-separated:
    /// <c>path[\tcondition[\tbiorep[\ttechrep[\tfraction]]]]</c>. stdin rather than argv because a
    /// real experiment has many runs and argv has a hard size ceiling — the same reason the PRIDE
    /// download verb takes its selection on stdin. Trailing design fields default the way
    /// MetaMorpheus defaults them with no experimental-design file: blank condition, each file its
    /// own biological replicate, fraction 0, technical replicate 0.
    /// </remarks>
    public static object FlashLfq(Program.Arguments arguments)
    {
        string psmPath = arguments.Required("psms");

        var flashParams = new FlashLfqParameters
        {
            Normalize = arguments.Flag("normalize"),
            PpmTolerance = arguments.OptionalDouble("ppm", 10.0),
            IsotopePpmTolerance = arguments.OptionalDouble("isotope-ppm", 5.0),
            Integrate = arguments.Flag("integrate"),
            MatchBetweenRuns = arguments.Flag("mbr"),
            MbrPpmTolerance = arguments.OptionalDouble("mbr-ppm", 10.0),
            MbrQValueThreshold = arguments.OptionalDouble("mbr-q", 0.05),
            UseSharedPeptidesForProteinQuant = arguments.Flag("shared-peptides"),
            BayesianProteinQuant = arguments.Flag("bayesian"),
            MaxThreads = arguments.OptionalInt("threads", -1),
            // The engine writes progress to the console when not silent. stdout carries only the JSON
            // envelope (the bridge's core contract), so it must stay quiet; diagnostics belong on
            // stderr, which the console redirect below guarantees regardless of this flag.
            Silent = true,
        };

        bool usePepQValue = arguments.Flag("use-pep-q");
        string? outputDirectory = arguments.Optional("out");

        List<SpectraFileInfo> spectraFiles = BuildSpectraFiles(ReadStdinLines());

        IQuantifiableResultFile resultFile;
        try
        {
            resultFile = FileReader.ReadQuantifiableResultFile(psmPath);
        }
        catch (FileNotFoundException)
        {
            throw new Program.UsageException($"PSM result file not found: '{psmPath}'.");
        }

        // Every identification must map to a provided run, or MakeIdentifications throws mid-stream
        // with a bare "Spectra file not found". Checking up front turns that into one clear message
        // naming exactly which mzML files are missing, before any spectra are read.
        var providedNames = spectraFiles.Select(f => f.FilenameWithoutExtension).ToHashSet(StringComparer.Ordinal);
        List<string> missing = resultFile.GetQuantifiableResults()
            .Select(r => r.FileName)
            .Distinct(StringComparer.Ordinal)
            .Where(name => !providedNames.Contains(name))
            .OrderBy(name => name, StringComparer.Ordinal)
            .ToList();
        if (missing.Count > 0)
            throw new Program.UsageException(
                $"The PSM file identifies peptides in {missing.Count} run(s) with no mzML provided: " +
                $"{string.Join(", ", missing.Take(5))}{(missing.Count > 5 ? ", …" : string.Empty)}. " +
                "Provide the matching mzML for each (matched by base file name).");

        List<Identification> identifications = resultFile.MakeIdentifications(spectraFiles, usePepQValue);

        // The engine reads the mzML runs by path through Readers' extension-dispatch factory; only
        // the requested-and-missing case above needs guarding. Any progress it would print is kept
        // off stdout by redirecting the console for the duration of the run.
        TextWriter originalOut = Console.Out;
        Console.SetOut(Console.Error);
        FlashLfqResults results;
        try
        {
            results = new FlashLfqEngine(flashParams, identifications).Run();
        }
        finally
        {
            Console.SetOut(originalOut);
        }

        if (!string.IsNullOrWhiteSpace(outputDirectory))
        {
            Directory.CreateDirectory(outputDirectory);
            results.WriteResults(
                Path.Combine(outputDirectory, "QuantifiedPeaks.tsv"),
                Path.Combine(outputDirectory, "QuantifiedPeptides.tsv"),
                Path.Combine(outputDirectory, "QuantifiedProteins.tsv"),
                flashParams.BayesianProteinQuant ? Path.Combine(outputDirectory, "BayesianProteinQuant.tsv") : null,
                silent: true);
        }

        return ToWire(psmPath, flashParams, identifications, spectraFiles, results, outputDirectory);
    }

    /// <summary>Reads stdin into non-blank trimmed lines.</summary>
    private static List<string> ReadStdinLines()
    {
        var lines = new List<string>();
        string? line;
        while ((line = Console.In.ReadLine()) != null)
        {
            if (!string.IsNullOrWhiteSpace(line))
                lines.Add(line);
        }
        return lines;
    }

    /// <summary>
    /// Turns the stdin spectra lines into <see cref="SpectraFileInfo"/>s, applying the
    /// MetaMorpheus-style defaults for any omitted experimental-design field.
    /// </summary>
    /// <remarks>
    /// Base file names must be unique: FlashLFQ matches an identification to its run by base name,
    /// so two runs sharing one would be indistinguishable to both the engine and the wire output.
    /// </remarks>
    internal static List<SpectraFileInfo> BuildSpectraFiles(IReadOnlyList<string> lines)
    {
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No spectra files were provided on stdin. Supply one mzML path per line: " +
                "'path[<TAB>condition[<TAB>biorep[<TAB>techrep[<TAB>fraction]]]]'.");

        var files = new List<SpectraFileInfo>(lines.Count);
        var seenNames = new HashSet<string>(StringComparer.Ordinal);

        for (int i = 0; i < lines.Count; i++)
        {
            string[] parts = lines[i].Split('\t');
            string path = parts[0].Trim();
            if (path.Length == 0)
                throw new Program.UsageException($"Spectra line {i + 1} has no file path.");

            string extension = Path.GetExtension(path);
            if (!extension.Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                throw new Program.UsageException(
                    $"Only mzML is supported for now; '{path}' has extension " +
                    $"'{(extension.Length == 0 ? "none" : extension)}'. Convert .raw/.d to mzML first.");

            if (!File.Exists(path))
                throw new Program.UsageException($"Spectra file not found: '{path}'.");

            string name = Path.GetFileNameWithoutExtension(path);
            if (!seenNames.Add(name))
                throw new Program.UsageException(
                    $"Two spectra files share the base name '{name}'. Base names must be unique — " +
                    "FlashLFQ matches identifications to runs by base name.");

            string condition = parts.Length > 1 ? parts[1].Trim() : string.Empty;
            int biologicalReplicate = ParseDesignField(parts, 2, i, "biorep", fallback: i);
            int technicalReplicate = ParseDesignField(parts, 3, i, "techrep", fallback: 0);
            int fraction = ParseDesignField(parts, 4, i, "fraction", fallback: 0);

            files.Add(new SpectraFileInfo(path, condition,
                biorep: biologicalReplicate, techrep: technicalReplicate, fraction: fraction));
        }

        return files;
    }

    private static int ParseDesignField(string[] parts, int index, int lineIndex, string field, int fallback)
    {
        if (parts.Length <= index || string.IsNullOrWhiteSpace(parts[index]))
            return fallback;
        if (!int.TryParse(parts[index].Trim(), out int value) || value < 0)
            throw new Program.UsageException(
                $"Spectra line {lineIndex + 1}: {field} must be a non-negative integer; got '{parts[index]}'.");
        return value;
    }

    /// <summary>Flattens the quantification results onto the wire.</summary>
    private static object ToWire(
        string psmPath,
        FlashLfqParameters flashParams,
        List<Identification> identifications,
        List<SpectraFileInfo> spectraFiles,
        FlashLfqResults results,
        string? outputDirectory)
    {
        return new
        {
            psm_file = Path.GetFullPath(psmPath),
            identification_count = identifications.Count,
            parameters = new
            {
                normalize = flashParams.Normalize,
                ppm_tolerance = flashParams.PpmTolerance,
                isotope_ppm_tolerance = flashParams.IsotopePpmTolerance,
                integrate = flashParams.Integrate,
                match_between_runs = flashParams.MatchBetweenRuns,
                mbr_ppm_tolerance = flashParams.MbrPpmTolerance,
                mbr_q_value_threshold = flashParams.MbrQValueThreshold,
                use_shared_peptides_for_protein_quant = flashParams.UseSharedPeptidesForProteinQuant,
                bayesian_protein_quant = flashParams.BayesianProteinQuant,
                max_threads = flashParams.MaxThreads,
            },
            spectra_files = spectraFiles.Select(file => new
            {
                file_name = file.FilenameWithoutExtension,
                full_path = file.FullFilePathWithExtension,
                condition = file.Condition,
                biological_replicate = file.BiologicalReplicate,
                technical_replicate = file.TechnicalReplicate,
                fraction = file.Fraction,
                peak_count = results.Peaks[file].Count,
                // MBR is off by default; when on, this is the count of peaks quantified from a run
                // where the peptide was never identified — the whole point of match-between-runs.
                mbr_peak_count = results.Peaks[file].Count(p => p.DetectionType == DetectionType.MBR),
            }).ToList(),
            peptide_count = results.PeptideModifiedSequences.Count,
            protein_count = results.ProteinGroups.Count,
            peptides = results.PeptideModifiedSequences.Values.Select(peptide => new
            {
                sequence = peptide.Sequence,
                base_sequence = peptide.BaseSequence,
                protein_groups = string.Join(";", peptide.ProteinGroups.Select(g => g.ProteinGroupName).Distinct()),
                intensities = spectraFiles.ToDictionary(
                    file => file.FilenameWithoutExtension, file => Finite(peptide.GetIntensity(file))),
                detection_types = spectraFiles.ToDictionary(
                    file => file.FilenameWithoutExtension, file => peptide.GetDetectionType(file).ToString()),
            }).ToList(),
            proteins = results.ProteinGroups.Values.Select(group => new
            {
                protein_group = group.ProteinGroupName,
                gene_name = group.GeneName,
                organism = group.Organism,
                // Protein intensities can be NaN by design: FlashLFQ's median-polish protein quant
                // marks a protein NaN when the peptide matrix is degenerate — several files reporting
                // the same intensity, or too few peptides per file to resolve (a real artifact,
                // documented in mzLib's own FlashLFQ tests). NaN is not valid JSON, so it crosses as
                // null — "could not be quantified" — rather than crashing serialization.
                intensities = spectraFiles.ToDictionary(
                    file => file.FilenameWithoutExtension, file => Finite(group.GetIntensity(file))),
            }).ToList(),
            // The chromatographic peaks — the ONLY surface that fully represents match-between-runs.
            // FlashLFQ's peptide roll-up (the `peptides` list above, which mirrors
            // QuantifiedPeptides.tsv) reports far fewer MBR entries than were actually transferred:
            // on the K562 pair, 140 MBR peaks but only 52 appear as MBR at the peptide level, and a
            // whole run's transfers can vanish there. So a caller building an MBR-inclusive matrix
            // must read these peaks, not the peptide intensities. This mirrors QuantifiedPeaks.tsv.
            peaks = spectraFiles.SelectMany(file => results.Peaks[file].Select(peak => new
            {
                file_name = file.FilenameWithoutExtension,
                base_sequence = peak.Identifications.FirstOrDefault()?.BaseSequence ?? string.Empty,
                sequence = peak.Identifications.FirstOrDefault()?.ModifiedSequence ?? string.Empty,
                intensity = Finite(peak.Intensity),
                detection_type = peak.DetectionType.ToString(),
                retention_time = peak.Apex is null ? (double?)null : Finite(peak.ApexRetentionTime),
                // >1 identification means the peak is ambiguous — more than one peptide could explain
                // it — and its intensity should be treated with care.
                num_identifications = peak.Identifications.Count,
                protein_groups = string.Join(";", peak.Identifications
                    .SelectMany(id => id.ProteinGroups).Select(g => g.ProteinGroupName).Distinct()),
            })).ToList(),
            output_directory = string.IsNullOrWhiteSpace(outputDirectory) ? null : Path.GetFullPath(outputDirectory),
        };
    }

    /// <summary>A finite double, or null for NaN/±∞ — the only doubles System.Text.Json rejects.</summary>
    private static double? Finite(double value) => double.IsFinite(value) ? value : null;
}
