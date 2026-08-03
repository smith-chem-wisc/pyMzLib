using System.Globalization;
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

    /// <summary>
    /// <c>quant median-polish --peptides PATH [--shared-peptides] [--out DIR]</c> — roll a
    /// <c>QuantifiedPeptides.tsv</c> up to protein intensities with FlashLFQ's median-polish
    /// algorithm, without re-running peak-finding.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The heavy lifting is again entirely mzLib's. This verb reads a peptide table that FlashLFQ
    /// itself wrote (the <c>Intensity_</c> and <c>Detection Type_</c> columns keyed by run), rebuilds
    /// the FlashLFQ <see cref="Peptide"/> and <see cref="ProteinGroup"/> object graph those columns
    /// describe, and calls <see cref="FlashLfqResults.CalculateProteinResultsMedianPolish"/> — the
    /// same method the full <c>quant flashlfq</c> workflow runs. The bridge re-implements none of the
    /// algorithm; it only parses the table and shapes the result onto the wire.
    /// </para>
    /// <para>
    /// The experimental design comes on <b>stdin</b>, one run per line, tab-separated:
    /// <c>run_base_name[\tcondition[\tbiorep[\ttechrep[\tfraction]]]]</c>, where the run base name
    /// matches an <c>Intensity_&lt;name&gt;</c> column. Median polish groups measurements by condition
    /// and biological replicate, so the design is how a caller tells it which columns are replicates of
    /// which sample. With no design each column becomes its own biological replicate with a blank
    /// condition — exactly how FlashLFQ writes protein intensities when given no experimental-design
    /// file. Unlike the <c>flashlfq</c> verb these names are not files on disk, so nothing is opened.
    /// </para>
    /// </remarks>
    public static object MedianPolish(Program.Arguments arguments)
    {
        string peptidesPath = arguments.Required("peptides");
        bool useSharedPeptides = arguments.Flag("shared-peptides");
        string? outputDirectory = arguments.Optional("out");

        if (!File.Exists(peptidesPath))
            throw new Program.UsageException($"Quantified peptides file not found: '{peptidesPath}'.");

        Dictionary<string, DesignEntry> design = ParseMedianPolishDesign(ReadStdinLines());

        PeptideTable table = ReadPeptideTable(peptidesPath);
        List<SpectraFileInfo> spectraFiles = BuildDesignedFiles(table.RunNames, design);

        // Rebuild the FlashLFQ object graph the peptide table describes. The results constructor
        // populates PeptideModifiedSequences and ProteinGroups from a set of identifications, so one
        // synthetic identification per row carries the sequence and its protein group(s) across; the
        // measured numbers (intensity, detection type) are then written onto the reconstructed
        // peptides, since those normally arrive from peaks that this path never sees.
        var proteinGroupsByName = new Dictionary<string, ProteinGroup>(StringComparer.Ordinal);
        var identifications = new List<Identification>(table.Rows.Count);
        foreach (PeptideRow row in table.Rows)
        {
            List<ProteinGroup> groups = ResolveProteinGroups(row, proteinGroupsByName);
            identifications.Add(new Identification(
                spectraFiles[0], row.BaseSequence, row.Sequence,
                monoisotopicMass: 0, ms2RetentionTimeInMinutes: 0, chargeState: 0,
                proteinGroups: groups, useForProteinQuant: true));
        }

        var results = new FlashLfqResults(spectraFiles, identifications);

        foreach (PeptideRow row in table.Rows)
        {
            if (!results.PeptideModifiedSequences.TryGetValue(row.Sequence, out Peptide? peptide))
                continue;
            for (int f = 0; f < spectraFiles.Count; f++)
            {
                peptide.SetIntensity(spectraFiles[f], row.Intensities[f]);
                peptide.SetDetectionType(spectraFiles[f], row.DetectionTypes[f]);
            }
        }

        results.CalculateProteinResultsMedianPolish(useSharedPeptides);

        if (!string.IsNullOrWhiteSpace(outputDirectory))
        {
            Directory.CreateDirectory(outputDirectory);
            results.WriteResults(
                peaksOutputPath: null,
                modPeptideOutputPath: null,
                proteinOutputPath: Path.Combine(outputDirectory, "QuantifiedProteins.tsv"),
                bayesianProteinQuantOutput: null,
                silent: true);
        }

        return MedianPolishToWire(peptidesPath, useSharedPeptides, spectraFiles, results, outputDirectory);
    }

    /// <summary>One run's experimental-design coordinates, as parsed from the median-polish stdin.</summary>
    private readonly record struct DesignEntry(string Condition, int BiologicalReplicate, int TechnicalReplicate, int Fraction);

    /// <summary>A run column, its measured sequence rows, and where each run sits in the row arrays.</summary>
    private sealed record PeptideTable(List<string> RunNames, List<PeptideRow> Rows);

    /// <summary>One row of a <c>QuantifiedPeptides.tsv</c>, its per-run arrays aligned to <see cref="PeptideTable.RunNames"/>.</summary>
    private sealed record PeptideRow(
        string Sequence,
        string BaseSequence,
        List<string> ProteinGroupNames,
        List<string> GeneNames,
        List<string> Organisms,
        double[] Intensities,
        DetectionType[] DetectionTypes);

    /// <summary>
    /// Parses the optional experimental design from stdin. Keyed by run base name so it can be looked
    /// up against the peptide table's <c>Intensity_</c> columns. An empty design is legal and means
    /// "apply the default": the caller supplied no stdin at all.
    /// </summary>
    private static Dictionary<string, DesignEntry> ParseMedianPolishDesign(IReadOnlyList<string> lines)
    {
        var design = new Dictionary<string, DesignEntry>(StringComparer.Ordinal);

        for (int i = 0; i < lines.Count; i++)
        {
            string[] parts = lines[i].Split('\t');
            string name = parts[0].Trim();
            if (name.Length == 0)
                throw new Program.UsageException($"Design line {i + 1} has no run name.");
            if (design.ContainsKey(name))
                throw new Program.UsageException($"Design names run '{name}' more than once.");

            string condition = parts.Length > 1 ? parts[1].Trim() : string.Empty;
            // biorep falls back to the design's line order, so a bare list of runs makes each its own
            // biological replicate — the same rule the flashlfq verb applies to bare spectra paths.
            int biorep = ParseDesignField(parts, 2, i, "biorep", fallback: i);
            int techrep = ParseDesignField(parts, 3, i, "techrep", fallback: 0);
            int fraction = ParseDesignField(parts, 4, i, "fraction", fallback: 0);

            design[name] = new DesignEntry(condition, biorep, techrep, fraction);
        }

        return design;
    }

    /// <summary>
    /// Builds the <see cref="SpectraFileInfo"/> per run column, applying the design where one was
    /// given and the default otherwise. When a design is present it must name every run in the table
    /// and no run absent from it — a name that matches no column is a typo that would otherwise
    /// silently do nothing, and a column with no design entry would quantify against an ambiguous
    /// replicate assignment.
    /// </summary>
    private static List<SpectraFileInfo> BuildDesignedFiles(List<string> runNames, Dictionary<string, DesignEntry> design)
    {
        if (design.Count > 0)
        {
            var runSet = new HashSet<string>(runNames, StringComparer.Ordinal);
            List<string> unmatched = design.Keys.Where(name => !runSet.Contains(name))
                .OrderBy(name => name, StringComparer.Ordinal).ToList();
            if (unmatched.Count > 0)
                throw new Program.UsageException(
                    $"The design names {unmatched.Count} run(s) with no matching Intensity_ column: " +
                    $"{string.Join(", ", unmatched.Take(5))}{(unmatched.Count > 5 ? ", …" : string.Empty)}. " +
                    "Design run names must match the peptide table's run columns.");

            List<string> undesigned = runNames.Where(name => !design.ContainsKey(name)).ToList();
            if (undesigned.Count > 0)
                throw new Program.UsageException(
                    $"The peptide table has {undesigned.Count} run(s) the design does not mention: " +
                    $"{string.Join(", ", undesigned.Take(5))}{(undesigned.Count > 5 ? ", …" : string.Empty)}. " +
                    "Give every run a design line, or supply no design at all to default each to its own replicate.");
        }

        var files = new List<SpectraFileInfo>(runNames.Count);
        for (int i = 0; i < runNames.Count; i++)
        {
            string name = runNames[i];
            DesignEntry entry = design.Count > 0
                ? design[name]
                : new DesignEntry(string.Empty, BiologicalReplicate: i, TechnicalReplicate: 0, Fraction: 0);
            // The run name is not a path — nothing is opened — so it stands in for the file path; the
            // base name recovered from it is what every roll-up and label keys on.
            files.Add(new SpectraFileInfo(name, entry.Condition,
                biorep: entry.BiologicalReplicate, techrep: entry.TechnicalReplicate, fraction: entry.Fraction));
        }

        return files;
    }

    /// <summary>
    /// Reads a <c>QuantifiedPeptides.tsv</c> into run names and rows. The columns are located by their
    /// header names — the layout FlashLFQ's <see cref="Peptide.TabSeparatedHeader"/> writes — so an
    /// added IsoTracker retention-time block, or a reordering, does not shift the parse.
    /// </summary>
    private static PeptideTable ReadPeptideTable(string path)
    {
        string[] lines = File.ReadAllLines(path);
        if (lines.Length == 0 || string.IsNullOrWhiteSpace(lines[0]))
            throw new Program.UsageException("The quantified peptides file is empty.");

        string[] header = lines[0].Split('\t');
        int seqIndex = RequireColumn(header, "Sequence");
        int baseIndex = RequireColumn(header, "Base Sequence");
        int proteinIndex = RequireColumn(header, "Protein Groups");
        int geneIndex = FindColumn(header, "Gene Names", "Gene Name");
        int organismIndex = FindColumn(header, "Organism");

        // Run columns are recovered from the Intensity_ headers, in file order; the matching
        // Detection Type_ column for each run is found by name so the two need not be adjacent.
        var runNames = new List<string>();
        var intensityColumns = new List<int>();
        var detectionColumnByRun = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int c = 0; c < header.Length; c++)
        {
            string name = header[c];
            if (name.StartsWith("Intensity_", StringComparison.Ordinal))
            {
                runNames.Add(name.Substring("Intensity_".Length));
                intensityColumns.Add(c);
            }
            else if (name.StartsWith("Detection Type_", StringComparison.Ordinal))
            {
                detectionColumnByRun[name.Substring("Detection Type_".Length)] = c;
            }
        }

        if (runNames.Count == 0)
            throw new Program.UsageException(
                "The quantified peptides file has no Intensity_ columns; it does not look like a FlashLFQ QuantifiedPeptides.tsv.");

        var rows = new List<PeptideRow>(lines.Length - 1);
        // Sequence keys the peptide graph, so two rows sharing one would have the second silently
        // overwrite the first's measurements. Rejecting says so instead of quantifying from half a
        // table.
        var seenSequences = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int r = 1; r < lines.Length; r++)
        {
            if (string.IsNullOrWhiteSpace(lines[r]))
                continue;
            string[] fields = lines[r].Split('\t');

            string sequence = Field(fields, seqIndex);
            if (sequence.Length == 0)
                continue;

            if (seenSequences.TryGetValue(sequence, out int firstSeenOn))
                throw new Program.UsageException(
                    $"Line {r + 1} repeats the sequence '{sequence}', first seen on line {firstSeenOn}. " +
                    "Each sequence must appear once; a FlashLFQ QuantifiedPeptides.tsv has one row per peptide.");
            seenSequences[sequence] = r + 1;

            var intensities = new double[runNames.Count];
            var detections = new DetectionType[runNames.Count];
            for (int f = 0; f < runNames.Count; f++)
            {
                intensities[f] = ParseIntensity(Field(fields, intensityColumns[f]), r + 1, "Intensity_" + runNames[f]);
                bool detectionColumnPresent = detectionColumnByRun.TryGetValue(runNames[f], out int dc);
                detections[f] = ResolveDetectionType(
                    detectionColumnPresent ? Field(fields, dc) : string.Empty,
                    detectionColumnPresent,
                    intensities[f],
                    r + 1,
                    "Detection Type_" + runNames[f]);
            }

            rows.Add(new PeptideRow(
                sequence,
                Field(fields, baseIndex),
                SplitGroups(Field(fields, proteinIndex)),
                SplitGroups(Field(fields, geneIndex)),
                SplitGroups(Field(fields, organismIndex)),
                intensities,
                detections));
        }

        return new PeptideTable(runNames, rows);
    }

    /// <summary>
    /// Turns a row's protein-group column into interned <see cref="ProteinGroup"/> objects. FlashLFQ
    /// <c>;</c>-joins the names (and the parallel gene/organism columns) when a peptide is shared, so
    /// the cell is split back into one group per name; gene and organism are paired positionally when
    /// their counts line up and otherwise carried whole. Interning by name keeps a single instance per
    /// group so peptide membership and the results' group map agree.
    /// </summary>
    private static List<ProteinGroup> ResolveProteinGroups(PeptideRow row, Dictionary<string, ProteinGroup> interned)
    {
        var groups = new List<ProteinGroup>(row.ProteinGroupNames.Count);
        bool genesAlign = row.GeneNames.Count == row.ProteinGroupNames.Count;
        bool organismsAlign = row.Organisms.Count == row.ProteinGroupNames.Count;

        for (int i = 0; i < row.ProteinGroupNames.Count; i++)
        {
            string name = row.ProteinGroupNames[i];
            if (!interned.TryGetValue(name, out ProteinGroup? group))
            {
                string gene = genesAlign ? row.GeneNames[i] : string.Join(";", row.GeneNames);
                string organism = organismsAlign ? row.Organisms[i] : string.Join(";", row.Organisms);
                group = new ProteinGroup(name, gene, organism);
                interned[name] = group;
            }
            groups.Add(group);
        }

        return groups;
    }

    /// <summary>Flattens the median-polish protein results onto the wire.</summary>
    private static object MedianPolishToWire(
        string peptidesPath,
        bool useSharedPeptides,
        List<SpectraFileInfo> spectraFiles,
        FlashLfqResults results,
        string? outputDirectory)
    {
        // Group runs into samples the way FlashLFQ's protein quant does — by condition, then
        // biological replicate. Ordered the same way the engine orders them
        // (CalculateProteinResultsMedianPolish walks GroupBy(Condition).OrderBy(Key) then
        // GroupBy(BiologicalReplicate).OrderBy(Key)) so that samples[i] here is the engine's sample i.
        // Intensities below are read per file rather than by index, so the order is not load-bearing
        // today; it is matched so it cannot quietly become wrong.
        //
        // A run's label is its own name when there is no design to label it with, and
        // "condition_biorep" once a real design groups runs.
        //
        // NOTE: this deliberately leads mzLib. FlashLFQ's own QuantifiedProteins.tsv applies the same
        // rule with the boolean inverted, so an unfractionated run is labelled by file name exactly
        // when a design exists and yields "Intensity__1" when one does not — smith-chem-wisc/mzLib#1128,
        // fixed by mzLib#1129. Until that lands and the pin moves, the labels here and the ones in a
        // file written by --out disagree for unfractionated data. The values do not: the engine sets a
        // sample's intensity on its first run and zeroes the rest, so both readings agree.
        bool unfractionated = spectraFiles.Select(f => f.Fraction).Distinct().Count() == 1;
        bool conditionsUndefined = spectraFiles.All(f => f.Condition == "Default")
            || spectraFiles.All(f => string.IsNullOrWhiteSpace(f.Condition));
        bool labelByRunName = conditionsUndefined && unfractionated;

        var samples = new List<(string Label, IGrouping<int, SpectraFileInfo> Files)>();
        foreach (var conditionGroup in spectraFiles.GroupBy(f => f.Condition).OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            foreach (var replicate in conditionGroup.GroupBy(f => f.BiologicalReplicate).OrderBy(g => g.Key))
            {
                SpectraFileInfo representative = replicate.First();
                // The run name was handed to SpectraFileInfo as its path, so it survives whole in
                // FullFilePathWithExtension. FilenameWithoutExtension would truncate it at the last
                // dot — "QC.2" becomes "QC" — which both mislabels the sample and can collide two
                // runs onto one key.
                string label = labelByRunName
                    ? representative.FullFilePathWithExtension
                    : representative.Condition + "_" + (representative.BiologicalReplicate + 1);
                samples.Add((label, replicate));
            }
        }

        return new
        {
            peptides_file = Path.GetFullPath(peptidesPath),
            parameters = new
            {
                use_shared_peptides_for_protein_quant = useSharedPeptides,
            },
            samples = samples.Select(s => new
            {
                label = s.Label,
                condition = s.Files.First().Condition,
                biological_replicate = s.Files.First().BiologicalReplicate,
            }).ToList(),
            peptide_count = results.PeptideModifiedSequences.Count,
            protein_count = results.ProteinGroups.Count,
            proteins = results.ProteinGroups.Values
                .OrderBy(g => g.ProteinGroupName, StringComparer.Ordinal)
                .Select(group => new
                {
                    protein_group = group.ProteinGroupName,
                    gene_name = group.GeneName,
                    organism = group.Organism,
                    // One run per sample: read the intensity straight off the representative. A sample
                    // aggregating several runs (fractions, technical replicates): sum them, as
                    // ProteinGroup.ToString does — the engine sets the value on the first run only, so
                    // the sum recovers it. NaN — FlashLFQ's "unquantifiable" — crosses as null.
                    intensities = samples.ToDictionary(
                        s => s.Label,
                        s => Finite(labelByRunName
                            ? group.GetIntensity(s.Files.First())
                            : s.Files.Sum(f => group.GetIntensity(f)))),
                }).ToList(),
            output_directory = string.IsNullOrWhiteSpace(outputDirectory) ? null : Path.GetFullPath(outputDirectory),
        };
    }

    private static int RequireColumn(string[] header, string name)
    {
        int index = FindColumn(header, name);
        if (index < 0)
            throw new Program.UsageException(
                $"The quantified peptides file is missing the required '{name}' column.");
        return index;
    }

    private static int FindColumn(string[] header, params string[] names)
    {
        foreach (string name in names)
        {
            int index = Array.FindIndex(header, h => string.Equals(h, name, StringComparison.Ordinal));
            if (index >= 0)
                return index;
        }
        return -1;
    }

    private static string Field(string[] fields, int index) =>
        index >= 0 && index < fields.Length ? fields[index].Trim() : string.Empty;

    /// <summary>Splits a <c>;</c>-joined cell into its parts, dropping blanks. An empty cell yields no parts.</summary>
    private static List<string> SplitGroups(string cell) =>
        cell.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries).ToList();

    /// <summary>
    /// Parses an intensity cell. A blank cell is no measurement (0) — FlashLFQ writes those. A cell
    /// that is present but unreadable is rejected: silently reading it as 0 would turn a corrupt or
    /// mis-delimited table into a table of "not measured", which is indistinguishable from a real
    /// result and would quantify proteins from data that was never there.
    /// </summary>
    private static double ParseIntensity(string cell, int lineNumber, string column)
    {
        if (cell.Length == 0)
            return 0.0;
        if (double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
            return value;

        throw new Program.UsageException(
            $"Line {lineNumber}, column '{column}': '{cell}' is not a number. " +
            "Intensity cells must be numeric (a blank cell means the peptide was not measured in that run).");
    }

    /// <summary>
    /// Resolves a run's detection type from its cell.
    /// </summary>
    /// <remarks>
    /// <para>
    /// When the table carries no <c>Detection Type_</c> column for the run at all, the type is
    /// inferred from the intensity — a measured run is treated as an MS/MS detection rather than
    /// silently excluding the peptide from protein quant, since
    /// <see cref="Peptide.UnambiguousPeptideQuant"/> requires at least one non-ambiguous detection.
    /// That is the one inference this verb makes, and it is made only where the table says nothing.
    /// </para>
    /// <para>
    /// A column that is present but holds an unrecognized value is rejected instead. Inferring there
    /// would overwrite what the table actually said with a guess, and an unreadable detection type is
    /// a sign the table is not the FlashLFQ output it claims to be.
    /// </para>
    /// </remarks>
    private static DetectionType ResolveDetectionType(
        string cell, bool columnPresent, double intensity, int lineNumber, string column)
    {
        if (!columnPresent || cell.Length == 0)
            return intensity > 0 ? DetectionType.MSMS : DetectionType.NotDetected;

        if (Enum.TryParse(cell.Trim(), out DetectionType parsed))
            return parsed;

        throw new Program.UsageException(
            $"Line {lineNumber}, column '{column}': '{cell}' is not a FlashLFQ detection type. " +
            $"Expected one of {string.Join(", ", Enum.GetNames<DetectionType>())}.");
    }

    /// <summary>Reads stdin into non-blank trimmed lines.</summary>
    /// <remarks>
    /// A UTF-8 byte-order mark is stripped from the first line if present. A caller that pipes a
    /// BOM-prefixed stream (some shells and editors add one) would otherwise carry the mark into the
    /// first field — a file path or run name that then matches nothing — so it is removed here once,
    /// where every stdin-consuming verb benefits.
    /// </remarks>
    private static List<string> ReadStdinLines()
    {
        var lines = new List<string>();
        string? line;
        bool first = true;
        while ((line = Console.In.ReadLine()) != null)
        {
            if (first)
            {
                line = line.TrimStart('﻿');
                first = false;
            }
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
