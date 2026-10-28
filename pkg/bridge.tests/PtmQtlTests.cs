using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.Json;
using MzLibUtil;
using Quantification.PtmQtl;
using StatisticalModels;

namespace MzLibBridge.Tests;

/// <summary>
/// The mzLib cases the <c>ptmqtl</c> fixtures are made of: mzLib's own hand-computed PtmQtl test cases
/// (<c>mzLib/Test/Quantification/PtmQtlTests.cs</c> at the pin), each on its own protein so that no case
/// changes another (occupancy and physical pairs are computed per protein). The committed input tables
/// under <c>pkg/python/tests/fixtures/ptmqtl/</c> are rendered from these lists and held to them by
/// <see cref="PtmQtlTests.Fixtures_AreMzLibsHandComputedCases"/>.
/// </summary>
internal static class PtmQtlCases
{
    internal const string Phos = "Common Biological:Phosphorylation on S";
    internal const string HexNAc = "Common Biological:HexNAc on S";
    internal const string Acet = "Common Biological:Acetylation on K";
    internal const string Gg = "Trypsin Digested:GG (Ubiquitination Site) on K";

    internal static PeptidoformObservation Obs(string run, string full, int start, double intensity, string protein,
        char? previous = null) =>
        new(run, full, protein, start, start + full.GetBaseSequenceFromFullSequence().Length - 1, intensity, previous);

    /// <summary>
    /// A phosphosite on <c>&lt;prefix&gt;SK</c> starting at <paramref name="start"/> (so at start + 3), as a modified
    /// and an unmodified form. Each protein of a co-varying case gets its own prefix: mzLib calls two sites whose
    /// covering peptides share a base sequence overlapping, and leaves the pair out of the BH family.
    /// </summary>
    private static void Site(List<PeptidoformObservation> list, string run, string protein, int start, double modified,
        double covering, string prefix = "PEP")
    {
        list.Add(Obs(run, $"{prefix}S[{Phos}]K", start, modified, protein));
        list.Add(Obs(run, $"{prefix}SK", start, covering - modified, protein));
    }

    private static readonly string[] Prefixes = { "PEP", "GGA", "ALA", "VVA", "TTA", "MMA", "QQA", "NNA", "DDA", "EEA", "FFA", "HHA", "IIA" };

    /// <summary>The single-run cases of mzLib's PtmQtlTests, one protein each, with mzLib's own run names.</summary>
    internal static List<PeptidoformObservation> Cases() => new()
    {
        // OccupancyIsModifiedOverCoveringIntensity: 30 / (30 + 70 + 100) = 0.15.
        Obs("r1", $"PEPS[{Phos}]K", 10, 30, "RATIO"),
        Obs("r1", "PEPSK", 10, 70, "RATIO"),
        Obs("r1", "AAPEPSKR", 8, 100, "RATIO"),
        // SeenOnlyModifiedIsOccupancyOneWithTheCeilingFlag.
        Obs("r1", $"PEPS[{Phos}]K", 10, 30, "CEILING"),
        // TheFourStatesAreKeptApart.
        Obs("floor", $"PEPS[{Phos}]K", 10, double.NaN, "STATES"),
        Obs("floor", "PEPSK", 10, 50, "STATES"),
        Obs("countOnly", $"PEPS[{Phos}]K", 10, double.NaN, "STATES"),
        Obs("notDetected", "PEPSK", 10, 40, "STATES"),
        Obs("quant", $"PEPS[{Phos}]K", 10, 5, "STATES"),
        Obs("quant", "PEPSK", 10, 45, "STATES"),
        // TwoModificationsAtOnePositionShareTheDenominator: 20 and 30 of 100.
        Obs("r1", $"PEPS[{Phos}]K", 10, 20, "SHARED_POSITION"),
        Obs("r1", $"PEPS[{HexNAc}]K", 10, 30, "SHARED_POSITION"),
        Obs("r1", "PEPSK", 10, 50, "SHARED_POSITION"),
        // ExclusionsAndFilterFollowDefOccPsms and PeptideNTerminalModificationsAndNonMetStartsAreNotProteinNTermini.
        Obs("r1", "PEPM[Common Variable:Oxidation on M]C[Common Fixed:Carbamidomethyl on C]K", 10, 10, "EXCLUDED"),
        Obs("r1", "[Common Artifact:Ammonia loss on C]CPEPK", 1, 10, "PEPTIDE_NTERM"),
        Obs("r1", "[UniProt:N-acetylserine on S]SPEPK", 2, 10, "NTERM_AFTER_MET", 'M'),
        Obs("r1", "[UniProt:N-acetylserine on S]SPEPK", 2, 10, "NTERM_AFTER_K", 'K'),
        Obs("r1", "[Renamed:Unregistered on C]CPEPK", 1, 10, "UNREGISTERED"),
        // PhysicalPairsCarryCoOccupancy: 25 of 100 covering both, identified together in two runs.
        Obs("r1", $"PEPS[{Phos}]K[{Acet}]R", 10, 25, "TWO_ON_ONE"),
        Obs("r1", $"PEPS[{Phos}]KR", 10, 25, "TWO_ON_ONE"),
        Obs("r1", "PEPSKR", 10, 50, "TWO_ON_ONE"),
        Obs("r2", $"PEPS[{Phos}]K[{Acet}]R", 10, double.NaN, "TWO_ON_ONE"),
    };

    /// <summary>FractionsOfOneSampleAreSummedBeforeTheRatio and AFormUnquantifiedInOneFractionAndCoveredInAnotherIsAFloorForTheSample.</summary>
    internal static List<PeptidoformObservation> Fractions() => new()
    {
        Obs("f1", $"PEPS[{Phos}]K", 10, 30, "SUMMED"),
        Obs("f1", "PEPSK", 10, 70, "SUMMED"),
        Obs("f2", $"PEPS[{Phos}]K", 10, 10, "SUMMED"),
        Obs("f2", "PEPSK", 10, 90, "SUMMED"),
        Obs("g1", "PEPSK", 10, 40, "SUMMED"),
        Obs("f1", $"PEPS[{Phos}]K", 10, double.NaN, "FLOOR_ACROSS"),
        Obs("f2", "PEPSK", 10, 50, "FLOOR_ACROSS"),
    };

    internal static readonly (string Run, string Sample)[] FractionSamples = { ("f1", "A"), ("f2", "A"), ("g1", "B") };

    /// <summary>
    /// RunLevelRemovalTakesOutABatchShift: twelve runs in two batches, every site 0.2 higher in the second.
    /// Occupancy in percent of a covering intensity of 100, so every intensity is exact.
    /// </summary>
    internal static List<PeptidoformObservation> Batch()
    {
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 12; r++)
        {
            string run = $"r{r:D2}";
            double batch = r < 6 ? 10 : 30;
            Site(list, run, "PAIR_S", 10, batch + noiseS[r % 6], 100, Prefixes[0]);
            Site(list, run, "PAIR_T", 40, batch + noiseT[r % 6], 100, Prefixes[1]);
            for (int i = 1; i <= 4; i++)
                Site(list, run, $"LEVEL_{i}", 4, batch + 3.5, 100, Prefixes[1 + i]);
        }
        return list;
    }

    /// <summary>
    /// GroupedRunLevelTakesOutAShiftInOneChemistry: GG sites shift with a batch, phosphosites do not.
    /// Occupancy in thousandths of a covering intensity of 1000.
    /// </summary>
    internal static List<PeptidoformObservation> Grouped()
    {
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        double[] phosphoLevel = { 300, 100, 250, 150, 200, 350, 120, 280, 180, 330, 220, 140 };
        var list = new List<PeptidoformObservation>();
        void Gg(string run, string protein, double modified, string prefix)
        {
            list.Add(Obs(run, $"{prefix}K[{PtmQtlCases.Gg}]R", 10, modified, protein));
            list.Add(Obs(run, $"{prefix}KR", 10, 1000 - modified, protein));
        }
        for (int r = 0; r < 12; r++)
        {
            string run = $"r{r:D2}";
            double batch = r < 6 ? 100 : 300;
            Gg(run, "GG_S", batch + 10 * noiseS[r % 6], Prefixes[0]);
            Gg(run, "GG_T", batch + 10 * noiseT[r % 6], Prefixes[1]);
            for (int i = 1; i <= 4; i++)
                Gg(run, $"GG_LEVEL_{i}", batch + 35, Prefixes[1 + i]);
            for (int i = 0; i < 7; i++)
                Site(list, run, $"PHOSPHO_{i}", 6, phosphoLevel[r] + (i == 0 ? 10 * noiseS[r % 6] : i), 1000, Prefixes[6 + i]);
        }
        return list;
    }

    /// <summary>CeilingCellsStayOutOfTypeAByDefault: opposed in five runs, both seen only modified in three.</summary>
    internal static List<PeptidoformObservation> Ceiling()
    {
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 5; r++)
        {
            Site(list, $"r{r}", "CEIL_S", 10, 10 * (r + 1), 100, Prefixes[0]);
            Site(list, $"r{r}", "CEIL_T", 40, 10 * (5 - r), 100, Prefixes[1]);
        }
        for (int r = 5; r < 8; r++)
        {
            list.Add(Obs($"r{r}", $"{Prefixes[0]}S[{Phos}]K", 10, 100, "CEIL_S"));
            list.Add(Obs($"r{r}", $"{Prefixes[1]}S[{Phos}]K", 40, 100, "CEIL_T"));
        }
        return list;
    }

    private static readonly double[] TraitLevels = { 0.25, 0.25, 1, 1, 2, 2 };
    private static readonly double[] TraitOffsets = { 0.1, -0.1, 0.05, -0.05, 0.08, -0.08 };

    private static double Expit(double x) => 1 / (1 + Math.Exp(-x));

    /// <summary>mzLib's TraitDesign: logit = -2 + 0.5 trait + replicate offset, two injections of six replicates.</summary>
    internal static double TraitFraction(int r, int injection) =>
        Expit(-2 + 0.5 * TraitLevels[r] + TraitOffsets[r] + (injection == 0 ? 0.01 : -0.01) + (injection == 1 ? 0.003 * (r % 3 - 1) : 0));

    /// <summary>
    /// TRAIT_SITE is TraitEffectIsRecoveredOnTheLogitScale; CENSORED_SITE is TraitFitUsesOnlyQuantifiedNonCeilingCells,
    /// whose replicate R5 is a floor in one injection and a ceiling in the other. Covering intensity 1000.
    /// </summary>
    internal static List<PeptidoformObservation> Trait()
    {
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 6; r++)
            for (int injection = 0; injection < 2; injection++)
            {
                string run = $"R{r}_{injection}";
                double modified = 1000 * TraitFraction(r, injection);
                Site(list, run, "TRAIT_SITE", 10, modified, 1000, Prefixes[0]);
                if (r == 5 && injection == 0)
                {
                    list.Add(Obs(run, $"{Prefixes[1]}S[{Phos}]K", 10, double.NaN, "CENSORED_SITE"));
                    list.Add(Obs(run, $"{Prefixes[1]}SK", 10, 1000 - modified, "CENSORED_SITE"));
                }
                else if (r == 5)
                    list.Add(Obs(run, $"{Prefixes[1]}S[{Phos}]K", 10, modified, "CENSORED_SITE"));
                else
                    Site(list, run, "CENSORED_SITE", 10, modified, 1000, Prefixes[1]);
            }
        return list;
    }

    /// <summary>
    /// BoundsEnterTypeAAsTiedRanks: two sites with a value in six of ten runs, and in the other four a
    /// covered-zero, a floor, a ceiling, or a count-only cell seen only modified (a bound) or in both forms (no
    /// rank). Occupancy in percent of a covering intensity of 100.
    /// </summary>
    internal static List<PeptidoformObservation> Bounds()
    {
        double[] valuesS = { 10, 20, 30, 40, 50, 60 }, valuesT = { 15, 10, 35, 30, 60, 45 };
        string s = Prefixes[0], t = Prefixes[1];
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 6; r++)
        {
            Site(list, $"r{r}", "BOUND_S", 10, valuesS[r], 100, s);
            Site(list, $"r{r}", "BOUND_T", 40, valuesT[r], 100, t);
        }
        // r6: s a covered-zero (only the unmodified form identified); t a ceiling (only the modified form quantified).
        list.Add(Obs("r6", $"{s}SK", 10, 100, "BOUND_S"));
        list.Add(Obs("r6", $"{t}S[{Phos}]K", 40, 100, "BOUND_T"));
        // r7: s a floor (modified identified, unquantified; unmodified quantified); t a covered-zero.
        list.Add(Obs("r7", $"{s}S[{Phos}]K", 10, double.NaN, "BOUND_S"));
        list.Add(Obs("r7", $"{s}SK", 10, 100, "BOUND_S"));
        list.Add(Obs("r7", $"{t}SK", 40, 100, "BOUND_T"));
        // r8: s a ceiling; t count-only with both forms identified, nothing quantified: no rank.
        list.Add(Obs("r8", $"{s}S[{Phos}]K", 10, 100, "BOUND_S"));
        list.Add(Obs("r8", $"{t}S[{Phos}]K", 40, double.NaN, "BOUND_T"));
        list.Add(Obs("r8", $"{t}SK", 40, double.NaN, "BOUND_T"));
        // r9: s count-only seen only modified (a high bound); t a floor.
        list.Add(Obs("r9", $"{s}S[{Phos}]K", 10, double.NaN, "BOUND_S"));
        list.Add(Obs("r9", $"{t}S[{Phos}]K", 40, double.NaN, "BOUND_T"));
        list.Add(Obs("r9", $"{t}SK", 40, 100, "BOUND_T"));
        return list;
    }

    /// <summary>
    /// APairWhoseSitesAreBothMostlyBoundsIsNotTested: twenty runs; s and t have a value in the same two runs and
    /// are covered-zeros elsewhere, u varies everywhere. Occupancy in percent of a covering intensity of 100.
    /// </summary>
    internal static List<PeptidoformObservation> MostlyBounds()
    {
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 20; r++)
        {
            string run = $"r{r:D2}";
            if (r < 2)
            {
                Site(list, run, "MOSTLY_S", 10, 10 + 10 * r, 100, Prefixes[0]);
                Site(list, run, "MOSTLY_T", 40, 20 + 10 * r, 100, Prefixes[1]);
            }
            else
            {
                list.Add(Obs(run, $"{Prefixes[0]}SK", 10, 100, "MOSTLY_S"));
                list.Add(Obs(run, $"{Prefixes[1]}SK", 40, 100, "MOSTLY_T"));
            }
            Site(list, run, "MOSTLY_U", 4, (7 * r) % 20 + 1, 100, Prefixes[2]);
        }
        return list;
    }

    /// <summary>The target PSMs per run of ARunCovariateIsPartialledOut: the run's depth.</summary>
    internal static readonly double[] DepthPsms = { 9000, 12000, 8000, 15000, 11000, 7000, 16000, 10000, 13000, 9500, 14000, 8500 };

    /// <summary>
    /// ARunCovariateIsPartialledOut: twelve runs in two batches, two sites that rise with the run's depth (target
    /// PSMs), and four background sites at the batch's centre. Occupancy in millionths of a covering intensity of
    /// 1e6, so every intensity is an exact integer.
    /// </summary>
    internal static List<PeptidoformObservation> Depth()
    {
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        var list = new List<PeptidoformObservation>();
        for (int r = 0; r < 12; r++)
        {
            string run = $"r{r:D2}";
            double batch = r < 6 ? 100_000 : 300_000;
            Site(list, run, "DEPTH_S", 10, batch + 10_000 * noiseS[r % 6] + DepthPsms[r], 1_000_000, Prefixes[0]);
            Site(list, run, "DEPTH_T", 40, batch + 10_000 * noiseT[r % 6] + 2 * DepthPsms[r], 1_000_000, Prefixes[1]);
            for (int i = 1; i <= 4; i++)
                Site(list, run, $"DEPTH_LEVEL_{i}", 4, batch + 35_000, 1_000_000, Prefixes[1 + i]);
        }
        return list;
    }

    internal static Dictionary<string, double> DepthMap() =>
        Enumerable.Range(0, 12).ToDictionary(r => $"r{r:D2}", r => DepthPsms[r]);

    internal static string DepthPsmsTsv() =>
        "run\ttarget_psms\n" + string.Concat(Enumerable.Range(0, 12).Select(r => $"r{r:D2}\t{R(DepthPsms[r])}\n"));

    /// <summary>Run, replicate, the trait (dose) and the injection covariate of TraitFitReportsCovariateEffects.</summary>
    internal static IEnumerable<(string Run, string Replicate, double Dose, double Injection)> TraitRuns() =>
        from r in Enumerable.Range(0, 6)
        from injection in new[] { 0, 1 }
        select ($"R{r}_{injection}", $"R{r}", TraitLevels[r], injection == 0 ? 1.0 : 0.0);

    internal static string R(double value) => value.ToString("R", CultureInfo.InvariantCulture);

    internal static string ObservationsTsv(IEnumerable<PeptidoformObservation> observations)
    {
        var text = new StringBuilder("run\tfull_sequence\tprotein_accession\tstart_residue\tend_residue\tintensity\tprevious_residue\n");
        foreach (var o in observations)
            text.Append($"{o.Run}\t{o.FullSequence}\t{o.ProteinAccession}\t{o.StartResidue}\t{o.EndResidue}\t" +
                        $"{(double.IsNaN(o.Intensity) ? "" : R(o.Intensity))}\t{o.PreviousResidue}\n");
        return text.ToString();
    }

    internal static string SamplesTsv() =>
        "run\tsample\n" + string.Concat(FractionSamples.Select(s => $"{s.Run}\t{s.Sample}\n"));

    internal static string TraitRunsTsv(bool covariate = true) =>
        (covariate ? "run\treplicate\tdose\tinjection\n" : "run\treplicate\tdose\n") +
        string.Concat(TraitRuns().Select(t => $"{t.Run}\t{t.Replicate}\t{R(t.Dose)}" + (covariate ? $"\t{R(t.Injection)}" : "") + "\n"));
}

/// <summary>
/// The <c>ptmqtl</c> verbs. Their numbers must be mzLib's <c>Quantification.PtmQtl</c>: each test runs a
/// verb on a table, makes the same mzLib call on the same records in memory, and holds every column of
/// the wire to it. Where mzLib's own test asserts a hand-computed value (30 / 200, a co-occupancy of
/// 0.25, a trait effect near 0.5), the wire is held to that value too.
/// </summary>
[TestFixture]
public class PtmQtlTests
{
    private string _temp = "";

    [SetUp]
    public void CreateTemp()
    {
        _temp = Path.Combine(Path.GetTempPath(), "mzlib-bridge-ptmqtl-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_temp);
    }

    [TearDown]
    public void DeleteTemp()
    {
        if (Directory.Exists(_temp))
            Directory.Delete(_temp, recursive: true);
    }

    // ---- helpers ---------------------------------------------------------------------------------

    private static string RepoRoot([CallerFilePath] string thisFile = "") =>
        Path.GetFullPath(Path.Combine(Path.GetDirectoryName(thisFile)!, "..", ".."));

    internal static string Fixtures => Path.Combine(RepoRoot(), "pkg", "python", "tests", "fixtures", "ptmqtl");

    private string Write(string name, string content)
    {
        string path = Path.Combine(_temp, name);
        File.WriteAllText(path, content);
        return path;
    }

    private string WriteObservations(IEnumerable<PeptidoformObservation> observations, string name = "observations.tsv") =>
        Write(name, PtmQtlCases.ObservationsTsv(observations));

    private static JsonElement Run(Func<Program.Arguments, object> verb, params string[] args) =>
        JsonDocument.Parse(JsonSerializer.Serialize(verb(new Program.Arguments(args)), Program.JsonOptions)).RootElement;

    private static string Usage(Func<Program.Arguments, object> verb, params string[] args) =>
        Assert.Throws<Program.UsageException>(() => verb(new Program.Arguments(args)))!.Message;

    private static List<JsonElement> Column(JsonElement root, string name) =>
        root.GetProperty("columns").GetProperty(name).EnumerateArray().ToList();

    private static double? D(JsonElement e) => e.ValueKind == JsonValueKind.Null ? null : e.GetDouble();

    private static string? S(JsonElement e) => e.ValueKind == JsonValueKind.Null ? null : e.GetString();

    private static int? I(JsonElement e) => e.ValueKind == JsonValueKind.Null ? null : e.GetInt32();

    private static string ReplaceLast(string text, string what, string with)
    {
        int at = text.LastIndexOf(what, StringComparison.Ordinal);
        return text[..at] + with + text[(at + what.Length)..];
    }

    private static List<string> Caveats(JsonElement root) =>
        root.GetProperty("caveats").EnumerateArray().Select(e => e.GetString()!).ToList();

    private static double? Finite(double value) => double.IsFinite(value) ? value : null;

    /// <summary>The row of the table whose site_key is <paramref name="key"/> and run is <paramref name="run"/>.</summary>
    private static int RowOf(JsonElement root, string key, string run)
    {
        var keys = Column(root, "site_key");
        var runs = Column(root, "run");
        var rows = Enumerable.Range(0, keys.Count).Where(i => keys[i].GetString() == key && runs[i].GetString() == run).ToList();
        Assert.That(rows, Has.Count.EqualTo(1), $"{key} in {run}");
        return rows[0];
    }

    /// <summary>Every column of an occupancy table equals mzLib's cells, in mzLib's order.</summary>
    private static void AssertOccupancyIs(JsonElement root, IReadOnlyList<SiteRunOccupancy> expected)
    {
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));
        Assert.That(root.GetProperty("column_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(PtmQtl.OccupancyColumns));
        for (int i = 0; i < expected.Count; i++)
        {
            SiteRunOccupancy c = expected[i];
            string at = $"row {i} ({c.Site.Key} in {c.Run})";
            Assert.That(Column(root, "site_key")[i].GetString(), Is.EqualTo(c.Site.Key), at);
            Assert.That(Column(root, "protein_accession")[i].GetString(), Is.EqualTo(c.Site.ProteinAccession), at);
            Assert.That(Column(root, "position")[i].GetInt32(), Is.EqualTo(c.Site.Position), at);
            Assert.That(Column(root, "residue")[i].GetString(), Is.EqualTo(c.Site.Residue.ToString()), at);
            Assert.That(Column(root, "modification")[i].GetString(), Is.EqualTo(c.Site.Modification), at);
            Assert.That(Column(root, "category")[i].GetString(), Is.EqualTo(c.Site.Category), at);
            Assert.That(Column(root, "run")[i].GetString(), Is.EqualTo(c.Run), at);
            Assert.That(Column(root, "state")[i].GetString(), Is.EqualTo(PtmQtl.StateName(c.State)), at);
            Assert.That(Column(root, "modified_intensity")[i].GetDouble(), Is.EqualTo(c.ModifiedIntensity), at);
            Assert.That(Column(root, "covering_intensity")[i].GetDouble(), Is.EqualTo(c.CoveringIntensity), at);
            Assert.That(D(Column(root, "fraction")[i]), Is.EqualTo(Finite(c.Fraction)), at);
            Assert.That(Column(root, "unmodified_quantified")[i].GetBoolean(), Is.EqualTo(c.UnmodifiedQuantified), at);
            Assert.That(I(Column(root, "covering_peptidoforms")[i]), Is.EqualTo(c.CoveringPeptidoforms), at);
            Assert.That(I(Column(root, "modified_peptidoforms")[i]), Is.EqualTo(c.ModifiedPeptidoforms), at);
        }
    }

    private static void AssertPairsAre(JsonElement root, IReadOnlyList<PtmPair> expected)
    {
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));
        for (int i = 0; i < expected.Count; i++)
        {
            PtmPair p = expected[i];
            string at = $"pair {i} ({p.SiteA.Key} / {p.SiteB.Key})";
            Assert.That(Column(root, "result_type")[i].GetString(), Is.EqualTo(PtmQtl.PairTypeName(p.ResultType)), at);
            Assert.That(Column(root, "site_a")[i].GetString(), Is.EqualTo(p.SiteA.Key), at);
            Assert.That(Column(root, "site_b")[i].GetString(), Is.EqualTo(p.SiteB.Key), at);
            Assert.That(Column(root, "same_protein")[i].GetBoolean(), Is.EqualTo(p.SameProtein), at);
            Assert.That(Column(root, "overlapping")[i].GetBoolean(), Is.EqualTo(p.Overlapping), at);
            Assert.That(D(Column(root, "statistic")[i]), Is.EqualTo(Finite(p.Statistic)), at);
            Assert.That(D(Column(root, "p_value")[i]), Is.EqualTo(Finite(p.PValue)), at);
            Assert.That(D(Column(root, "q")[i]), Is.EqualTo(Finite(p.Q)), at);
            Assert.That(Column(root, "fdr_family")[i].GetString(), Is.EqualTo(p.FdrFamily), at);
            Assert.That(Column(root, "n")[i].GetInt32(), Is.EqualTo(p.N), at);
            Assert.That(Column(root, "statistic_n")[i].GetInt32(), Is.EqualTo(p.StatisticN), at);
            Assert.That(S(Column(root, "spearman_method")[i]),
                Is.EqualTo(p.SpearmanMethod is SpearmanPValueMethod m ? Statistics.SpearmanMethodName(m) : null), at);
            Assert.That(Column(root, "run_level_removed")[i].GetBoolean(), Is.EqualTo(p.RunLevelRemoved), at);
            Assert.That(S(Column(root, "test_status")[i]),
                Is.EqualTo(p.TestStatus is PairTestStatus t ? PtmQtl.TestStatusName(t) : null), at);
            Assert.That(Column(root, "bounds_included")[i].GetBoolean(), Is.EqualTo(p.BoundsIncluded), at);
            Assert.That(Column(root, "run_covariate_removed")[i].GetBoolean(), Is.EqualTo(p.RunCovariateRemoved), at);
        }
    }

    // ---- the committed fixtures --------------------------------------------------------------------

    private static bool Regenerate => Environment.GetEnvironmentVariable("PYMZLIB_REGENERATE_PTMQTL_FIXTURES") == "1";

    private static void AssertFixture(string name, string expected)
    {
        string path = Path.Combine(Fixtures, name);
        if (Regenerate)
        {
            Directory.CreateDirectory(Fixtures);
            File.WriteAllText(path, expected, new UTF8Encoding(false));
        }
        Assert.That(File.Exists(path), Is.True, $"{path} is missing; regenerate with PYMZLIB_REGENERATE_PTMQTL_FIXTURES=1");
        Assert.That(File.ReadAllText(path).Replace("\r\n", "\n"), Is.EqualTo(expected), name);
    }

    /// <summary>
    /// The input tables the Python docstrings, the ptmqtl guide and the recorded fixtures use are exactly
    /// mzLib's hand-computed cases above. PYMZLIB_REGENERATE_PTMQTL_FIXTURES=1 rewrites them.
    /// </summary>
    [Test]
    public void Fixtures_AreMzLibsHandComputedCases()
    {
        AssertFixture("cases.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Cases()));
        AssertFixture("fractions.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Fractions()));
        AssertFixture("fractions_samples.tsv", PtmQtlCases.SamplesTsv());
        AssertFixture("batch.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Batch()));
        AssertFixture("trait.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Trait()));
        AssertFixture("trait_runs.tsv", PtmQtlCases.TraitRunsTsv());
        AssertFixture("bounds.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Bounds()));
        AssertFixture("depth.tsv", PtmQtlCases.ObservationsTsv(PtmQtlCases.Depth()));
        AssertFixture("depth_psms.tsv", PtmQtlCases.DepthPsmsTsv());
    }

    /// <summary>
    /// trait_occupancy.tsv, the pivot trait-effects reads in the docs, is what <c>ptmqtl occupancy --out</c>
    /// writes from trait.tsv today, byte for byte.
    /// </summary>
    [Test]
    public void TraitOccupancyFixture_IsWhatOccupancyOutWrites()
    {
        string observations = Path.Combine(Fixtures, "trait.tsv");
        if (!File.Exists(observations))
            Assert.Ignore("trait.tsv is generated by Fixtures_AreMzLibsHandComputedCases");
        string written = Path.Combine(_temp, "trait_occupancy.tsv");
        Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", observations, "--out", written);
        AssertFixture("trait_occupancy.tsv", File.ReadAllText(written));
    }

    // ---- ptmqtl occupancy --------------------------------------------------------------------------

    [Test]
    public void Occupancy_ReproducesMzLibsHandComputedCases_ThroughTheWire()
    {
        var observations = PtmQtlCases.Cases();
        JsonElement root = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", WriteObservations(observations));
        AssertOccupancyIs(root, SiteOccupancyCalculator.Calculate(observations));

        string phos = PtmQtlCases.Phos;
        int ratio = RowOf(root, $"RATIO:S13:{phos}", "r1");
        Assert.That(D(Column(root, "fraction")[ratio]), Is.EqualTo(30.0 / 200).Within(1e-15));
        Assert.That(Column(root, "unmodified_quantified")[ratio].GetBoolean(), Is.True);
        int ceiling = RowOf(root, $"CEILING:S13:{phos}", "r1");
        Assert.That(D(Column(root, "fraction")[ceiling]), Is.EqualTo(1.0));
        Assert.That(Column(root, "unmodified_quantified")[ceiling].GetBoolean(), Is.False, "seen only modified: a ceiling");

        string State(string run) => Column(root, "state")[RowOf(root, $"STATES:S13:{phos}", run)].GetString()!;
        Assert.That(State("floor"), Is.EqualTo("floor"));
        Assert.That(State("countOnly"), Is.EqualTo("count_only"));
        Assert.That(State("notDetected"), Is.EqualTo("not_detected"));
        Assert.That(State("quant"), Is.EqualTo("quantified"));
        Assert.That(D(Column(root, "fraction")[RowOf(root, $"STATES:S13:{phos}", "floor")]), Is.Null, "a floor is not a zero");
        Assert.That(D(Column(root, "fraction")[RowOf(root, $"STATES:S13:{phos}", "quant")]), Is.EqualTo(0.1).Within(1e-15));

        double shared = D(Column(root, "fraction")[RowOf(root, $"SHARED_POSITION:S13:{phos}", "r1")])!.Value
                      + D(Column(root, "fraction")[RowOf(root, $"SHARED_POSITION:S13:{PtmQtlCases.HexNAc}", "r1")])!.Value;
        Assert.That(shared, Is.EqualTo(0.5).Within(1e-15), "two modifications at one position share the denominator");

        var proteins = Column(root, "protein_accession").Select(e => e.GetString()).Distinct().ToList();
        Assert.That(proteins, Does.Contain("NTERM_AFTER_MET").And.Contain("UNREGISTERED"));
        Assert.That(proteins, Has.None.EqualTo("EXCLUDED").And.None.EqualTo("PEPTIDE_NTERM").And.None.EqualTo("NTERM_AFTER_K"),
            "Common Variable/Fixed, a peptide N-terminal modification and a residue-2 start after K are not sites");
        Assert.That(Column(root, "position")[RowOf(root, "NTERM_AFTER_MET:S2:UniProt:N-acetylserine on S", "r1")].GetInt32(), Is.EqualTo(2),
            "a protein N-terminal modification sits at the residue it is on");

        Assert.That(root.GetProperty("observation_count").GetInt32(), Is.EqualTo(observations.Count));
        Assert.That(root.GetProperty("samples_file").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("sample_count").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("modifications").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("written").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("state_counts").GetProperty("floor").GetInt32(), Is.EqualTo(1));
        Assert.That(root.GetProperty("ignored_columns").GetArrayLength(), Is.EqualTo(0));
        Assert.That(Caveats(root), Has.Some.StartsWith("Intensity occupancy (mzLib's SiteOccupancyCalculator"));
        Assert.That(Caveats(root), Has.Some.StartsWith("2 cells are floor or not_detected"));
    }

    [Test]
    public void Occupancy_WithSamples_SumsEachSamplesFractionsBeforeTheRatio()
    {
        var observations = PtmQtlCases.Fractions();
        string path = WriteObservations(observations);
        string samples = Write("samples.tsv", PtmQtlCases.SamplesTsv());
        JsonElement root = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--samples", samples);
        var map = PtmQtlCases.FractionSamples.ToDictionary(s => s.Run, s => s.Sample);
        AssertOccupancyIs(root, SiteOccupancyCalculator.Calculate(observations, map));

        string phos = PtmQtlCases.Phos;
        int a = RowOf(root, $"SUMMED:S13:{phos}", "A");
        Assert.That(D(Column(root, "fraction")[a]), Is.EqualTo(40.0 / 200).Within(1e-15), "(30 + 10) / (100 + 100), not the mean of 0.3 and 0.1");
        Assert.That(Column(root, "state")[RowOf(root, $"SUMMED:S13:{phos}", "B")].GetString(), Is.EqualTo("not_detected"));
        Assert.That(Column(root, "state")[RowOf(root, $"FLOOR_ACROSS:S13:{phos}", "A")].GetString(), Is.EqualTo("floor"),
            "identified unquantified in one fraction, covered in another: a floor for the sample");
        Assert.That(root.GetProperty("run_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "A", "B" }));
        Assert.That(root.GetProperty("run_count").GetInt32(), Is.EqualTo(3));
        Assert.That(root.GetProperty("sample_count").GetInt32(), Is.EqualTo(2));
        Assert.That(root.GetProperty("samples_file").GetString(), Is.EqualTo(Path.GetFullPath(samples)));
        Assert.That(Caveats(root), Has.Some.StartsWith("With --samples"));

        // Per run, the same table gives each fraction its own ratio.
        JsonElement perRun = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path);
        Assert.That(D(Column(perRun, "fraction")[RowOf(perRun, $"SUMMED:S13:{phos}", "f1")]), Is.EqualTo(0.3).Within(1e-15));
        Assert.That(Column(perRun, "state")[RowOf(perRun, $"FLOOR_ACROSS:S13:{phos}", "f1")].GetString(), Is.EqualTo("count_only"));
    }

    [Test]
    public void Occupancy_Filters_AreMzLibsIncludeModification()
    {
        var observations = PtmQtlCases.Cases();
        string path = WriteObservations(observations);
        string phos = PtmQtlCases.Phos;

        JsonElement byName = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path,
            "--modifications", $"{phos}\tUniProt:N-acetylserine on S");
        AssertOccupancyIs(byName, SiteOccupancyCalculator.Calculate(observations,
            includeModification: m => m == phos || m == "UniProt:N-acetylserine on S"));
        Assert.That(byName.GetProperty("modifications").EnumerateArray().Select(e => e.GetString()),
            Is.EqualTo(new[] { phos, "UniProt:N-acetylserine on S" }));

        JsonElement byCategory = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--categories", "Common Biological");
        AssertOccupancyIs(byCategory, SiteOccupancyCalculator.Calculate(observations,
            includeModification: m => m.StartsWith("Common Biological:", StringComparison.Ordinal)));
        Assert.That(Column(byCategory, "category").Select(e => e.GetString()).Distinct(), Is.EqualTo(new[] { "Common Biological" }));

        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--modifications", phos, "--categories", "x"),
            Does.Contain("not both"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--categories", $"Common Biological\t"),
            Does.Contain("blank name"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--categories"),
            Does.Contain("without a value"));
    }

    [Test]
    public void Occupancy_Out_WritesThePlainTableAndReturnsNoColumns()
    {
        var observations = PtmQtlCases.Cases();
        string path = WriteObservations(observations);
        string output = Path.Combine(_temp, "nested", "occupancy.tsv");
        JsonElement root = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--out", output);

        var expected = SiteOccupancyCalculator.Calculate(observations);
        Assert.That(root.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));
        JsonElement written = root.GetProperty("written");
        Assert.That(written.GetProperty("path").GetString(), Is.EqualTo(Path.GetFullPath(output)));
        Assert.That(written.GetProperty("format").GetString(), Is.EqualTo("tsv"));
        Assert.That(written.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));

        string[] lines = File.ReadAllText(output).Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.That(lines[0], Is.EqualTo(string.Join('\t', PtmQtl.OccupancyColumns)));
        Assert.That(lines, Has.Length.EqualTo(expected.Count + 1));
        string[] floor = lines.Single(l => l.Contains("\tfloor\t")).Split('\t');
        Assert.That(floor[Array.IndexOf(PtmQtl.OccupancyColumns, "fraction")], Is.Empty, "null is an empty cell");
        Assert.That(lines.Count(l => l.Contains("\tfalse\t")), Is.GreaterThan(0));

        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--out", Path.Combine(_temp, "o.csv")),
            Does.Contain(".tsv"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--out", path),
            Does.Contain("must differ from --observations"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", path, "--out"), Does.Contain("without a value"));
    }

    [Test]
    public void Occupancy_RefusesMalformedObservations_AsUsageErrors()
    {
        const string header = "run\tfull_sequence\tprotein_accession\tstart_residue\tend_residue\tintensity\n";
        string Refused(string body, params string[] extra) =>
            Usage(PtmQtl.Occupancy, new[] { "ptmqtl", "occupancy", "--observations", Write($"bad{Guid.NewGuid():N}.tsv", body) }.Concat(extra).ToArray());

        Assert.That(Refused("run\tfull_sequence\tprotein_accession\tstart_residue\tintensity\nr1\tPEPSK\tP1\t10\t5\n"),
            Does.Contain("no column end_residue"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\tten\t14\t5\n"), Does.Contain("start_residue must be an integer"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\tlots\n"), Does.Contain("is not a number"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\n"), Does.Contain("has 5 fields; the header has 6"));
        Assert.That(Refused(header), Does.Contain("has a header but no rows"));
        Assert.That(Refused(header.TrimEnd('\n') + "\tprevious_residue\nr1\tPEPSK\tP1\t10\t14\t5\tMK\n"), Does.Contain("one residue letter"));
        // mzLib's own refusals cross as usage errors, in mzLib's words.
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t20\t5\n"), Does.Contain("has 5 residues but spans 10-20"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\t5\nr1\tPEPSK\tP1\t10\t14\t5\n"), Does.Contain("not one per PSM"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\t-1\n"), Does.Contain("finite and non-negative"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\t5\n", "--samples", Write("s.tsv", "run\tsample\nr2\tA\n")),
            Does.Contain("Run r1 has no sample"));
        Assert.That(Refused(header + "r1\tPEPSK\tP1\t10\t14\t5\n", "--samples", Write("s2.tsv", "run\tsample\nr1\n")),
            Does.Contain("it needs two"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", Path.Combine(_temp, "absent.tsv")),
            Does.Contain("was not found"));
        Assert.That(Usage(PtmQtl.Occupancy, "ptmqtl", "occupancy"), Does.Contain("Missing required option --observations"));
    }

    [Test]
    public void Occupancy_ReportsColumnsItDidNotRead()
    {
        const string body = "run\tfull_sequence\tprotein_accession\tstart_residue\tend_residue\tintensity\tprevious_residu\tscore\n" +
                            "r1\tPEPS[Common Biological:Phosphorylation on S]K\tP1\t10\t14\t5\tK\t30\n";
        JsonElement root = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", Write("typo.tsv", body));
        Assert.That(root.GetProperty("ignored_columns").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "previous_residu", "score" }));
        Assert.That(Caveats(root), Has.Some.Contains("previous_residu, score"));
    }

    [Test]
    public void WireNames_CoverEveryValueMzLibHas()
    {
        Assert.That(Enum.GetValues<OccupancyState>().Select(PtmQtl.StateName), Is.EqualTo(new[] { "quantified", "floor", "count_only", "not_detected" }));
        Assert.That(Enum.GetValues<PairResultType>().Select(PtmQtl.PairTypeName), Is.EquivalentTo(new[] { "physical", "covarying" }));
        Assert.That(PtmQtl.TraitStatusName("Fitted"), Is.EqualTo("fitted"));
        Assert.That(PtmQtl.TraitStatusName("BelowSupport"), Is.EqualTo("below_support"));
        foreach (FeatureFitStatus status in new[] { FeatureFitStatus.TooFewObservations, FeatureFitStatus.RankDeficient,
                     FeatureFitStatus.TooFewGroups, FeatureFitStatus.VarianceRatioAtLimit, FeatureFitStatus.NotConverged })
            Assert.That(PtmQtl.TraitStatuses, Does.Contain(PtmQtl.TraitStatusName(status.ToString())));
        Assert.Throws<InvalidOperationException>(() => PtmQtl.TraitStatusName("Separated"), "a logistic-only outcome MixedModel never reports");
        Assert.Throws<InvalidOperationException>(() => PtmQtl.TraitStatusName("SomethingNew"));
        Assert.That(Enum.GetValues<PairTestStatus>().Select(PtmQtl.TestStatusName),
            Is.EqualTo(new[] { "tested", "not_estimable", "too_few_shared_runs", "both_sites_mostly_bounds" }));
        Assert.Throws<InvalidOperationException>(() => PtmQtl.TestStatusName((PairTestStatus)99));
        Assert.That(Enum.GetValues<OccupancyScale>().Select(PtmQtl.ScaleName), Is.EqualTo(new[] { "natural_log_odds", "log2_odds" }));
        Assert.Throws<InvalidOperationException>(() => PtmQtl.ScaleName((OccupancyScale)99));
        Assert.That(Program.Verbs, Is.SupersetOf(new[] { "ptmqtl occupancy", "ptmqtl pairs", "ptmqtl trait-effects" }));
    }

    // ---- ptmqtl pairs ------------------------------------------------------------------------------

    [Test]
    public void Pairs_Physical_ReproducesMzLibsCoOccupancy()
    {
        var observations = PtmQtlCases.Cases();
        JsonElement root = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", WriteObservations(observations), "--type", "physical");
        var expected = PtmPairEngine.Physical(observations);
        AssertPairsAre(root, expected);

        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(1));
        Assert.That(D(Column(root, "statistic")[0]), Is.EqualTo(0.25).Within(1e-15), "25 of 100 covering both positions");
        Assert.That(Column(root, "n")[0].GetInt32(), Is.EqualTo(2), "identified together in two runs");
        Assert.That(Column(root, "statistic_n")[0].GetInt32(), Is.EqualTo(1), "quantified in one");
        Assert.That(D(Column(root, "p_value")[0]), Is.Null);
        Assert.That(Column(root, "fdr_family")[0].GetString(), Is.EqualTo("P:intra"));
        foreach (string covaryingOnly in new[] { "covarying_count", "covarying_run_count", "min_quantified_fraction", "exclude_ceiling",
                     "remove_run_level", "run_level_group", "min_runs_for_run_level" })
            Assert.That(root.GetProperty(covaryingOnly).ValueKind, Is.EqualTo(JsonValueKind.Null), covaryingOnly);
        Assert.That(root.GetProperty("physical_count").GetInt32(), Is.EqualTo(1));
        Assert.That(Caveats(root), Has.Some.StartsWith("physical (mzLib's PtmPairEngine.Physical"));
    }

    [Test]
    public void Pairs_Covarying_BatchShift_IsMzLibsSpearmanWithAndWithoutTheRunLevel()
    {
        var observations = PtmQtlCases.Batch();
        string path = WriteObservations(observations);
        var occupancy = SiteOccupancyCalculator.Calculate(observations);
        string s = $"PAIR_S:S13:{PtmQtlCases.Phos}", t = $"PAIR_T:S43:{PtmQtlCases.Phos}";
        int PairRow(JsonElement root) => Enumerable.Range(0, root.GetProperty("row_count").GetInt32())
            .Single(i => Column(root, "site_a")[i].GetString() == s && Column(root, "site_b")[i].GetString() == t);

        JsonElement plain = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying");
        AssertPairsAre(plain, PtmPairEngine.CoVarying(occupancy, observations));
        Assert.That(D(Column(plain, "statistic")[PairRow(plain)]), Is.GreaterThan(0.7), "the batch alone correlates the two sites");
        Assert.That(plain.GetProperty("min_quantified_fraction").GetDouble(), Is.EqualTo(PtmPairEngine.DefaultMinQuantifiedFraction));
        Assert.That(plain.GetProperty("exclude_ceiling").ValueKind, Is.EqualTo(JsonValueKind.Null),
            "absent, mzLib resolves its nullable excludeCeiling itself");
        Assert.That(plain.GetProperty("remove_run_level").GetBoolean(), Is.False);
        Assert.That(plain.GetProperty("run_level_group").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(plain.GetProperty("covarying_run_count").GetInt32(), Is.EqualTo(12));
        Assert.That(plain.GetProperty("physical_count").ValueKind, Is.EqualTo(JsonValueKind.Null));

        JsonElement adjusted = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--remove-run-level");
        AssertPairsAre(adjusted, PtmPairEngine.CoVarying(occupancy, observations, removeRunLevel: true));
        int row = PairRow(adjusted);
        Assert.That(Math.Abs(D(Column(adjusted, "statistic")[row])!.Value), Is.LessThan(0.1), "given the batch, the pair is unrelated");
        Assert.That(Column(adjusted, "n")[row].GetInt32(), Is.EqualTo(12));
        Assert.That(Column(adjusted, "run_level_removed")[row].GetBoolean(), Is.True);
        Assert.That(adjusted.GetProperty("run_level_group").GetString(), Is.EqualTo("none"));
        Assert.That(adjusted.GetProperty("min_runs_for_run_level").GetInt32(), Is.EqualTo(PtmPairEngine.MinRunsForRunLevel));
        Assert.That(adjusted.GetProperty("min_sites_for_run_level").GetInt32(), Is.EqualTo(PtmPairEngine.MinSitesForRunLevel));
        Assert.That(Caveats(adjusted), Has.Some.StartsWith("With --remove-run-level"));

        // Both types in one call: physical rows first (none here, no peptidoform carries two sites), then co-varying.
        JsonElement both = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "both");
        Assert.That(both.GetProperty("physical_count").GetInt32(), Is.EqualTo(0));
        Assert.That(both.GetProperty("covarying_count").GetInt32(), Is.EqualTo(plain.GetProperty("row_count").GetInt32()));
    }

    [Test]
    public void Pairs_RunLevelGroups_AreMzLibsRunLevelGroupFunc()
    {
        var observations = PtmQtlCases.Grouped();
        string path = WriteObservations(observations);
        var occupancy = SiteOccupancyCalculator.Calculate(observations);
        string Chemistry(ModificationSite site) => site.Modification.Contains("GG") ? "GG" : "phospho";
        var perChemistry = PtmPairEngine.CoVarying(occupancy, observations, removeRunLevel: true, runLevelGroup: Chemistry);

        string groups = Write("groups.tsv", $"modification\tgroup\n{PtmQtlCases.Gg}\tGG\n{PtmQtlCases.Phos}\tphospho\n");
        JsonElement byFile = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--remove-run-level", "--run-level-groups", groups);
        AssertPairsAre(byFile, perChemistry);
        Assert.That(byFile.GetProperty("run_level_group").GetString(), Is.EqualTo("file"));
        Assert.That(byFile.GetProperty("run_level_groups_file").GetString(), Is.EqualTo(Path.GetFullPath(groups)));

        // The two chemistries here are two categories and two modification names, so all three rules agree.
        foreach (string rule in new[] { "category", "modification" })
        {
            JsonElement byRule = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
                "--remove-run-level", "--run-level-group", rule);
            AssertPairsAre(byRule, perChemistry);
            Assert.That(byRule.GetProperty("run_level_group").GetString(), Is.EqualTo(rule));
        }

        string s = $"GG_S:K13:{PtmQtlCases.Gg}", t = $"GG_T:K13:{PtmQtlCases.Gg}";
        int row = Enumerable.Range(0, byFile.GetProperty("row_count").GetInt32())
            .Single(i => Column(byFile, "site_a")[i].GetString() == s && Column(byFile, "site_b")[i].GetString() == t);
        Assert.That(Math.Abs(D(Column(byFile, "statistic")[row])!.Value), Is.LessThan(0.1), "given the GG level, the GG sites are unrelated");

        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--remove-run-level",
                "--run-level-groups", Write("partial.tsv", $"modification\tgroup\n{PtmQtlCases.Gg}\tGG\n")),
            Does.Contain("no group for 1 modification(s)").And.Contain(PtmQtlCases.Phos));
    }

    [Test]
    public void Pairs_KeepCeilingAndMinQuantifiedFraction_AreMzLibsArguments()
    {
        var observations = PtmQtlCases.Ceiling();
        string path = WriteObservations(observations);
        var occupancy = SiteOccupancyCalculator.Calculate(observations);

        JsonElement excluded = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--min-quantified-fraction", "0.5");
        AssertPairsAre(excluded, PtmPairEngine.CoVarying(occupancy, observations, 0.5));
        Assert.That(Column(excluded, "n")[0].GetInt32(), Is.EqualTo(5));
        Assert.That(D(Column(excluded, "statistic")[0]), Is.EqualTo(-1).Within(1e-12));
        Assert.That(excluded.GetProperty("min_quantified_fraction").GetDouble(), Is.EqualTo(0.5));

        JsonElement kept = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--min-quantified-fraction", "0.5", "--keep-ceiling");
        AssertPairsAre(kept, PtmPairEngine.CoVarying(occupancy, observations, 0.5, excludeCeiling: false));
        Assert.That(Column(kept, "n")[0].GetInt32(), Is.EqualTo(8));
        Assert.That(D(Column(kept, "statistic")[0]), Is.EqualTo(0.5).Within(1e-12), "the shared ceilings correlate the sites through detection");
        Assert.That(kept.GetProperty("exclude_ceiling").GetBoolean(), Is.False);

        // At mzLib's default of 0.7 neither site is quantified (non-ceiling) in enough of the 8 runs.
        JsonElement guarded = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying");
        Assert.That(guarded.GetProperty("row_count").GetInt32(), Is.EqualTo(0));
    }

    [Test]
    public void Pairs_IncludeBounds_EntersOneFormCellsAsTiedRanks()
    {
        // mzLib's BoundsEnterTypeAAsTiedRanks, through the observations.
        var observations = PtmQtlCases.Bounds();
        string path = WriteObservations(observations);
        var occupancy = SiteOccupancyCalculator.Calculate(observations);

        // Without bounds, each site has a value in 6 of 10 runs: below 70 %, so no pair.
        JsonElement plain = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying");
        Assert.That(plain.GetProperty("row_count").GetInt32(), Is.EqualTo(0));
        Assert.That(plain.GetProperty("include_bounds").GetBoolean(), Is.False);
        Assert.That(plain.GetProperty("min_untied_share").ValueKind, Is.EqualTo(JsonValueKind.Null));

        JsonElement bounded = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--include-bounds");
        var expected = PtmPairEngine.CoVarying(occupancy, observations, includeBounds: true);
        AssertPairsAre(bounded, expected);
        Assert.That(bounded.GetProperty("row_count").GetInt32(), Is.EqualTo(1));
        Assert.That(Column(bounded, "n")[0].GetInt32(), Is.EqualTo(9), "r8 is left out: t's count-only cell there has no rank");
        Assert.That(Column(bounded, "test_status")[0].GetString(), Is.EqualTo("tested"));
        Assert.That(Column(bounded, "bounds_included")[0].GetBoolean(), Is.True);
        // mzLib's own expectation: low bounds tie below every fraction and high bounds above.
        double lo = -5, hi = 9;
        var vs = new[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, lo, lo, hi, hi };
        var vt = new[] { 0.15, 0.1, 0.35, 0.3, 0.6, 0.45, hi, lo, double.NaN, lo };
        Assert.That(D(Column(bounded, "statistic")[0]), Is.EqualTo(SpearmanCorrelation.Correlate(vs, vt).Rho).Within(1e-12));
        Assert.That(bounded.GetProperty("include_bounds").GetBoolean(), Is.True);
        Assert.That(bounded.GetProperty("min_untied_share").GetDouble(), Is.EqualTo(PtmPairEngine.MinUntiedShare));
        Assert.That(bounded.GetProperty("exclude_ceiling").ValueKind, Is.EqualTo(JsonValueKind.Null), "mzLib decided: a ceiling is a bound");
        Assert.That(Caveats(bounded), Has.Some.StartsWith("With --include-bounds, a pair is tested only if"));
        Assert.That(Caveats(bounded), Has.Some.Contains("enters as a bound"));

        // --keep-ceiling sends excludeCeiling false, which mzLib accepts with bounds: the same pair.
        JsonElement kept = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--include-bounds", "--keep-ceiling");
        AssertPairsAre(kept, PtmPairEngine.CoVarying(occupancy, observations, excludeCeiling: false, includeBounds: true));
        Assert.That(D(Column(kept, "statistic")[0]), Is.EqualTo(D(Column(bounded, "statistic")[0])));
        Assert.That(kept.GetProperty("exclude_ceiling").GetBoolean(), Is.False);
    }

    [Test]
    public void Pairs_IncludeBounds_LeavesAPairOfMostlyBoundSitesUntested()
    {
        // mzLib's APairWhoseSitesAreBothMostlyBoundsIsNotTested, through the observations.
        var observations = PtmQtlCases.MostlyBounds();
        var occupancy = SiteOccupancyCalculator.Calculate(observations);
        JsonElement root = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", WriteObservations(observations),
            "--type", "covarying", "--include-bounds");
        AssertPairsAre(root, PtmPairEngine.CoVarying(occupancy, observations, includeBounds: true));

        string phos = PtmQtlCases.Phos;
        int Row(string a, string b) => Enumerable.Range(0, root.GetProperty("row_count").GetInt32())
            .Single(i => Column(root, "site_a")[i].GetString() == a && Column(root, "site_b")[i].GetString() == b);
        int st = Row($"MOSTLY_S:S13:{phos}", $"MOSTLY_T:S43:{phos}");
        Assert.That(Column(root, "test_status")[st].GetString(), Is.EqualTo("both_sites_mostly_bounds"));
        Assert.That(Column(root, "n")[st].GetInt32(), Is.EqualTo(20));
        Assert.That(D(Column(root, "statistic")[st]), Is.Null);
        Assert.That(D(Column(root, "q")[st]), Is.Null, "left out of the adjustment");
        Assert.That(Column(root, "spearman_method")[st].GetString(), Is.EqualTo("not_estimable"));
        int su = Row($"MOSTLY_S:S13:{phos}", $"MOSTLY_U:S7:{phos}");
        Assert.That(Column(root, "test_status")[su].GetString(), Is.EqualTo("tested"), "one site with values in most runs is enough");
        Assert.That(D(Column(root, "q")[su]), Is.Not.Null);
        Assert.That(Caveats(root), Has.Some.EqualTo(
            "1 co-varying pairs were not tested (1 both_sites_mostly_bounds): their statistic, p_value and q are null, and spearman_method is not_estimable."));
    }

    [Test]
    public void Pairs_RunCovariate_IsMzLibsPartialCorrelationGivenTheRunsValue()
    {
        // mzLib's ARunCovariateIsPartialledOut, through the observations: two sites that rise with each run's depth.
        var observations = PtmQtlCases.Depth();
        string path = WriteObservations(observations);
        var occupancy = SiteOccupancyCalculator.Calculate(observations);
        string psms = Write("psms.tsv", PtmQtlCases.DepthPsmsTsv());
        string s = $"DEPTH_S:S13:{PtmQtlCases.Phos}", t = $"DEPTH_T:S43:{PtmQtlCases.Phos}";
        int PairRow(JsonElement root) => Enumerable.Range(0, root.GetProperty("row_count").GetInt32())
            .Single(i => Column(root, "site_a")[i].GetString() == s && Column(root, "site_b")[i].GetString() == t);

        JsonElement given = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--run-covariate", psms);
        AssertPairsAre(given, PtmPairEngine.CoVarying(occupancy, observations, runCovariate: PtmQtlCases.DepthMap()));
        int row = PairRow(given);
        double[] noiseS = { 1, 4, 2, 6, 3, 5 }, noiseT = { 6, 2, 1, 4, 3, 5 };
        var vs = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseS[r % 6] + PtmQtlCases.DepthPsms[r] * 1e-6).ToArray();
        var vt = Enumerable.Range(0, 12).Select(r => (r < 6 ? 0.1 : 0.3) + 0.01 * noiseT[r % 6] + PtmQtlCases.DepthPsms[r] * 2e-6).ToArray();
        var expected = SpearmanCorrelation.PartialCorrelate(vs, vt, PtmQtlCases.DepthPsms);
        Assert.That(D(Column(given, "statistic")[row]), Is.EqualTo(expected.Rho).Within(1e-12), "mzLib's hand expectation");
        Assert.That(D(Column(given, "p_value")[row]), Is.EqualTo(expected.PValue).Within(1e-12));
        Assert.That(Column(given, "run_covariate_removed")[row].GetBoolean(), Is.True);
        Assert.That(Column(given, "run_level_removed")[row].GetBoolean(), Is.False);
        Assert.That(given.GetProperty("run_covariate_file").GetString(), Is.EqualTo(Path.GetFullPath(psms)));
        Assert.That(given.GetProperty("run_covariate_name").GetString(), Is.EqualTo("target_psms"));
        Assert.That(given.GetProperty("run_covariate_runs_missing").GetArrayLength(), Is.EqualTo(0));
        Assert.That(given.GetProperty("min_runs_for_run_level").GetInt32(), Is.EqualTo(PtmPairEngine.MinRunsForRunLevel),
            "a pair partial on a covariate is held to the same shared runs");
        Assert.That(given.GetProperty("min_sites_for_run_level").ValueKind, Is.EqualTo(JsonValueKind.Null));
        Assert.That(Caveats(given), Has.Some.StartsWith("With --run-covariate, rho and p_value are the partial Spearman correlation given each run's target_psms"));

        // With the run level too: both.
        JsonElement both = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", psms, "--remove-run-level");
        AssertPairsAre(both, PtmPairEngine.CoVarying(occupancy, observations, removeRunLevel: true, runCovariate: PtmQtlCases.DepthMap()));

        // A run with no value is left out of every pair and reported; below MinRunsForRunLevel the pair is not tested.
        string lines = PtmQtlCases.DepthPsmsTsv();
        string missing = Write("missing.tsv", string.Join('\n', lines.Split('\n').Where(l => !l.StartsWith("r00"))));
        JsonElement withoutOne = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--run-covariate", missing);
        Assert.That(Column(withoutOne, "n")[PairRow(withoutOne)].GetInt32(), Is.EqualTo(11));
        Assert.That(withoutOne.GetProperty("run_covariate_runs_missing").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "r00" }));
        Assert.That(Caveats(withoutOne), Has.Some.StartsWith("1 run(s) have no target_psms in the --run-covariate file"));
        string eight = Write("eight.tsv", string.Join('\n', lines.Split('\n').Take(9)) + "\n");
        JsonElement few = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying", "--run-covariate", eight);
        Assert.That(Column(few, "test_status")[PairRow(few)].GetString(), Is.EqualTo("too_few_shared_runs"));
        Assert.That(D(Column(few, "p_value")[PairRow(few)]), Is.Null);

        // mzLib's refusal of a value that is not finite, in mzLib's words.
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("nan.tsv", "run\ttarget_psms\nr00\tNA\n")), Does.Contain("Every run covariate value must be finite"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("one.tsv", "run\nr00\n")), Does.Contain("needs a header of two fields"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("wide.tsv", "run\tdepth\nr00\t1\t2\n")), Does.Contain("has 3 fields"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("twice.tsv", "run\tdepth\nr00\t1\nr00\t2\n")), Does.Contain("more than once"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("blankrun.tsv", "run\tdepth\n\t1\n")), Does.Contain("empty run"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("empty.tsv", "run\tdepth\n")), Does.Contain("has a header but no runs"));
        Assert.That(Usage(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", path, "--type", "covarying",
            "--run-covariate", Write("word.tsv", "run\tdepth\nr00\tdeep\n")), Does.Contain("is not a number"));
    }

    [Test]
    public void Occupancy_CountsPeptidoformsPerCell_AsMzLibDoes()
    {
        // mzLib's ACoveredZeroNeedsOnlyAnIdentificationAndCellsCountPeptidoforms.
        string phos = PtmQtlCases.Phos;
        var observations = new List<PeptidoformObservation>
        {
            PtmQtlCases.Obs("idOnly", "PEPSK", 10, double.NaN, "P1"),
            PtmQtlCases.Obs("quant", $"PEPS[{phos}]K", 10, 5, "P1"), PtmQtlCases.Obs("quant", "PEPSK", 10, 45, "P1"),
            PtmQtlCases.Obs("quant", "AAPEPSKR", 8, 50, "P1"),
            PtmQtlCases.Obs("onlyModified", $"PEPS[{phos}]K", 10, double.NaN, "P1"),
            PtmQtlCases.Obs("bothForms", $"PEPS[{phos}]K", 10, double.NaN, "P1"), PtmQtlCases.Obs("bothForms", "PEPSK", 10, double.NaN, "P1"),
            PtmQtlCases.Obs("elsewhere", "GGGGK", 30, 20, "P1"),
        };
        JsonElement root = Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", WriteObservations(observations));
        AssertOccupancyIs(root, SiteOccupancyCalculator.Calculate(observations));
        (string State, int? Covering, int? Modified) Cell(string run)
        {
            int i = RowOf(root, $"P1:S13:{phos}", run);
            return (Column(root, "state")[i].GetString()!, I(Column(root, "covering_peptidoforms")[i]), I(Column(root, "modified_peptidoforms")[i]));
        }
        Assert.That(Cell("idOnly"), Is.EqualTo(("not_detected", (int?)1, (int?)0)), "a covered-zero needs only an identification");
        Assert.That(Cell("quant"), Is.EqualTo(("quantified", (int?)3, (int?)1)));
        Assert.That(Cell("onlyModified"), Is.EqualTo(("count_only", (int?)1, (int?)1)));
        Assert.That(Cell("bothForms"), Is.EqualTo(("count_only", (int?)2, (int?)1)));
        Assert.That(Caveats(root), Has.Some.StartsWith("covering_peptidoforms and modified_peptidoforms count identifications"));
    }

    [Test]
    public void Pairs_WithSamples_FeedsCombinedObservationsToPhysicalAndCoVarying()
    {
        string gg = PtmQtlCases.Gg, phos = PtmQtlCases.Phos;
        var observations = new List<PeptidoformObservation>
        {
            // CombinedObservationsFeedSameMoleculePairsPerSample.
            PtmQtlCases.Obs("f1", $"PEPS[{phos}]K[{gg}]R", 10, 20, "P1"), PtmQtlCases.Obs("f1", "PEPSKR", 10, 80, "P1"),
            PtmQtlCases.Obs("f2", $"PEPS[{phos}]K[{gg}]R", 10, 10, "P1"), PtmQtlCases.Obs("f2", "PEPSKR", 10, 90, "P1"),
            PtmQtlCases.Obs("g1", $"PEPS[{phos}]K[{gg}]R", 10, double.NaN, "P1"),
        };
        string samples = Write("samples.tsv", "run\tsample\nf1\tA\nf2\tA\ng1\tB\n");
        JsonElement root = Run(PtmQtl.Pairs, "ptmqtl", "pairs", "--observations", WriteObservations(observations),
            "--type", "both", "--samples", samples);
        var combined = SiteOccupancyCalculator.CombineObservations(observations, new Dictionary<string, string> { ["f1"] = "A", ["f2"] = "A", ["g1"] = "B" });
        var expected = PtmPairEngine.Physical(combined).Concat(PtmPairEngine.CoVarying(SiteOccupancyCalculator.Calculate(combined), combined)).ToList();
        AssertPairsAre(root, expected);
        Assert.That(Column(root, "n")[0].GetInt32(), Is.EqualTo(2), "two samples, not three runs");
        Assert.That(D(Column(root, "statistic")[0]), Is.EqualTo(30.0 / 200).Within(1e-15));
        Assert.That(root.GetProperty("run_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "A", "B" }));
    }

    [Test]
    public void Pairs_RefusesOptionsThatDoNotApply()
    {
        string path = WriteObservations(PtmQtlCases.Batch());
        string[] Call(params string[] extra) => new[] { "ptmqtl", "pairs", "--observations", path }.Concat(extra).ToArray();

        Assert.That(Usage(PtmQtl.Pairs, Call()), Does.Contain("Missing required option --type"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "all")), Does.Contain("must be one of physical, covarying, both"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "physical", "--remove-run-level")), Does.Contain("shapes co-varying pairs only"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "physical", "--min-quantified-fraction", "0.5")), Does.Contain("--min-quantified-fraction shapes"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--min-quantified-fraction", "0")), Does.Contain("above 0 and at most 1"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--min-quantified-fraction", "1.5")), Does.Contain("above 0 and at most 1"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--min-quantified-fraction")), Does.Contain("without a value"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--run-level-group", "category")), Does.Contain("give --remove-run-level as well"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--remove-run-level", "--run-level-group", "chemistry")),
            Does.Contain("must be one of category, modification"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--remove-run-level", "--run-level-group", "category",
            "--run-level-groups", path)), Does.Contain("not both"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--keep-ceiling", "false")), Does.Contain("takes no value"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "physical", "--include-bounds")), Does.Contain("--include-bounds shapes co-varying pairs only"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "physical", "--run-covariate", path)), Does.Contain("--run-covariate shapes co-varying pairs only"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--include-bounds", "yes")), Does.Contain("takes no value"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--run-covariate")), Does.Contain("without a value"));
        Assert.That(Usage(PtmQtl.Pairs, Call("--type", "covarying", "--run-covariate", Path.Combine(_temp, "absent.tsv"))), Does.Contain("was not found"));
    }

    [Test]
    public void CoVaryingDefaults_AreReadFromMzLibsDeclaration()
    {
        Assert.That(PtmQtl.CoVaryingDefault<double>("minQuantifiedFraction"), Is.EqualTo(PtmPairEngine.DefaultMinQuantifiedFraction));
        Assert.That(PtmQtl.CoVaryingDefault<double>("minQuantifiedFraction"), Is.EqualTo(0.7), "the value the spec states at the pin");
        Assert.That(typeof(PtmPairEngine).GetMethod(nameof(PtmPairEngine.CoVarying))!.GetParameters()
            .Single(p => p.Name == "excludeCeiling").DefaultValue, Is.Null, "mzLib resolves it inside CoVarying, so the bridge sends null");
        Assert.That(PtmQtl.CoVaryingDefault<bool>("removeRunLevel"), Is.False);
        Assert.That(PtmQtl.CoVaryingDefault<bool>("includeBounds"), Is.False);
        Assert.That(typeof(PtmPairEngine).GetMethod(nameof(PtmPairEngine.CoVarying))!.GetParameters()
            .Single(p => p.Name == "runCovariate").DefaultValue, Is.Null, "absent, the bridge sends none");
    }

    // ---- ptmqtl trait-effects ----------------------------------------------------------------------

    private string TraitOccupancy()
    {
        string output = Path.Combine(_temp, "trait_occupancy.tsv");
        Run(PtmQtl.Occupancy, "ptmqtl", "occupancy", "--observations", WriteObservations(PtmQtlCases.Trait(), "trait.tsv"), "--out", output);
        return output;
    }

    private static (Dictionary<string, double> Trait, Dictionary<string, string> Replicate, Dictionary<string, double[]> Covariates) TraitMaps() =>
        (PtmQtlCases.TraitRuns().ToDictionary(t => t.Run, t => t.Dose),
         PtmQtlCases.TraitRuns().ToDictionary(t => t.Run, t => t.Replicate),
         PtmQtlCases.TraitRuns().ToDictionary(t => t.Run, t => new[] { t.Injection }));

    private static void AssertEffectsAre(JsonElement root, IReadOnlyList<SiteTraitEffect> expected)
    {
        Assert.That(root.GetProperty("row_count").GetInt32(), Is.EqualTo(expected.Count));
        for (int i = 0; i < expected.Count; i++)
        {
            SiteTraitEffect e = expected[i];
            string at = e.Site.Key;
            Assert.That(Column(root, "site_key")[i].GetString(), Is.EqualTo(e.Site.Key), at);
            Assert.That(Column(root, "protein_accession")[i].GetString(), Is.EqualTo(e.Site.ProteinAccession), at);
            Assert.That(Column(root, "position")[i].GetInt32(), Is.EqualTo(e.Site.Position), at);
            Assert.That(Column(root, "residue")[i].GetString(), Is.EqualTo(e.Site.Residue.ToString()), at);
            Assert.That(Column(root, "modification")[i].GetString(), Is.EqualTo(e.Site.Modification), at);
            Assert.That(Column(root, "status")[i].GetString(), Is.EqualTo(PtmQtl.TraitStatusName(e.Status)), at);
            Assert.That(D(Column(root, "effect")[i]), Is.EqualTo(Finite(e.Effect)), at);
            Assert.That(D(Column(root, "standard_error")[i]), Is.EqualTo(Finite(e.StandardError)), at);
            Assert.That(D(Column(root, "df")[i]), Is.EqualTo(Finite(e.DegreesOfFreedom)), at);
            Assert.That(D(Column(root, "p_value")[i]), Is.EqualTo(Finite(e.PValue)), at);
            Assert.That(D(Column(root, "q")[i]), Is.EqualTo(Finite(e.Q)), at);
            Assert.That(Column(root, "runs")[i].GetInt32(), Is.EqualTo(e.Runs), at);
            Assert.That(Column(root, "replicates")[i].GetInt32(), Is.EqualTo(e.Replicates), at);
            Assert.That(D(Column(root, "median_fraction")[i]), Is.EqualTo(Finite(e.MedianFraction)), at);
            Assert.That(Column(root, "covariate_effects")[i].EnumerateArray().Select(D), Is.EqualTo(e.CovariateEffects.Select(Finite)), at);
            Assert.That(D(Column(root, "replicate_variance")[i]), Is.EqualTo(Finite(e.ReplicateVariance)), at);
            Assert.That(D(Column(root, "residual_variance")[i]), Is.EqualTo(Finite(e.ResidualVariance)), at);
        }
    }

    [Test]
    public void TraitEffects_ThroughTheOccupancyTable_AreMzLibsFitOnTheSameCells()
    {
        string occupancyPath = TraitOccupancy();
        var occupancy = SiteOccupancyCalculator.Calculate(PtmQtlCases.Trait());
        var (trait, replicate, covariates) = TraitMaps();

        JsonElement plain = Run(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", occupancyPath,
            "--runs", Write("runs.tsv", PtmQtlCases.TraitRunsTsv(covariate: false)));
        var expectedPlain = SiteTraitEffects.Fit(occupancy, trait, replicate);
        AssertEffectsAre(plain, expectedPlain);
        int site = Column(plain, "site_key").FindIndex(e => e.GetString()!.StartsWith("TRAIT_SITE:"));
        Assert.That(Column(plain, "status")[site].GetString(), Is.EqualTo("fitted"));
        Assert.That(D(Column(plain, "effect")[site]), Is.EqualTo(0.5).Within(0.1), "mzLib's TraitEffectIsRecoveredOnTheLogitScale");
        Assert.That(D(Column(plain, "df")[site]), Is.EqualTo(4), "trait constant within a replicate: G - 1 - 1");
        Assert.That(Column(plain, "runs")[site].GetInt32(), Is.EqualTo(12));
        Assert.That(Column(plain, "replicates")[site].GetInt32(), Is.EqualTo(6));
        int censored = Column(plain, "site_key").FindIndex(e => e.GetString()!.StartsWith("CENSORED_SITE:"));
        Assert.That(Column(plain, "runs")[censored].GetInt32(), Is.EqualTo(10), "the floor and the ceiling leave the fit");
        Assert.That(Column(plain, "replicates")[censored].GetInt32(), Is.EqualTo(5));
        Assert.That(plain.GetProperty("trait_name").GetString(), Is.EqualTo("dose"));
        Assert.That(plain.GetProperty("covariate_names").GetArrayLength(), Is.EqualTo(0));
        Assert.That(plain.GetProperty("min_runs").GetInt32(), Is.EqualTo(new SiteTraitOptions().MinRuns));
        Assert.That(plain.GetProperty("min_replicates").GetInt32(), Is.EqualTo(new SiteTraitOptions().MinReplicates));
        Assert.That(plain.GetProperty("logit_epsilon").GetDouble(), Is.EqualTo(new SiteTraitOptions().LogitEpsilon));
        Assert.That(plain.GetProperty("exclude_ceiling").GetBoolean(), Is.True);
        Assert.That(plain.GetProperty("cell_count").GetInt32(), Is.EqualTo(occupancy.Count));
        Assert.That(plain.GetProperty("status_counts").GetProperty("fitted").GetInt32(), Is.EqualTo(2));

        JsonElement withCovariate = Run(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", occupancyPath,
            "--runs", Write("runs_cov.tsv", PtmQtlCases.TraitRunsTsv()));
        AssertEffectsAre(withCovariate, SiteTraitEffects.Fit(occupancy, trait, replicate, covariates, new[] { "injection" }));
        Assert.That(D(Column(withCovariate, "covariate_effects")[site][0]), Is.EqualTo(0.02).Within(0.005),
            "mzLib's TraitFitReportsCovariateEffects: the injection offset, +0.01 against -0.01");
        Assert.That(withCovariate.GetProperty("covariate_names").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "injection" }));
    }

    [Test]
    public void TraitEffects_Options_AreSiteTraitOptions()
    {
        string occupancyPath = TraitOccupancy();
        string runs = Write("runs.tsv", PtmQtlCases.TraitRunsTsv());
        var occupancy = SiteOccupancyCalculator.Calculate(PtmQtlCases.Trait());
        var (trait, replicate, covariates) = TraitMaps();
        string[] Call(params string[] extra) => new[] { "ptmqtl", "trait-effects", "--occupancy", occupancyPath, "--runs", runs }.Concat(extra).ToArray();

        JsonElement strict = Run(PtmQtl.TraitEffects, Call("--min-replicates", "6"));
        AssertEffectsAre(strict, SiteTraitEffects.Fit(occupancy, trait, replicate, covariates, new[] { "injection" },
            new SiteTraitOptions { MinReplicates = 6 }));
        Assert.That(Column(strict, "status").Select(e => e.GetString()), Is.EqualTo(new[] { "below_support", "fitted" }));
        Assert.That(strict.GetProperty("min_replicates").GetInt32(), Is.EqualTo(6));
        Assert.That(strict.GetProperty("status_counts").GetProperty("below_support").GetInt32(), Is.EqualTo(1));
        Assert.That(Caveats(strict), Has.Some.StartsWith("1 of 2 sites were not fitted (1 below_support)"));

        JsonElement every = Run(PtmQtl.TraitEffects, Call("--min-runs", "12", "--logit-epsilon", "1e-12", "--keep-ceiling"));
        AssertEffectsAre(every, SiteTraitEffects.Fit(occupancy, trait, replicate, covariates, new[] { "injection" },
            new SiteTraitOptions { MinRuns = 12, LogitEpsilon = 1e-12, ExcludeCeiling = false }));
        Assert.That(every.GetProperty("exclude_ceiling").GetBoolean(), Is.False);
        Assert.That(every.GetProperty("logit_epsilon").GetDouble(), Is.EqualTo(1e-12));

        Assert.That(Usage(PtmQtl.TraitEffects, Call("--logit-epsilon", "0.5")), Does.Contain("above 0 and below 0.5"));
        Assert.That(Usage(PtmQtl.TraitEffects, Call("--logit-epsilon", "0")), Does.Contain("above 0 and below 0.5"));
        Assert.That(Usage(PtmQtl.TraitEffects, Call("--min-runs", "-1")), Does.Contain("zero or more"));
        Assert.That(Usage(PtmQtl.TraitEffects, Call("--min-runs", "six")), Does.Contain("must be an integer"));
        Assert.That(Usage(PtmQtl.TraitEffects, Call("--keep-ceiling", "yes")), Does.Contain("takes no value"));
    }

    [Test]
    public void TraitEffects_Scale_IsSiteTraitOptionsScale()
    {
        string occupancyPath = TraitOccupancy();
        string runs = Write("runs.tsv", PtmQtlCases.TraitRunsTsv());
        var occupancy = SiteOccupancyCalculator.Calculate(PtmQtlCases.Trait());
        var (trait, replicate, covariates) = TraitMaps();
        string[] Call(params string[] extra) => new[] { "ptmqtl", "trait-effects", "--occupancy", occupancyPath, "--runs", runs }.Concat(extra).ToArray();

        JsonElement natural = Run(PtmQtl.TraitEffects, Call());
        Assert.That(natural.GetProperty("scale").GetString(), Is.EqualTo("natural_log_odds"), "absent, mzLib's default");
        Assert.That(Caveats(natural), Has.Some.StartsWith("scale is natural_log_odds"));

        JsonElement log2 = Run(PtmQtl.TraitEffects, Call("--scale", "log2_odds"));
        AssertEffectsAre(log2, SiteTraitEffects.Fit(occupancy, trait, replicate, covariates, new[] { "injection" },
            new SiteTraitOptions { Scale = OccupancyScale.Log2Odds }));
        Assert.That(log2.GetProperty("scale").GetString(), Is.EqualTo("log2_odds"));
        Assert.That(Caveats(log2), Has.Some.StartsWith("scale is log2_odds"));
        // mzLib's TraitFitOnLog2OddsIsTheNaturalLogFitDividedByLn2: a change of unit, not a refit.
        for (int i = 0; i < log2.GetProperty("row_count").GetInt32(); i++)
        {
            Assert.That(D(Column(log2, "effect")[i]), Is.EqualTo(D(Column(natural, "effect")[i]) / Math.Log(2)));
            Assert.That(D(Column(log2, "p_value")[i]), Is.EqualTo(D(Column(natural, "p_value")[i])));
            Assert.That(D(Column(log2, "df")[i]), Is.EqualTo(D(Column(natural, "df")[i])));
        }

        JsonElement explicitNatural = Run(PtmQtl.TraitEffects, Call("--scale", "natural_log_odds"));
        Assert.That(Column(explicitNatural, "effect").Select(D), Is.EqualTo(Column(natural, "effect").Select(D)));

        Assert.That(Usage(PtmQtl.TraitEffects, Call("--scale", "log10_odds")), Does.Contain("must be one of natural_log_odds, log2_odds"));
        Assert.That(Usage(PtmQtl.TraitEffects, Call("--scale")), Does.Contain("without a value"));
    }

    [Test]
    public void TraitEffects_ReadsTheReportedFraction_AndRefusesARepeatedCell()
    {
        // A stored cell whose intensities are rounded: the fraction column wins, as mzLib's ReportedFraction does.
        var header = string.Join('\t', PtmQtl.OccupancyColumns);
        var lines = new StringBuilder(header + "\n");
        var cells = new List<SiteRunOccupancy>();
        var site = new ModificationSite("P1", 13, 'S', PtmQtlCases.Phos);
        foreach (var (run, rep, i) in PtmQtlCases.TraitRuns().Select((t, i) => (t.Run, t.Replicate, i)))
        {
            double fraction = 0.1 + 0.01 * i;
            cells.Add(new SiteRunOccupancy(site, run, OccupancyState.Quantified, 1000, 5000, true) { ReportedFraction = fraction });
            // The peptidoform counts left blank: not supplied, as a stored table may not have them.
            lines.Append($"{site.Key}\tP1\t13\tS\t{PtmQtlCases.Phos}\tCommon Biological\t{run}\tquantified\t1000\t5000\t{PtmQtlCases.R(fraction)}\ttrue\t\t\n");
        }
        string runs = Write("runs.tsv", PtmQtlCases.TraitRunsTsv(covariate: false));
        var (trait, replicate, _) = TraitMaps();
        JsonElement root = Run(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", Write("stored.tsv", lines.ToString()), "--runs", runs);
        AssertEffectsAre(root, SiteTraitEffects.Fit(cells, trait, replicate));
        Assert.That(D(Column(root, "median_fraction")[0]), Is.EqualTo(0.155).Within(1e-12), "the reported fractions, not 1000 / 5000");

        string twice = lines.ToString() + lines.ToString().Split('\n')[1] + "\n";
        Assert.That(Usage(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", Write("twice.tsv", twice), "--runs", runs),
            Does.Contain("appears twice in run R0_0"));
    }

    [Test]
    public void TraitEffects_RefusesMalformedTables_AsUsageErrors()
    {
        string occupancyPath = TraitOccupancy();
        string Refused(string runsBody) =>
            Usage(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", occupancyPath, "--runs", Write($"r{Guid.NewGuid():N}.tsv", runsBody));

        Assert.That(Refused("run\treplicate\nR0_0\tR0\n"), Does.Contain("needs the run, its replicate and the trait"));
        Assert.That(Refused("run\treplicate\tdose\nR0_0\tR0\t1\nR0_0\tR0\t1\n"), Does.Contain("more than once"));
        Assert.That(Refused("run\treplicate\tdose\tdose\nR0_0\tR0\t1\t2\n"), Does.Contain("'dose' appears more than once"));
        Assert.That(Refused("run\treplicate\tdose\nR0_0\tR0\tlow\n"), Does.Contain("is not a number"));
        // mzLib's refusals, in mzLib's words: a trait run with no replicate, a covariate that is not finite.
        Assert.That(Refused(PtmQtlCases.TraitRunsTsv(covariate: false).Replace("R3_1\tR3\t", "R3_1\t\t")), Does.Contain("has no biological replicate"));
        Assert.That(Refused(PtmQtlCases.TraitRunsTsv().Replace("\t2\t1\n", "\t2\t\n")), Does.Contain("finite covariate values"));
        // A trait that does not vary leaves no full-rank design: mzLib refuses it.
        string constant = "run\treplicate\tdose\n" + string.Concat(PtmQtlCases.TraitRuns().Select(t => $"{t.Run}\t{t.Replicate}\t1\n"));
        Assert.That(Refused(constant), Does.Contain("not of full column rank"));

        // A run with a blank trait is not used, as mzLib reads a run absent from the trait.
        string blank = PtmQtlCases.TraitRunsTsv(covariate: false).Replace("R0_0\tR0\t0.25", "R0_0\tR0\t");
        JsonElement root = Run(PtmQtl.TraitEffects, "ptmqtl", "trait-effects", "--occupancy", occupancyPath, "--runs", Write("blank.tsv", blank));
        Assert.That(root.GetProperty("trait_run_count").GetInt32(), Is.EqualTo(11));
        Assert.That(root.GetProperty("run_count").GetInt32(), Is.EqualTo(12));
        Assert.That(root.GetProperty("occupancy_runs_not_used").EnumerateArray().Select(e => e.GetString()), Is.EqualTo(new[] { "R0_0" }));
        Assert.That(Caveats(root), Has.Some.StartsWith("1 run(s) of the occupancy table have no trait value"));

        string occupancyText = File.ReadAllText(occupancyPath);
        string BadOccupancy(string text) => Usage(PtmQtl.TraitEffects, "ptmqtl", "trait-effects",
            "--occupancy", Write($"o{Guid.NewGuid():N}.tsv", text), "--runs", Write("runs_ok.tsv", PtmQtlCases.TraitRunsTsv()));
        Assert.That(BadOccupancy(occupancyText.Replace("\tquantified\t", "\tmeasured\t")), Does.Contain("must be one of quantified, floor, count_only, not_detected"));
        Assert.That(BadOccupancy(occupancyText.Replace("\ttrue\t", "\tyes\t")), Does.Contain("must be true or false"));
        Assert.That(BadOccupancy(ReplaceLast(occupancyText, "\t1\n", "\tone\n")), Does.Contain("modified_peptidoforms must be an integer"));
        Assert.That(BadOccupancy(occupancyText.Replace("unmodified_quantified", "unmodified")), Does.Contain("no column unmodified_quantified"));
    }
}
