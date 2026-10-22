using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// <c>stats combine-p</c>: one p-value per feature from several independent tests of it, by mzLib's
/// <see cref="PValueCombination"/> (mzLib #1367, 1.0.594). The bridge groups the stdin lines by
/// feature and makes one mzLib call per feature; every number is mzLib's.
/// </summary>
internal static partial class Statistics
{
    /// <summary>Row columns of the <c>stats combine-p</c> table, in order.</summary>
    internal static readonly string[] CombineColumns = { "feature", "studies", "omitted", "statistic", "p_value" };

    /// <summary>
    /// <c>stats combine-p --method fisher|stouffer</c>, the tests on stdin as
    /// <c>feature\tp_value[\tweight[\teffect]]</c>, one test per line.
    /// </summary>
    /// <remarks>
    /// <para>
    /// No header line; blank lines are skipped. A feature's lines need not be adjacent: features are
    /// reported in order of first appearance, and each is combined by its own mzLib call.
    /// </para>
    /// <para>
    /// <b>p_value</b>: blank, NA or NaN is a test that could not be run. mzLib omits it and does not
    /// count it in <c>studies</c>; the bridge reports it in <c>omitted</c>.
    /// </para>
    /// <para>
    /// <b>weight</b> (Stouffer only): on every line or on none. A blank or NA third field is no
    /// weight. mzLib needs a weight for every test of a feature, those it omits included.
    /// </para>
    /// <para>
    /// <b>effect</b>: a fourth field on every line or on none. With it, each p-value is read as
    /// TWO-SIDED and converted by <see cref="PValueCombination.OneSided"/> to the one-sided p-value
    /// for effect &gt; 0 before combining, as mzLib's own remarks direct for tests whose effects may
    /// point either way. A blank, NA, NaN or 0 effect has no direction; mzLib returns NaN for it and
    /// the test is omitted.
    /// </para>
    /// </remarks>
    public static object CombineP(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "method");
        // No default: mzLib has two rules and no preferred one, so the caller names it.
        string methodText = arguments.Required("method");
        PValueCombinationMethod method = methodText switch
        {
            "fisher" => PValueCombinationMethod.Fisher,
            "stouffer" => PValueCombinationMethod.Stouffer,
            _ => throw new Program.UsageException(
                $"Option --method must be 'fisher' or 'stouffer'; got '{methodText}'."),
        };

        List<string> lines = Program.ReadStdinLines();
        if (lines.Count == 0)
            throw new Program.UsageException(
                "No tests were given on stdin. Give one per line: feature<TAB>p_value, optionally <TAB>weight and <TAB>effect.");

        var order = new List<string>();
        var byFeature = new Dictionary<string, List<(double P, double? Weight, double Effect)>>(StringComparer.Ordinal);
        bool? signed = null;
        bool? weighted = null;
        for (int i = 0; i < lines.Count; i++)
        {
            string[] parts = lines[i].Split('\t');
            if (parts.Length < 2 || parts.Length > 4)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin has {parts.Length} tab-separated field{(parts.Length == 1 ? "" : "s")}; " +
                    "expected 2 to 4: feature, p_value, then optionally weight and effect.");
            string feature = parts[0].Trim();
            if (feature.Length == 0)
                throw new Program.UsageException($"Line {i + 1} of stdin has an empty feature name.");

            double p = ParseOptional(parts[1], $"Line {i + 1} p_value");

            double? weight = null;
            if (parts.Length >= 3)
            {
                double w = ParseOptional(parts[2], $"Line {i + 1} weight");
                if (!double.IsNaN(w))
                    weight = w;
            }
            bool lineWeighted = weight.HasValue;
            if (weighted is null)
                weighted = lineWeighted;
            else if (weighted != lineWeighted)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin {(lineWeighted ? "gives" : "has no")} weight, but line 1 {(weighted.Value ? "gives one" : "has none")}. " +
                    "Give a weight on every line or on none.");

            bool lineSigned = parts.Length == 4;
            if (signed is null)
                signed = lineSigned;
            else if (signed != lineSigned)
                throw new Program.UsageException(
                    $"Line {i + 1} of stdin {(lineSigned ? "gives" : "has no")} effect field, but line 1 {(signed.Value ? "gives one" : "has none")}. " +
                    "Give an effect on every line (a blank one for a test with no direction) or on none.");
            double effect = lineSigned ? ParseOptional(parts[3], $"Line {i + 1} effect") : double.NaN;

            if (!byFeature.TryGetValue(feature, out var tests))
            {
                tests = new List<(double, double?, double)>();
                byFeature[feature] = tests;
                order.Add(feature);
            }
            tests.Add((p, weight, effect));
        }

        if (weighted == true && method == PValueCombinationMethod.Fisher)
            throw new Program.UsageException(
                "Fisher's method takes no weights (mzLib's PValueCombination.Fisher has none). Use --method stouffer, or leave the weight field blank.");

        var table = new Sdrf.ColumnTable(CombineColumns);
        int omittedTotal = 0, nothingCombined = 0, contradictory = 0, infinite = 0;
        foreach (string feature in order)
        {
            var tests = byFeature[feature];
            CombinedPValue combined;
            try
            {
                double[] p = signed == true
                    ? tests.Select(t => PValueCombination.OneSided(t.P, t.Effect)).ToArray()
                    : tests.Select(t => t.P).ToArray();
                combined = method == PValueCombinationMethod.Fisher
                    ? PValueCombination.Fisher(p)
                    : PValueCombination.Stouffer(p, weighted == true ? tests.Select(t => t.Weight!.Value).ToArray() : null);
            }
            catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
            {
                // A p-value outside [0, 1] or a weight that is not finite and positive: mzLib's
                // refusal of the caller's input, so a usage error naming the feature.
                throw new Program.UsageException(
                    $"Feature '{feature}': {invalid.Message.Split(" (Parameter")[0].Split(Environment.NewLine)[0]}");
            }

            int omitted = tests.Count - combined.Studies;
            omittedTotal += omitted;
            if (combined.Studies == 0) nothingCombined++;
            else if (double.IsNaN(combined.PValue)) contradictory++;
            else if (double.IsInfinity(combined.Statistic)) infinite++;
            table.Add(feature, combined.Studies, omitted, Finite(combined.Statistic), Finite(combined.PValue));
        }

        return new
        {
            method = CombinationName(method),
            weighted = weighted == true,
            signed = signed == true,
            line_count = lines.Count,
            feature_count = order.Count,
            row_count = table.RowCount,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = CombineCaveats(method, signed == true, omittedTotal, nothingCombined, contradictory, infinite),
        };
    }

    /// <summary>
    /// mzLib's <see cref="PValueCombinationMethod"/> as the wire writes it. Every value is named here,
    /// so a value mzLib adds fails a test on the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string CombinationName(PValueCombinationMethod method) => method switch
    {
        PValueCombinationMethod.Fisher => "fisher",
        PValueCombinationMethod.Stouffer => "stouffer",
        _ => throw new InvalidOperationException(
            $"PValueCombinationMethod.{method} has no wire name. Name it in Statistics.CombinationName."),
    };

    /// <summary>A numeric field that may be missing: blank, NA and NaN are NaN; anything else unparseable is refused.</summary>
    private static double ParseOptional(string raw, string what)
    {
        double value = ParseCell(raw, out bool missing);
        if (!missing && double.IsNaN(value))
            throw new Program.UsageException(
                $"{what} is not a number: '{raw.Trim()}'. Leave it blank, or write NA, for a missing value.");
        return value;
    }

    private static List<string> CombineCaveats(PValueCombinationMethod method, bool signed, int omitted,
        int nothingCombined, int contradictory, int infinite)
    {
        var caveats = new List<string>
        {
            "Both rules assume the tests are independent and test the same one-sided hypothesis. Fisher is driven by the smallest p-values; Stouffer weighs every test, and with weights (for example the square root of each study's sample size) lets larger studies count more.",
            "p_value is the upper tail of the statistic, so it tests one direction. Two-sided p-values given without their effects are combined as if one-sided, and opposite effects then reinforce each other instead of cancelling: give each test's effect for that.",
            "Neither rule estimates an effect size; pool the estimates with stats meta for that. p_value is not adjusted across features; send the p_value column to stats adjust.",
        };
        if (signed)
        {
            caveats.Add(
                "Each p-value was read as two-sided and converted with the sign of its effect (mzLib's PValueCombination.OneSided): p/2 for a positive effect, 1 - p/2 for a negative one. p_value therefore tests effect > 0" +
                (method == PValueCombinationMethod.Stouffer ? ", and a p_value near 1 means a consistently negative effect. " : ". ") +
                "A test whose effect is 0, blank or NA has no direction and is omitted.");
            if (method == PValueCombinationMethod.Stouffer)
                caveats.Add(
                    "For a negative effect whose two-sided p-value is below about 2e-16, 1 - p/2 rounds to 1, so that test's z is minus infinity and decides the feature (p_value 1).");
        }
        if (omitted > 0)
            caveats.Add(
                $"{omitted} test{(omitted == 1 ? " was" : "s were")} omitted: a blank, NA or NaN p-value" +
                (signed ? ", or an effect that is 0, blank or NA" : "") + ". studies counts only the tests combined.");
        if (nothingCombined > 0)
            caveats.Add($"{nothingCombined} feature{(nothingCombined == 1 ? " has" : "s have")} no test left to combine; statistic and p_value are null.");
        if (contradictory > 0)
            caveats.Add(
                $"{contradictory} feature{(contradictory == 1 ? " holds" : "s hold")} both a p-value of 0 and a p-value of 1, whose z values are plus and minus infinity; statistic and p_value are null.");
        if (infinite > 0)
            caveats.Add(
                $"{infinite} feature{(infinite == 1 ? " has" : "s have")} an infinite statistic, reported as null: a p-value of 0" +
                (method == PValueCombinationMethod.Stouffer ? " (p_value 0) or of 1 (p_value 1)" : " (p_value 0)") + " among the tests.");
        return caveats;
    }
}
