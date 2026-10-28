using System.Globalization;
using System.Reflection;
using System.Text;
using Quantification.PtmQtl;
using StatisticalModels;

namespace MzLibBridge;

/// <summary>
/// The <c>ptmqtl</c> verbs: intensity-based site occupancy per run, PTM-PTM pairs, and the effect of a
/// per-run trait on a site's occupancy, by mzLib's <c>Quantification.PtmQtl</c> (mzLib #1369, first in
/// 1.0.594; fractions combined per sample #1430 and the run level #1446, first in 1.0.595). The bridge
/// reads plain tab-separated tables into mzLib's own plain records and shapes the answer; every
/// number is mzLib's.
/// </summary>
/// <remarks>
/// <para>
/// <b>The input is mzLib's own observation record, as a table.</b> <see cref="PeptidoformObservation"/>
/// is plain values on purpose ("so any caller, a result file, a stored catalog, a language binding, can
/// supply it"). mzLib has no reader that builds it from MetaMorpheus or FlashLFQ output, and the bridge
/// does not write one: that adapter would decide which PSMs count and how a shared peptide is placed,
/// which is mzLib's to decide (bridge PRINCIPLE.md). It is tracked as a mainland item.
/// </para>
/// <para>
/// <b>One pivot table.</b> <c>ptmqtl occupancy --out</c> writes the occupancy table that
/// <c>ptmqtl trait-effects --occupancy</c> reads, so the two compose through a file with no glue in any
/// binding. <c>ptmqtl pairs</c> takes the observations themselves, because mzLib's
/// <see cref="PtmPairEngine.Physical"/> needs the peptidoforms and <see cref="PtmPairEngine.CoVarying"/>
/// needs them to decide which pairs overlap; it computes the occupancy inside with the same mzLib call.
/// </para>
/// <para>
/// <b>No <c>--threads</c>.</b> None of these mzLib calls takes a degree of parallelism.
/// <see cref="SiteTraitEffects.Fit"/> runs <see cref="MixedModel.Fit"/> at its own default, whose output
/// mzLib documents as independent of the thread count.
/// </para>
/// </remarks>
internal static class PtmQtl
{
    /// <summary>Columns of the observations table, as <see cref="PeptidoformObservation"/> names its members.</summary>
    internal static readonly string[] ObservationColumns =
        { "run", "full_sequence", "protein_accession", "start_residue", "end_residue", "intensity", "previous_residue" };

    /// <summary>Row columns of <c>ptmqtl occupancy</c>, in order; also the header <c>--out</c> writes.</summary>
    internal static readonly string[] OccupancyColumns =
    {
        "site_key", "protein_accession", "position", "residue", "modification", "category", "run", "state",
        "modified_intensity", "covering_intensity", "fraction", "unmodified_quantified",
    };

    /// <summary>Row columns of <c>ptmqtl pairs</c>, in order.</summary>
    internal static readonly string[] PairColumns =
    {
        "result_type", "site_a", "site_b", "same_protein", "overlapping", "statistic", "p_value", "q",
        "fdr_family", "n", "statistic_n", "spearman_method", "run_level_removed",
    };

    /// <summary>Row columns of <c>ptmqtl trait-effects</c>, in order.</summary>
    internal static readonly string[] TraitEffectColumns =
    {
        "site_key", "protein_accession", "position", "residue", "modification", "status", "effect",
        "standard_error", "df", "p_value", "q", "runs", "replicates", "median_fraction", "covariate_effects",
        "replicate_variance", "residual_variance",
    };

    private static readonly string[] PairTypes = { "physical", "covarying", "both" };
    private static readonly string[] RunLevelGroupRules = { "category", "modification" };

    // ---- ptmqtl occupancy ------------------------------------------------------------------------

    /// <summary>
    /// <c>ptmqtl occupancy --observations FILE [--samples FILE] [--modifications NAMES | --categories NAMES] [--out FILE.tsv]</c>.
    /// </summary>
    public static object Occupancy(Program.Arguments arguments)
    {
        string? outputPath = OutputPath(arguments);
        string observationsPath = RequireFile(arguments, "observations");
        if (outputPath is not null && SamePath(outputPath, observationsPath))
            throw new Program.UsageException(
                $"Option --out must differ from --observations: writing to '{observationsPath}' would overwrite the input file.");
        string? samplesPath = OptionalFile(arguments, "samples");
        ModificationFilter filter = ModificationFilter.From(arguments);

        ObservationTable table = ReadObservations(observationsPath);
        Dictionary<string, string>? samples = samplesPath is null ? null : ReadSamples(samplesPath);

        IReadOnlyList<SiteRunOccupancy> occupancy = FromMzLib(() => samples is null
            ? SiteOccupancyCalculator.Calculate(table.Rows, includeModification: filter.Include)
            : SiteOccupancyCalculator.Calculate(table.Rows, samples, filter.Include));

        var columns = new Sdrf.ColumnTable(OccupancyColumns);
        foreach (SiteRunOccupancy cell in occupancy)
            columns.Add(OccupancyRow(cell));

        object? written = null;
        if (outputPath is not null)
            written = WriteTsv(outputPath, columns);

        var stateCounts = Enum.GetValues<OccupancyState>().ToDictionary(StateName, s => occupancy.Count(o => o.State == s));
        return new
        {
            observations_file = Path.GetFullPath(observationsPath),
            samples_file = samplesPath is null ? null : Path.GetFullPath(samplesPath),
            observation_count = table.Rows.Count,
            run_count = table.Rows.Select(o => o.Run).Distinct(StringComparer.Ordinal).Count(),
            sample_count = samples is null ? (int?)null : samples.Values.Distinct(StringComparer.Ordinal).Count(),
            run_names = occupancy.Select(o => o.Run).Distinct(StringComparer.Ordinal).OrderBy(r => r, StringComparer.Ordinal).ToList(),
            site_count = occupancy.Select(o => o.Site).Distinct().Count(),
            state_counts = stateCounts,
            modifications = filter.Modifications,
            categories = filter.Categories,
            ignored_columns = table.IgnoredColumns,
            row_count = columns.RowCount,
            column_names = columns.Names,
            columns = written is null ? columns.Columns() : null,
            written,
            caveats = OccupancyCaveats(samples is not null, table.IgnoredColumns, stateCounts),
        };
    }

    private static object?[] OccupancyRow(SiteRunOccupancy cell) =>
        new object?[]
        {
            cell.Site.Key, cell.Site.ProteinAccession, cell.Site.Position, cell.Site.Residue.ToString(),
            cell.Site.Modification, cell.Site.Category, cell.Run, StateName(cell.State),
            cell.ModifiedIntensity, cell.CoveringIntensity, Finite(cell.Fraction), cell.UnmodifiedQuantified,
        };

    /// <summary>
    /// mzLib's <see cref="OccupancyState"/> as the wire writes it. Every value is named here, so a value
    /// mzLib adds fails a test on the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string StateName(OccupancyState state) => state switch
    {
        OccupancyState.Quantified => "quantified",
        OccupancyState.Floor => "floor",
        OccupancyState.CountOnly => "count_only",
        OccupancyState.NotDetected => "not_detected",
        _ => throw new InvalidOperationException($"OccupancyState.{state} has no wire name. Name it in PtmQtl.StateName."),
    };

    private static List<string> OccupancyCaveats(bool combined, List<string> ignored, Dictionary<string, int> stateCounts)
    {
        var caveats = new List<string>
        {
            "Intensity occupancy (mzLib's SiteOccupancyCalculator, QuantProject DEF-OCC-INT, the estimator MetaMorpheus writes as IntensityOccupancy_): for a site in a run, the summed intensity of the quantified peptidoforms carrying the modification there, over the summed intensity of every quantified peptidoform covering the position. Supply MS/MS identifications only: a match-between-runs transfer must not be a row.",
            "fraction is null unless state is quantified. floor (the modified form was identified but not quantified while the site was covered) and not_detected (covered, modified form not identified) are censored, never 0; count_only means nothing covering the site was quantified.",
            "A fraction of 1 with unmodified_quantified false is a ceiling: the site was seen only modified, and the unmodified form may be present below detection.",
            "position is 1-based in the protein, and a protein N-terminal modification sits at the residue it is on (1, or 2 after the initiator Met is removed). readers read-occupancy counts the other way (0 is the N-terminus), so do not join the two tables on position.",
            "Excluded by mzLib, whatever --modifications or --categories say: the categories Common Variable and Common Fixed, peptide-terminal modifications, and every C-terminal modification (the protein's length is not supplied). An N-terminal modification is a site only on a peptide starting at residue 1, or at 2 when previous_residue is M or blank; whether it is protein- or peptide-terminal comes from mzLib's modification registry, and a name the registry does not know is treated as protein N-terminal.",
            "A peptide shared by several proteins is one row per protein, and each protein gets its full intensity. MetaMorpheus instead gives a shared peptide a proportional share per protein group, so its IntensityOccupancy_ can differ.",
        };
        if (combined)
            caveats.Add(
                "With --samples, each peptidoform's intensity is summed over a sample's runs before the ratio (mzLib's CombineObservations), so run holds the sample. Map only fractions of one sample together: technical replicates are separate measurements and keep their own key.");
        if (ignored.Count > 0)
            caveats.Add(
                $"{ignored.Count} column(s) of the observations file were not read, because mzLib's observation has no such field: {string.Join(", ", ignored)}. Check the spelling if one was meant to be previous_residue.");
        int censored = stateCounts[StateName(OccupancyState.Floor)] + stateCounts[StateName(OccupancyState.NotDetected)];
        if (censored > 0)
            caveats.Add($"{censored} cells are floor or not_detected: censored values with no fraction. Averaging them in as 0 would bias occupancy down.");
        return caveats;
    }

    // ---- ptmqtl pairs ----------------------------------------------------------------------------

    /// <summary>
    /// <c>ptmqtl pairs --observations FILE --type physical|covarying|both [--samples FILE]
    /// [--modifications NAMES | --categories NAMES] [--min-quantified-fraction F] [--keep-ceiling]
    /// [--remove-run-level [--run-level-group category|modification | --run-level-groups FILE]]</c>.
    /// </summary>
    public static object Pairs(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "type");
        string type = arguments.Required("type").Trim();
        if (!PairTypes.Contains(type, StringComparer.Ordinal))
            throw new Program.UsageException(
                $"Option --type must be one of {string.Join(", ", PairTypes)}; got '{type}'. physical: two modifications on one peptidoform; covarying: two sites' occupancies across runs.");
        bool physical = type is "physical" or "both";
        bool covarying = type is "covarying" or "both";

        string observationsPath = RequireFile(arguments, "observations");
        string? samplesPath = OptionalFile(arguments, "samples");
        ModificationFilter filter = ModificationFilter.From(arguments);

        bool fractionGiven = arguments.WasProvided("min-quantified-fraction");
        bool keepCeiling = Flag(arguments, "keep-ceiling");
        bool removeRunLevel = Flag(arguments, "remove-run-level");
        bool groupGiven = arguments.WasProvided("run-level-group");
        bool groupsGiven = arguments.WasProvided("run-level-groups");
        if (!covarying)
        {
            string? stray = new[] { ("min-quantified-fraction", fractionGiven), ("keep-ceiling", keepCeiling),
                    ("remove-run-level", removeRunLevel), ("run-level-group", groupGiven), ("run-level-groups", groupsGiven) }
                .Where(o => o.Item2).Select(o => o.Item1).FirstOrDefault();
            if (stray is not null)
                throw new Program.UsageException(
                    $"Option --{stray} shapes co-varying pairs only; use --type covarying or both, or omit it.");
        }
        if ((groupGiven || groupsGiven) && !removeRunLevel)
            throw new Program.UsageException(
                "Options --run-level-group and --run-level-groups choose whose run level is removed; give --remove-run-level as well, or omit them.");
        if (groupGiven && groupsGiven)
            throw new Program.UsageException(
                "Give --run-level-group or --run-level-groups, not both: one rule says which group a site's run level comes from.");

        double? minFraction = null;
        if (fractionGiven)
        {
            minFraction = arguments.OptionalDouble("min-quantified-fraction", double.NaN);
            // mzLib's own condition (PtmPairEngine.CoVarying), checked here so a caller's mistake is a
            // usage error naming the option, not mzLib's ArgumentOutOfRangeException.
            if (!(minFraction > 0 && minFraction <= 1))
                throw new Program.UsageException(
                    $"Option --min-quantified-fraction must be above 0 and at most 1; got {minFraction.Value.ToString(CultureInfo.InvariantCulture)}.");
        }
        string? groupRule = null;
        if (groupGiven)
        {
            RequireValueIfProvided(arguments, "run-level-group");
            groupRule = arguments.Optional("run-level-group")!.Trim();
            if (!RunLevelGroupRules.Contains(groupRule, StringComparer.Ordinal))
                throw new Program.UsageException(
                    $"Option --run-level-group must be one of {string.Join(", ", RunLevelGroupRules)}; got '{groupRule}'. For any other grouping, map each modification to a group in a file and pass --run-level-groups.");
        }
        string? groupsPath = groupsGiven ? RequireFile(arguments, "run-level-groups") : null;

        ObservationTable table = ReadObservations(observationsPath);
        Dictionary<string, string>? samples = samplesPath is null ? null : ReadSamples(samplesPath);
        Dictionary<string, string>? groups = groupsPath is null ? null : ReadModificationGroups(groupsPath);

        // As mzLib's CombineObservations says: feed the combined observations to Physical and to Calculate,
        // so both see samples rather than fractions.
        IReadOnlyList<PeptidoformObservation> observations = samples is null
            ? table.Rows
            : FromMzLib(() => SiteOccupancyCalculator.CombineObservations(table.Rows, samples));

        var pairs = new List<PtmPair>();
        int? physicalCount = null, covaryingCount = null, covaryingRunCount = null;
        bool? excludeCeiling = null;
        // Absent, mzLib's own default (false at the pin), read off the method like the others.
        bool remove = covarying && (removeRunLevel || CoVaryingDefault<bool>("removeRunLevel"));
        if (physical)
        {
            IReadOnlyList<PtmPair> found = FromMzLib(() => PtmPairEngine.Physical(observations, filter.Include));
            physicalCount = found.Count;
            pairs.AddRange(found);
        }
        if (covarying)
        {
            IReadOnlyList<SiteRunOccupancy> occupancy =
                FromMzLib(() => SiteOccupancyCalculator.Calculate(observations, includeModification: filter.Include));
            if (groups is not null)
            {
                List<string> unmapped = occupancy.Select(o => o.Site.Modification).Distinct(StringComparer.Ordinal)
                    .Where(m => !groups.ContainsKey(m)).OrderBy(m => m, StringComparer.Ordinal).ToList();
                if (unmapped.Count > 0)
                    throw new Program.UsageException(
                        $"The --run-level-groups file gives no group for {unmapped.Count} modification(s) of these sites: {string.Join("; ", unmapped)}. Give every modification a group.");
            }
            Func<ModificationSite, string>? runLevelGroup = groupRule switch
            {
                "category" => site => site.Category,
                "modification" => site => site.Modification,
                _ => groups is null ? null : site => groups[site.Modification],
            };
            // Absent options take mzLib's own defaults, read off the method, so the bridge sets none.
            minFraction ??= CoVaryingDefault<double>("minQuantifiedFraction");
            // mzLib's excludeCeiling is nullable and resolves its own default inside CoVarying, so absent
            // stays null on the wire too: the bridge does not restate the rule to name the value.
            excludeCeiling = keepCeiling ? false : null;
            IReadOnlyList<PtmPair> found = FromMzLib(() => PtmPairEngine.CoVarying(occupancy, observations,
                minFraction.Value, excludeCeiling, remove, runLevelGroup));
            covaryingCount = found.Count;
            covaryingRunCount = occupancy.Select(o => o.Run).Distinct(StringComparer.Ordinal).Count();
            pairs.AddRange(found);
        }

        var columns = new Sdrf.ColumnTable(PairColumns);
        foreach (PtmPair pair in pairs)
            columns.Add(
                PairTypeName(pair.ResultType), pair.SiteA.Key, pair.SiteB.Key, pair.SameProtein, pair.Overlapping,
                Finite(pair.Statistic), Finite(pair.PValue), Finite(pair.Q), pair.FdrFamily, pair.N, pair.StatisticN,
                pair.SpearmanMethod is SpearmanPValueMethod method ? Statistics.SpearmanMethodName(method) : null,
                pair.RunLevelRemoved);

        return new
        {
            observations_file = Path.GetFullPath(observationsPath),
            samples_file = samplesPath is null ? null : Path.GetFullPath(samplesPath),
            type,
            observation_count = table.Rows.Count,
            run_count = table.Rows.Select(o => o.Run).Distinct(StringComparer.Ordinal).Count(),
            sample_count = samples is null ? (int?)null : samples.Values.Distinct(StringComparer.Ordinal).Count(),
            run_names = observations.Select(o => o.Run).Distinct(StringComparer.Ordinal).OrderBy(r => r, StringComparer.Ordinal).ToList(),
            covarying_run_count = covaryingRunCount,
            min_quantified_fraction = covarying ? minFraction : null,
            exclude_ceiling = excludeCeiling,
            remove_run_level = covarying ? remove : (bool?)null,
            run_level_group = !remove ? null : groupRule ?? (groups is null ? "none" : "file"),
            run_level_groups_file = groupsPath is null ? null : Path.GetFullPath(groupsPath),
            min_runs_for_run_level = remove ? PtmPairEngine.MinRunsForRunLevel : (int?)null,
            min_sites_for_run_level = remove ? PtmPairEngine.MinSitesForRunLevel : (int?)null,
            modifications = filter.Modifications,
            categories = filter.Categories,
            ignored_columns = table.IgnoredColumns,
            physical_count = physicalCount,
            covarying_count = covaryingCount,
            row_count = columns.RowCount,
            column_names = columns.Names,
            columns = columns.Columns(),
            caveats = PairCaveats(physical, covarying, remove, samples is not null, table.IgnoredColumns, pairs),
        };
    }

    /// <summary>mzLib's <see cref="PairResultType"/> as the wire writes it.</summary>
    internal static string PairTypeName(PairResultType type) => type switch
    {
        PairResultType.P => "physical",
        PairResultType.A => "covarying",
        _ => throw new InvalidOperationException($"PairResultType.{type} has no wire name. Name it in PtmQtl.PairTypeName."),
    };

    /// <summary>
    /// The default mzLib declares for a parameter of <see cref="PtmPairEngine.CoVarying"/>, so an option the
    /// caller leaves out is mzLib's value at the pin, read rather than retyped.
    /// </summary>
    internal static T CoVaryingDefault<T>(string parameter) =>
        (T)typeof(PtmPairEngine).GetMethod(nameof(PtmPairEngine.CoVarying))!.GetParameters()
            .Single(p => p.Name == parameter).DefaultValue!;

    private static List<string> PairCaveats(bool physical, bool covarying, bool removeRunLevel, bool combined,
        List<string> ignored, List<PtmPair> pairs)
    {
        var caveats = new List<string>();
        if (physical)
            caveats.Add(
                "physical (mzLib's PtmPairEngine.Physical, type P): both modifications were identified on one peptidoform, so on one molecule. statistic is the median over runs of the co-occupancy, the intensity of peptidoforms carrying both over the intensity of peptidoforms covering both positions, on the statistic_n runs where both were quantified; n counts the runs where the doubly modified form was identified. There is no test: p_value and q are null.");
        if (covarying)
        {
            caveats.Add(
                "covarying (mzLib's PtmPairEngine.CoVarying, type A): Spearman's rho of two sites' occupancies across runs, on the runs where both are quantified. Floors and not_detected cells are not values and are left out, and so, unless --keep-ceiling, are ceilings: two unrelated sites both read 1 in low-load runs, which would correlate them through detection alone.");
            caveats.Add(
                "A site enters co-varying pairs only if it is quantified in at least min_quantified_fraction of covarying_run_count, the runs (or samples) where any site was identified or covered.");
            caveats.Add(
                "q is Benjamini-Hochberg within fdr_family (the type, then intra or inter protein). An overlapping pair, two sites one peptide covers, is reported but left out of the family (q null): shared peptides correlate it by construction. q is not a target-decoy q-value.");
            caveats.Add(
                "p_value is Spearman's two-sided p, exact for 9 or fewer runs without ties, otherwise asymptotic (t on n - 2 degrees of freedom), as stats correlate reports it.");
        }
        if (removeRunLevel)
            caveats.Add(
                $"With --remove-run-level, rho and p_value are the partial Spearman correlation given each run's level: the median occupancy, in that run, of the other tested sites of the pair's group (one group, or one per --run-level-group). A run whose group has fewer than {PtmPairEngine.MinSitesForRunLevel} other sites there has no level, and a pair sharing fewer than {PtmPairEngine.MinRunsForRunLevel} runs with a level is not_estimable. Use it for a batch or loading shift only: where the run level follows the design, removing it tests the correlation within conditions instead, which is a different question.");
        if (combined)
            caveats.Add(
                "With --samples, each peptidoform's intensity is summed over a sample's runs first (mzLib's CombineObservations), so pairs are computed across samples. Map only fractions of one sample together, never technical replicates.");
        caveats.Add(
            "The observations follow ptmqtl occupancy's rules: MS/MS identifications only, one row per (run, peptidoform, protein), and the same exclusions (Common Variable, Common Fixed, peptide-terminal and C-terminal modifications).");
        if (ignored.Count > 0)
            caveats.Add(
                $"{ignored.Count} column(s) of the observations file were not read, because mzLib's observation has no such field: {string.Join(", ", ignored)}.");
        int notEstimable = pairs.Count(p => p.SpearmanMethod == SpearmanPValueMethod.NotEstimable);
        if (notEstimable > 0)
            caveats.Add($"{notEstimable} co-varying pairs are not_estimable: their statistic, p_value and q are null.");
        return caveats;
    }

    // ---- ptmqtl trait-effects --------------------------------------------------------------------

    /// <summary>
    /// <c>ptmqtl trait-effects --occupancy FILE --runs FILE [--min-runs N] [--min-replicates N]
    /// [--logit-epsilon E] [--keep-ceiling]</c>.
    /// </summary>
    public static object TraitEffects(Program.Arguments arguments)
    {
        string occupancyPath = RequireFile(arguments, "occupancy");
        string runsPath = RequireFile(arguments, "runs");

        int? minRuns = OptionalCount(arguments, "min-runs");
        int? minReplicates = OptionalCount(arguments, "min-replicates");
        double? epsilon = null;
        if (arguments.WasProvided("logit-epsilon"))
        {
            epsilon = arguments.OptionalDouble("logit-epsilon", double.NaN);
            // mzLib clamps each fraction to [e, 1 - e] (Math.Clamp), which needs e below 0.5; at 0 a
            // fraction of 1 has no logit.
            if (!(epsilon > 0 && epsilon < 0.5))
                throw new Program.UsageException(
                    $"Option --logit-epsilon must be above 0 and below 0.5; got {epsilon.Value.ToString(CultureInfo.InvariantCulture)}.");
        }
        bool keepCeiling = Flag(arguments, "keep-ceiling");

        // Only what the caller gave is set; the rest stays at mzLib's SiteTraitOptions defaults, and with
        // nothing given no options object is passed at all.
        SiteTraitOptions? options = null;
        if (minRuns is not null || minReplicates is not null || epsilon is not null || keepCeiling)
        {
            options = new SiteTraitOptions();
            if (minRuns is int r) options = options with { MinRuns = r };
            if (minReplicates is int g) options = options with { MinReplicates = g };
            if (epsilon is double e) options = options with { LogitEpsilon = e };
            if (keepCeiling) options = options with { ExcludeCeiling = false };
        }
        SiteTraitOptions used = options ?? new SiteTraitOptions();

        OccupancyTable cells = ReadOccupancy(occupancyPath);
        RunTable runs = ReadRuns(runsPath);

        IReadOnlyList<SiteTraitEffect> effects = FromMzLib(() => SiteTraitEffects.Fit(cells.Cells, runs.Trait, runs.Replicate,
            runs.CovariateNames.Count > 0 ? runs.Covariates : null,
            runs.CovariateNames.Count > 0 ? runs.CovariateNames : null, options));

        var columns = new Sdrf.ColumnTable(TraitEffectColumns);
        foreach (SiteTraitEffect effect in effects)
            columns.Add(
                effect.Site.Key, effect.Site.ProteinAccession, effect.Site.Position, effect.Site.Residue.ToString(),
                effect.Site.Modification, TraitStatusName(effect.Status), Finite(effect.Effect), Finite(effect.StandardError),
                Finite(effect.DegreesOfFreedom), Finite(effect.PValue), Finite(effect.Q), effect.Runs, effect.Replicates,
                Finite(effect.MedianFraction), effect.CovariateEffects.Select(Finite).ToList(),
                Finite(effect.ReplicateVariance), Finite(effect.ResidualVariance));

        var statusCounts = TraitStatuses.ToDictionary(s => s, s => 0);
        foreach (SiteTraitEffect effect in effects)
            statusCounts[TraitStatusName(effect.Status)]++;
        List<string> notUsed = cells.Cells.Select(c => c.Run).Distinct(StringComparer.Ordinal)
            .Where(run => !runs.Trait.ContainsKey(run)).OrderBy(run => run, StringComparer.Ordinal).ToList();

        return new
        {
            occupancy_file = Path.GetFullPath(occupancyPath),
            runs_file = Path.GetFullPath(runsPath),
            trait_name = runs.TraitName,
            covariate_names = runs.CovariateNames,
            cell_count = cells.Cells.Count,
            run_count = runs.RunCount,
            trait_run_count = runs.Trait.Count,
            occupancy_runs_not_used = notUsed,
            min_runs = used.MinRuns,
            min_replicates = used.MinReplicates,
            logit_epsilon = used.LogitEpsilon,
            exclude_ceiling = used.ExcludeCeiling,
            status_counts = statusCounts,
            ignored_columns = cells.IgnoredColumns,
            row_count = columns.RowCount,
            column_names = columns.Names,
            columns = columns.Columns(),
            caveats = TraitCaveats(used, statusCounts, effects.Count, notUsed, cells.IgnoredColumns),
        };
    }

    /// <summary>Every value of the trait-effects status column, in the order status_counts lists them.</summary>
    internal static readonly string[] TraitStatuses =
    {
        "fitted", "below_support", "too_few_observations", "rank_deficient", "too_few_groups",
        "variance_ratio_at_limit", "not_converged",
    };

    /// <summary>
    /// mzLib's <see cref="SiteTraitEffect.Status"/> text as the wire writes it: <c>Fitted</c>,
    /// <c>BelowSupport</c>, or the <see cref="FeatureFitStatus"/> name of a mixed-model fit that could not be
    /// computed, in the snake_case <c>stats fit</c> uses for those names. A text mzLib adds fails a test on
    /// the repin instead of reaching the wire under its C# name.
    /// </summary>
    internal static string TraitStatusName(string status)
    {
        string name = status switch
        {
            "Fitted" => "fitted",
            "BelowSupport" => "below_support",
            _ => Enum.TryParse(status, out FeatureFitStatus fit) && fit != FeatureFitStatus.Fitted
                ? Statistics.StatusName(fit)
                : throw new InvalidOperationException($"SiteTraitEffect status '{status}' has no wire name. Name it in PtmQtl.TraitStatusName."),
        };
        return TraitStatuses.Contains(name)
            ? name
            : throw new InvalidOperationException($"SiteTraitEffect status '{status}' is not one MixedModel reports. Add it to PtmQtl.TraitStatuses.");
    }

    private static List<string> TraitCaveats(SiteTraitOptions used, Dictionary<string, int> statusCounts, int siteCount,
        List<string> notUsed, List<string> ignored)
    {
        var caveats = new List<string>
        {
            "Per site, mzLib's SiteTraitEffects fits logit(occupancy) = b0 + b_trait * trait + sum of b_k * covariate_k + u_replicate + e, a random-intercept linear mixed model by REML (mzLib's MixedModel), so repeated runs of one biological replicate are not counted as independent. effect is b_trait, tested by a Wald t with nlme's containment degrees of freedom: a trait constant within each replicate is tested at the between-replicate df.",
            "Only quantified cells are responses. floor, count_only and not_detected cells are left out, never entered as 0: this is the intensity part of the question, and it says nothing about whether detecting the modification changes with the trait.",
            used.ExcludeCeiling
                ? "Ceiling cells (a fraction of 1 with unmodified_quantified false) are left out (mzLib's SiteTraitOptions.ExcludeCeiling); --keep-ceiling keeps them."
                : "With --keep-ceiling, ceiling cells (a fraction of 1 seen only modified) are responses, clamped below 1 before the logit.",
            $"Each fraction is clamped to [logit_epsilon, 1 - logit_epsilon] before the logit, so a fraction of 0 or 1 has a finite logit. effect is in logit occupancy per unit of the trait as given: rescale the trait to change the unit.",
            $"A site is fitted only with at least min_runs runs and min_replicates distinct replicates contributing a response; otherwise it is below_support with null statistics. Sites with no usable response in a run of the runs table do not appear.",
            "q is Benjamini-Hochberg over the sites with a p-value in this call. It is not a target-decoy q-value, and it is not adjusted across datasets.",
        };
        if (notUsed.Count > 0)
            caveats.Add(
                $"{notUsed.Count} run(s) of the occupancy table have no trait value in the runs table and were not used: {string.Join(", ", notUsed.Take(10))}{(notUsed.Count > 10 ? ", ..." : "")}.");
        int notFitted = siteCount - statusCounts["fitted"];
        if (notFitted > 0)
            caveats.Add(
                $"{notFitted} of {siteCount} sites were not fitted ({string.Join(", ", statusCounts.Where(kv => kv.Key != "fitted" && kv.Value > 0).Select(kv => $"{kv.Value} {kv.Key}"))}); their statistics are null and they are not in the Benjamini-Hochberg family.");
        if (ignored.Count > 0)
            caveats.Add($"{ignored.Count} column(s) of the occupancy file were not read: {string.Join(", ", ignored)}.");
        return caveats;
    }

    // ---- input tables ----------------------------------------------------------------------------

    private sealed record ObservationTable(List<PeptidoformObservation> Rows, List<string> IgnoredColumns);

    private sealed record OccupancyTable(List<SiteRunOccupancy> Cells, List<string> IgnoredColumns);

    private sealed record RunTable(int RunCount, string TraitName, List<string> CovariateNames,
        Dictionary<string, double> Trait, Dictionary<string, string> Replicate, Dictionary<string, double[]> Covariates);

    /// <summary>
    /// The observations: a TSV with a header row naming its columns, in any order, as
    /// <see cref="ObservationColumns"/>. <c>previous_residue</c> may be left out; any other column is
    /// not read and is reported. A blank, NA or NaN intensity is NaN, mzLib's "identified, not quantified".
    /// </summary>
    private static ObservationTable ReadObservations(string path)
    {
        (List<string> header, List<string[]> rows) = ReadNamedTsv(path, "observations");
        int Column(string name) => header.IndexOf(name);
        string[] required = ObservationColumns.Where(c => c != "previous_residue").ToArray();
        RequireColumns(header, required, "observations", path);
        int run = Column("run"), sequence = Column("full_sequence"), protein = Column("protein_accession"),
            start = Column("start_residue"), end = Column("end_residue"), intensity = Column("intensity"),
            previous = Column("previous_residue");

        var observations = new List<PeptidoformObservation>(rows.Count);
        for (int r = 0; r < rows.Count; r++)
        {
            string[] row = rows[r];
            string where = $"Line {r + 2} of the observations file";
            observations.Add(new PeptidoformObservation(
                row[run].Trim(), row[sequence].Trim(), row[protein].Trim(),
                ParseInt(row[start], $"{where}, start_residue"), ParseInt(row[end], $"{where}, end_residue"),
                ParseOptional(row[intensity], $"{where}, intensity"),
                previous < 0 ? null : ParseResidue(row[previous], $"{where}, previous_residue", optional: true)));
        }
        List<string> ignored = header.Where(h => !ObservationColumns.Contains(h, StringComparer.Ordinal)).ToList();
        return new ObservationTable(observations, ignored);
    }

    /// <summary>The columns <c>trait-effects</c> reads from an occupancy table; <c>fraction</c> is optional.</summary>
    private static readonly string[] OccupancyInputColumns =
    {
        "protein_accession", "position", "residue", "modification", "run", "state", "modified_intensity",
        "covering_intensity", "unmodified_quantified",
    };

    /// <summary>
    /// An occupancy table as <c>ptmqtl occupancy --out</c> writes it, read back into mzLib's
    /// <see cref="SiteRunOccupancy"/>. <c>site_key</c> and <c>category</c> are derived by mzLib and not
    /// read; a <c>fraction</c> cell becomes <see cref="SiteRunOccupancy.ReportedFraction"/>, mzLib's slot for
    /// a fraction its producer wrote exactly. A (site, run) given twice is refused: mzLib's trait fit would
    /// silently use the first.
    /// </summary>
    private static OccupancyTable ReadOccupancy(string path)
    {
        (List<string> header, List<string[]> rows) = ReadNamedTsv(path, "occupancy");
        RequireColumns(header, OccupancyInputColumns, "occupancy", path);
        int Column(string name) => header.IndexOf(name);
        int protein = Column("protein_accession"), position = Column("position"), residue = Column("residue"),
            modification = Column("modification"), run = Column("run"), state = Column("state"),
            modified = Column("modified_intensity"), covering = Column("covering_intensity"),
            unmodified = Column("unmodified_quantified"), fraction = Column("fraction");

        var cells = new List<SiteRunOccupancy>(rows.Count);
        var seen = new HashSet<(ModificationSite, string)>();
        for (int r = 0; r < rows.Count; r++)
        {
            string[] row = rows[r];
            string where = $"Line {r + 2} of the occupancy file";
            var site = new ModificationSite(row[protein].Trim(), ParseInt(row[position], $"{where}, position"),
                ParseResidue(row[residue], $"{where}, residue", optional: false).GetValueOrDefault(), row[modification].Trim());
            string runName = row[run].Trim();
            if (site.ProteinAccession.Length == 0 || site.Modification.Length == 0 || runName.Length == 0)
                throw new Program.UsageException($"{where} needs a protein_accession, a modification and a run.");
            if (!seen.Add((site, runName)))
                throw new Program.UsageException(
                    $"{where}: {site.Key} appears twice in run {runName}. Supply one row per (site, run), as ptmqtl occupancy writes it.");
            double? reported = fraction < 0 ? null : NullIfMissing(ParseOptional(row[fraction], $"{where}, fraction"));
            cells.Add(new SiteRunOccupancy(site, runName, ParseState(row[state], $"{where}, state"),
                ParseRequired(row[modified], $"{where}, modified_intensity"), ParseRequired(row[covering], $"{where}, covering_intensity"),
                ParseBool(row[unmodified], $"{where}, unmodified_quantified")) { ReportedFraction = reported });
        }
        var known = OccupancyColumns.ToHashSet(StringComparer.Ordinal);
        return new OccupancyTable(cells, header.Where(h => !known.Contains(h)).ToList());
    }

    /// <summary>
    /// The runs table: a header row, then one row per run, the run first, its biological replicate second,
    /// the trait third (its header names the trait), and any further columns covariates (their headers name
    /// them). A run with a blank trait is left out of the trait, which mzLib reads as "not used".
    /// </summary>
    private static RunTable ReadRuns(string path)
    {
        List<string[]> rows = ReadTsv(path, "runs");
        string[] header = rows[0].Select(h => h.Trim()).ToArray();
        if (header.Length < 3)
            throw new Program.UsageException(
                $"The runs file '{path}' has {header.Length} column(s); it needs the run, its replicate and the trait, then any covariates.");
        string traitName = header[2];
        List<string> covariateNames = header.Skip(3).ToList();
        RequireUnique(header.Skip(2).ToList(), "column", "runs header");

        var trait = new Dictionary<string, double>(StringComparer.Ordinal);
        var replicate = new Dictionary<string, string>(StringComparer.Ordinal);
        var covariates = new Dictionary<string, double[]>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            string where = $"Line {r + 1} of the runs file";
            if (row.Length != header.Length)
                throw new Program.UsageException($"{where} has {row.Length} fields; the header has {header.Length}.");
            string run = row[0].Trim();
            if (run.Length == 0)
                throw new Program.UsageException($"{where} has an empty run name.");
            if (replicate.ContainsKey(run))
                throw new Program.UsageException($"Run '{run}' appears more than once in the runs file.");
            replicate[run] = row[1].Trim();
            double value = ParseOptional(row[2], $"{where}, {traitName}");
            if (!double.IsNaN(value))
                trait[run] = value;
            covariates[run] = covariateNames.Select((name, c) => ParseOptional(row[3 + c], $"{where}, {name}")).ToArray();
        }
        if (replicate.Count == 0)
            throw new Program.UsageException($"The runs file '{path}' has a header but no runs.");
        return new RunTable(replicate.Count, traitName, covariateNames, trait, replicate, covariates);
    }

    /// <summary>A run-to-sample map: a header row, then <c>run&lt;TAB&gt;sample</c> per run.</summary>
    private static Dictionary<string, string> ReadSamples(string path) => ReadPairs(path, "samples", "run", "sample");

    /// <summary>A modification-to-group map: a header row, then <c>modification&lt;TAB&gt;group</c>.</summary>
    private static Dictionary<string, string> ReadModificationGroups(string path) =>
        ReadPairs(path, "run-level-groups", "modification", "group");

    private static Dictionary<string, string> ReadPairs(string path, string what, string key, string value)
    {
        List<string[]> rows = ReadTsv(path, what);
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        for (int r = 1; r < rows.Count; r++)
        {
            string[] row = rows[r];
            if (row.Length != 2)
                throw new Program.UsageException(
                    $"Line {r + 1} of the --{what} file has {row.Length} fields; it needs two: {key}, then {value}.");
            string k = row[0].Trim(), v = row[1].Trim();
            if (k.Length == 0 || v.Length == 0)
                throw new Program.UsageException($"Line {r + 1} of the --{what} file has an empty {(k.Length == 0 ? key : value)}.");
            if (map.ContainsKey(k))
                throw new Program.UsageException($"The {key} '{k}' appears more than once in the --{what} file.");
            map[k] = v;
        }
        if (map.Count == 0)
            throw new Program.UsageException($"The --{what} file '{path}' has a header but no rows.");
        return map;
    }

    /// <summary>
    /// <c>--modifications</c> or <c>--categories</c>, the wire form of mzLib's <c>includeModification</c>
    /// filter: exact modification names (<c>Category:IdWithMotif</c>), or categories as mzLib's
    /// <see cref="ModificationSite.Category"/> reads them, separated by tabs (a name can hold spaces and
    /// commas, never a tab). Neither: no filter, as mzLib's default.
    /// </summary>
    private sealed record ModificationFilter(List<string>? Modifications, List<string>? Categories)
    {
        public Func<string, bool>? Include => Modifications is { } names
            ? names.ToHashSet(StringComparer.Ordinal).Contains
            : Categories is { } categories
                ? name => categories.Contains(new ModificationSite("", 0, ' ', name).Category, StringComparer.Ordinal)
                : null;

        public static ModificationFilter From(Program.Arguments arguments)
        {
            bool modifications = arguments.WasProvided("modifications"), categories = arguments.WasProvided("categories");
            if (modifications && categories)
                throw new Program.UsageException(
                    "Give --modifications or --categories, not both: mzLib takes one filter on the modification name.");
            return new ModificationFilter(modifications ? Names(arguments, "modifications") : null,
                categories ? Names(arguments, "categories") : null);
        }

        private static List<string> Names(Program.Arguments arguments, string option)
        {
            RequireValueIfProvided(arguments, option);
            List<string> names = arguments.Optional(option)!.Split('\t').Select(n => n.Trim()).ToList();
            if (names.Any(n => n.Length == 0))
                throw new Program.UsageException($"Option --{option} holds a blank name; separate names with a tab and leave none empty.");
            return names;
        }
    }

    // ---- parsing helpers -------------------------------------------------------------------------

    /// <summary>A TSV whose header names its columns: the trimmed header, and the data rows checked for width.</summary>
    private static (List<string> Header, List<string[]> Rows) ReadNamedTsv(string path, string what)
    {
        List<string[]> rows = ReadTsv(path, what);
        List<string> header = rows[0].Select(h => h.Trim()).ToList();
        RequireUnique(header, "column", $"{what} header");
        for (int r = 1; r < rows.Count; r++)
            if (rows[r].Length != header.Count)
                throw new Program.UsageException(
                    $"Line {r + 1} of the {what} file has {rows[r].Length} fields; the header has {header.Count}.");
        if (rows.Count == 1)
            throw new Program.UsageException($"The {what} file '{path}' has a header but no rows.");
        return (header, rows.Skip(1).ToList());
    }

    private static void RequireColumns(List<string> header, IEnumerable<string> required, string what, string path)
    {
        List<string> missing = required.Where(c => !header.Contains(c, StringComparer.Ordinal)).ToList();
        if (missing.Count > 0)
            throw new Program.UsageException(
                $"The {what} file '{path}' has no column {string.Join(", ", missing)}. Its header has: {string.Join(", ", header)}.");
    }

    private static List<string[]> ReadTsv(string path, string what)
    {
        var rows = File.ReadAllLines(path)
            .Select((line, i) => i == 0 ? line.TrimStart('﻿') : line)
            .Where(line => line.Trim().Length > 0)
            .Select(line => line.TrimEnd('\r').Split('\t'))
            .ToList();
        if (rows.Count == 0)
            throw new Program.UsageException($"The {what} file '{path}' is empty.");
        return rows;
    }

    private static void RequireUnique(List<string> names, string what, string where)
    {
        if (names.Any(n => n.Length == 0))
            throw new Program.UsageException($"An empty {what} name in the {where}.");
        string? duplicate = names.GroupBy(n => n, StringComparer.Ordinal).FirstOrDefault(g => g.Count() > 1)?.Key;
        if (duplicate != null)
            throw new Program.UsageException($"The {what} '{duplicate}' appears more than once in the {where}.");
    }

    private static int ParseInt(string raw, string what) =>
        int.TryParse(raw.Trim(), NumberStyles.Integer, CultureInfo.InvariantCulture, out int value)
            ? value
            : throw new Program.UsageException($"{what} must be an integer; got '{raw.Trim()}'.");

    /// <summary>A number, or NaN for blank, NA or NaN (any case); anything else is refused.</summary>
    private static double ParseOptional(string raw, string what)
    {
        string cell = raw.Trim();
        if (cell.Length == 0 || cell.Equals("NA", StringComparison.OrdinalIgnoreCase) || cell.Equals("NaN", StringComparison.OrdinalIgnoreCase))
            return double.NaN;
        return double.TryParse(cell, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
            ? value
            : throw new Program.UsageException($"{what}: '{cell}' is not a number. Leave it blank, or write NA, for a missing value.");
    }

    private static double ParseRequired(string raw, string what)
    {
        double value = ParseOptional(raw, what);
        return double.IsFinite(value) ? value : throw new Program.UsageException($"{what} must be a finite number; got '{raw.Trim()}'.");
    }

    private static double? NullIfMissing(double value) => double.IsNaN(value) ? null : value;

    private static char? ParseResidue(string raw, string what, bool optional)
    {
        string cell = raw.Trim();
        if (cell.Length == 0 && optional)
            return null;
        return cell.Length == 1 ? cell[0] : throw new Program.UsageException($"{what} must be one residue letter; got '{cell}'.");
    }

    private static OccupancyState ParseState(string raw, string what)
    {
        string cell = raw.Trim();
        foreach (OccupancyState state in Enum.GetValues<OccupancyState>())
            if (StateName(state) == cell)
                return state;
        throw new Program.UsageException(
            $"{what} must be one of {string.Join(", ", Enum.GetValues<OccupancyState>().Select(StateName))}; got '{cell}'.");
    }

    private static bool ParseBool(string raw, string what) =>
        raw.Trim().ToLowerInvariant() switch
        {
            "true" => true,
            "false" => false,
            _ => throw new Program.UsageException($"{what} must be true or false; got '{raw.Trim()}'."),
        };

    // ---- options ---------------------------------------------------------------------------------

    private static string RequireFile(Program.Arguments arguments, string name)
    {
        RequireValueIfProvided(arguments, name);
        string path = arguments.Required(name);
        if (!File.Exists(path))
            throw new Program.UsageException($"The --{name} file was not found: '{path}'.");
        return path;
    }

    private static string? OptionalFile(Program.Arguments arguments, string name) =>
        arguments.WasProvided(name) ? RequireFile(arguments, name) : null;

    private static void RequireValueIfProvided(Program.Arguments arguments, string name)
    {
        if (arguments.Flag(name))
            throw new Program.UsageException($"Option --{name} was given without a value.");
    }

    /// <summary>A flag, refused when written with a value: <c>--keep-ceiling false</c> would otherwise be ignored.</summary>
    private static bool Flag(Program.Arguments arguments, string name)
    {
        if (arguments.Optional(name) is { } stray)
            throw new Program.UsageException($"Option --{name} takes no value (got '{stray}'); give it alone, or leave it out.");
        return arguments.Flag(name);
    }

    private static int? OptionalCount(Program.Arguments arguments, string name)
    {
        if (!arguments.WasProvided(name))
            return null;
        int value = arguments.OptionalInt(name, 0);
        return value >= 0 ? value : throw new Program.UsageException($"Option --{name} must be zero or more; got {value}.");
    }

    /// <summary>The <c>--out</c> path, checked before anything is read: a <c>.tsv</c> (PYB-1), or null.</summary>
    private static string? OutputPath(Program.Arguments arguments)
    {
        RequireValueIfProvided(arguments, "out");
        string? path = arguments.Optional("out");
        if (path is not null && string.IsNullOrWhiteSpace(path))
            throw new Program.UsageException("Option --out needs a file path ending in .tsv.");
        Program.RequireTsvOutput(path);
        return path;
    }

    private static bool SamePath(string a, string b) =>
        string.Equals(Path.GetFullPath(a), Path.GetFullPath(b), StringComparison.OrdinalIgnoreCase);

    /// <summary>
    /// An mzLib call whose <see cref="ArgumentException"/> is the caller's input refused (a repeated row,
    /// coordinates that do not fit the sequence, a run with no sample, a trait or covariate that is not
    /// finite, a design that is not of full rank): reclassified as a usage error, as <c>stats fit</c> does.
    /// </summary>
    private static T FromMzLib<T>(Func<T> call)
    {
        try
        {
            return call();
        }
        catch (ArgumentException invalid) when (invalid is not ArgumentNullException)
        {
            throw new Program.UsageException(invalid.Message.Split(" (Parameter")[0]);
        }
    }

    // ---- output ----------------------------------------------------------------------------------

    /// <summary>
    /// The table as plain tab-separated text: one header row, then one row per cell, null as an empty cell,
    /// numbers in invariant culture at round-trip precision, booleans as true/false. No value can hold a tab
    /// or a line break (each came out of a tab-separated input), so nothing is quoted, and
    /// <see cref="ReadOccupancy"/> reads it back exactly.
    /// </summary>
    private static object WriteTsv(string outputPath, Sdrf.ColumnTable table)
    {
        string full = Path.GetFullPath(outputPath);
        string? directory = Path.GetDirectoryName(full);
        if (!string.IsNullOrEmpty(directory))
            Directory.CreateDirectory(directory);
        Dictionary<string, List<object?>> columns = table.Columns();
        var text = new StringBuilder();
        text.Append(string.Join('\t', table.Names)).Append('\n');
        for (int r = 0; r < table.RowCount; r++)
            text.Append(string.Join('\t', table.Names.Select(n => Cell(columns[n][r])))).Append('\n');
        File.WriteAllText(full, text.ToString(), new UTF8Encoding(false));
        return new { path = full, format = "tsv", row_count = table.RowCount };
    }

    private static string Cell(object? value) => value switch
    {
        null => "",
        double number => number.ToString("R", CultureInfo.InvariantCulture),
        int number => number.ToString(CultureInfo.InvariantCulture),
        bool flag => flag ? "true" : "false",
        _ => value.ToString() ?? "",
    };

    private static double? Finite(double value) => double.IsFinite(value) ? value : null;
}
