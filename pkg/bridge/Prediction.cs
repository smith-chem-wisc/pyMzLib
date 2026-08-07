using System.Globalization;
using System.Reflection;
using Chromatography.RetentionTimePrediction;
using PredictionClients.Koina.AbstractClasses;

namespace MzLibBridge;

/// <summary>
/// Peptide-property prediction: retention time, fragment intensities, collisional cross-section,
/// and detectability, from mzLib's Koina clients and its local analytic predictors.
/// </summary>
/// <remarks>
/// <para>
/// mzLib ships clients for <b>40 published models</b> on the
/// <a href="https://koina.wilhelmlab.org/">Koina</a> inference server, across five families, plus
/// three predictors that run locally. This exposes the Koina families and the two local predictors
/// that are pure arithmetic.
/// </para>
/// <para>
/// <b>The local Chronologer predictor is deliberately not exposed.</b> It is a TorchSharp neural
/// network: x64-only, extracting hundreds of megabytes of weights to a shared temp path at
/// construction, and racing with any concurrent process that does the same. pyMzLib publishes an
/// arm64 macOS wheel, so exposing it would either break that wheel or ship a verb that fails on it.
/// The same model is reachable over Koina as <c>Chronologer_RT</c>, which needs no native code —
/// with the caveat that the two report <i>different units</i>: the network model returns absolute
/// retention time and the local one returns % acetonitrile. That divergence is precisely why this
/// verb reports a unit as a value rather than leaving it to prose.
/// </para>
/// <para>
/// <b>Every verb here is stateless.</b> mzLib's <c>IPredictor</c> is not: <c>Predict</c> overwrites
/// <c>ModelInputs</c>, <c>ValidInputsMask</c> and <c>Predictions</c> on the instance, and the
/// library documents its models as not thread-safe. A wire contract that exposed that would be
/// unusable from three languages, so a model is constructed, used once, and discarded within a
/// single call, and everything the caller could want comes back from that call.
/// </para>
/// <para>
/// <b>Inputs arrive on stdin as a tab-separated table</b>, not on argv. A real prediction run is
/// thousands of peptides and argv has a hard ceiling of roughly 32 KB — the same reasoning that
/// put PRIDE's explicit file selection on stdin. A header line names the columns, so the same
/// transport serves families that need one field and families that need five. <c>--sequence</c>
/// remains for trying a single peptide without composing a table.
/// </para>
/// </remarks>
internal static class Prediction
{
    /// <summary>
    /// <c>predict models [--family F]</c> — every Koina model mzLib can call, with its constraints.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Enumerated by reflection over the loaded assembly rather than transcribed, for the same
    /// reason <c>readers formats</c> enumerates <c>SupportedFileType</c>: a published catalogue that
    /// is maintained by hand drifts from what the library can actually call.
    /// </para>
    /// <para>
    /// <b>mzLib has no runtime registry of its own.</b> The equivalent reflection lives in mzLib's
    /// <i>test</i> project (<c>KoinaModelDiscoveryTests.ConcreteKoinaModelTypes</c>), so a native C#
    /// consumer — a MetaMorpheus dialog offering a model list, say — has exactly the gap this fills
    /// and cannot fill it without copying test code. That makes promoting this into mzLib proper a
    /// mainland concern rather than a binding one; it is logged on <c>bridge/UPSTREAM.md</c>, and
    /// what lives here is the waystation until it lands.
    /// </para>
    /// </remarks>
    public static object Models(Program.Arguments arguments)
    {
        string? family = arguments.Optional("family");
        if (arguments.WasProvided("family") && string.IsNullOrWhiteSpace(family))
            throw new Program.UsageException("Option --family was given but is empty; omit it to list every family.");

        var models = new List<object>();
        foreach (Type type in ConcreteModelTypes())
        {
            object? model = TryConstruct(type);
            if (model is null)
            {
                // A model whose constructor changed shape is a broken mzLib build, not a reason to
                // take the whole catalogue down — the catalogue is the verb you would use to
                // diagnose exactly that.
                models.Add(new { type = type.Name, family = FamilyOf(type), error = "could not be constructed with its default arguments" });
                continue;
            }

            string modelFamily = FamilyOf(type);
            if (family is not null && !string.Equals(modelFamily, family, StringComparison.OrdinalIgnoreCase))
                continue;

            models.Add(Describe(model, type, modelFamily));
        }

        return new
        {
            model_count = models.Count,
            families = FamilyVerbs,
            // Named so a caller can see who they are calling, and throttle politely.
            service = PredictionClients.Koina.Client.HTTP.ModelsURL,
            models,
        };
    }

    /// <summary>
    /// <c>predict retention-time --model M</c> — predicted elution, one row per peptide.
    /// </summary>
    /// <remarks>
    /// <b>The unit is not the same for every model, and is reported as a value.</b> Most Koina RT
    /// models are trained on <i>indexed</i> retention time — iRT, a dimensionless scale anchored to
    /// standard peptides — while <c>Chronologer_RT</c> returns absolute minutes. mzLib carries the
    /// distinction as a bare <c>IsIndexed</c> boolean on each prediction; a boolean is not a unit,
    /// so it crosses here as <c>retention_time_unit</c>, alongside the same <c>"minutes"</c> /
    /// <c>"unknown"</c> vocabulary the readers verbs use. An iRT value plotted against a gradient
    /// in minutes is the mistake this exists to prevent.
    /// </remarks>
    public static object RetentionTime(Program.Arguments arguments)
    {
        var model = Build<RetentionTimeModel>(arguments);
        List<InputRow> rows = ReadInput(arguments, "sequence");

        List<PeptideRTPrediction> predictions = model.Predict(
            rows.Select(row => new RetentionTimePredictionInput(row.Required("sequence"))).ToList());

        return Table(
            arguments, model.ModelName, predictions.Count,
            RetentionTimeColumns,
            predictions,
            extra: new
            {
                // "indexed_retention_time" rather than a bool: iRT is dimensionless and anchored to
                // standard peptides, so it is a different quantity from minutes, not a scaling of it.
                retention_time_unit = predictions.Count > 0 && predictions[0].IsIndexed == true
                    ? "indexed_retention_time"
                    : predictions.Count > 0 && predictions[0].IsIndexed == false
                        ? "minutes"
                        : "unknown",
                caveats = RetentionTimeCaveats(model),
            });
    }

    /// <summary>
    /// <c>predict fragments --model M</c> — predicted MS2 fragment m/z and relative intensity.
    /// </summary>
    /// <remarks>
    /// <b>The arrays are ragged, and that is mzLib's answer rather than a truncation.</b> Koina
    /// returns a fixed-width grid with <c>-1</c> marking ions that cannot exist for a given peptide
    /// and charge, and mzLib drops those before handing the result back. So a peptide's three
    /// arrays are as long as its <i>possible</i> ions, not as long as the model's nominal ion
    /// count, and two peptides in one call will differ. Each row therefore carries its own arrays.
    /// </remarks>
    public static object Fragments(Program.Arguments arguments)
    {
        var model = Build<FragmentIntensityModel>(arguments);
        List<InputRow> rows = ReadInput(arguments, "sequence", "precursor_charge");

        List<PeptideFragmentIntensityPrediction> predictions = model.Predict(
            rows.Select(row => new FragmentIntensityPredictionInput(
                row.Required("sequence"),
                row.RequiredInt("precursor_charge"),
                row.OptionalInt("collision_energy"),
                row.Optional("instrument_type"),
                row.Optional("fragmentation_type"))).ToList());

        return Table(
            arguments, model.ModelName, predictions.Count,
            FragmentColumns,
            predictions,
            extra: new
            {
                intensity_scale = "relative",
                caveats = new[]
                {
                    "fragment_intensity is RELATIVE, on Koina's own 0-1 scale, and is not comparable " +
                    "with a measured intensity or between models.",
                    "The three fragment arrays are RAGGED: Koina returns a fixed-width grid with -1 " +
                    "marking ions that cannot exist for a peptide, and mzLib drops those, so each " +
                    "row's arrays are as long as its possible ions rather than the model's nominal " +
                    "ion count. Index them per row, never as a rectangle.",
                    "Predictions are a model's opinion, not a measurement. Nothing here has been " +
                    "matched against a spectrum.",
                },
            });
    }

    /// <summary>
    /// <c>predict ccs --model M</c> — predicted collisional cross-section, in square angstroms.
    /// </summary>
    public static object Ccs(Program.Arguments arguments)
    {
        var model = Build<CollisionalCrossSectionModel>(arguments);
        List<InputRow> rows = ReadInput(arguments, "sequence", "precursor_charge");

        List<PeptideCCSPrediction> predictions = model.Predict(
            rows.Select(row => new CCSPredictionInput(
                row.Required("sequence"), row.RequiredInt("precursor_charge"))).ToList());

        return Table(
            arguments, model.ModelName, predictions.Count,
            CcsColumns,
            predictions,
            extra: new
            {
                // Stated as a value because the alternative unit a reader will assume is 1/K0, the
                // reduced mobility a timsTOF actually reports, and the conversion between them
                // needs the drift-gas parameters mzLib does not carry.
                collisional_cross_section_unit = "square_angstroms",
                caveats = new[]
                {
                    "collisional_cross_section is in SQUARE ANGSTROMS, not 1/K0. Converting to the " +
                    "reduced mobility a timsTOF reports needs drift-gas temperature and pressure, " +
                    "which mzLib does not carry, so no conversion is offered here.",
                },
            });
    }

    /// <summary>
    /// <c>predict detectability --model M</c> — predicted flyability, as four class probabilities.
    /// </summary>
    public static object Detectability(Program.Arguments arguments)
    {
        var model = Build<DetectabilityModel>(arguments);
        List<InputRow> rows = ReadInput(arguments, "sequence");

        List<PeptideDetectabilityPrediction> predictions = model.Predict(
            rows.Select(row => new DetectabilityPredictionInput(row.Required("sequence"))).ToList());

        return Table(
            arguments, model.ModelName, predictions.Count,
            DetectabilityColumns,
            predictions,
            extra: new
            {
                caveats = new[]
                {
                    "The four probabilities are a distribution over classes and sum to 1 for each " +
                    "peptide. They are not an expected intensity and not a detection probability.",
                },
            });
    }

    /// <summary>
    /// <c>predict crosslink-fragments --model M</c> — MS2 intensities for a crosslinked pair.
    /// </summary>
    /// <remarks>
    /// <b>This family takes a different sequence language from every other verb here.</b> The other
    /// four accept mzLib's <c>FullSequence</c> notation and convert it; the crosslink models reject
    /// it and require raw UNIMOD brackets — <c>K[UNIMOD:1896]</c>. That is mzLib's constraint, not
    /// a choice made here, and it is stated in the response's caveats because a same-named
    /// <c>sequence</c> column meaning two different things is exactly the trap a wire contract
    /// should not hide.
    /// </remarks>
    public static object CrosslinkFragments(Program.Arguments arguments)
    {
        var model = Build<CrosslinkFragmentIntensityModel>(arguments);
        List<InputRow> rows = ReadInput(arguments, "alpha_sequence", "precursor_charge");

        List<CrosslinkFragmentIntensityPrediction> predictions = model.Predict(
            rows.Select(row => new CrosslinkIntensityPredictionInput(
                row.Required("alpha_sequence"),
                row.Optional("beta_sequence"),
                row.RequiredInt("precursor_charge"),
                row.OptionalInt("collision_energy"))).ToList());

        return Table(
            arguments, model.ModelName, predictions.Count,
            CrosslinkColumns,
            predictions,
            extra: new
            {
                intensity_scale = "relative",
                caveats = new[]
                {
                    "alpha_sequence and beta_sequence must use RAW UNIMOD notation for this family " +
                    "- 'K[UNIMOD:1896]' - not mzLib's FullSequence notation, which the crosslink " +
                    "models reject. Every other predict verb accepts mzLib notation and converts it. " +
                    "Same column name, two different input languages.",
                    "The three fragment arrays are RAGGED, as for predict fragments.",
                },
            });
    }

    // ---------------------------------------------------------------------------------------
    // Column sets
    // ---------------------------------------------------------------------------------------

    /// <summary>The retentiontime view's columns, named so a test can exercise the
    /// same projection the verb uses rather than a copy of it.</summary>
    internal static readonly (string Name, Func<PeptideRTPrediction, object?> Read)[] RetentionTimeColumns =
    [
                ("sequence", p => p.FullSequence),
                ("validated_sequence", p => p.ValidatedFullSequence),
                ("retention_time", p => p.PredictedRetentionTime),
                ("warning", p => p.Warning?.Message),
    ];

    /// <summary>The fragment view's columns, named so a test can exercise the
    /// same projection the verb uses rather than a copy of it.</summary>
    internal static readonly (string Name, Func<PeptideFragmentIntensityPrediction, object?> Read)[] FragmentColumns =
    [
                ("sequence", p => p.FullSequence),
                ("validated_sequence", p => p.ValidatedFullSequence),
                ("precursor_charge", p => p.PrecursorCharge),
                ("fragment_annotations", p => p.FragmentAnnotations),
                ("fragment_mz", p => p.FragmentMZs),
                ("fragment_intensity", p => p.FragmentIntensities),
                ("warning", p => p.Warning?.Message),
    ];

    /// <summary>The ccs view's columns, named so a test can exercise the
    /// same projection the verb uses rather than a copy of it.</summary>
    internal static readonly (string Name, Func<PeptideCCSPrediction, object?> Read)[] CcsColumns =
    [
                ("sequence", p => p.FullSequence),
                ("validated_sequence", p => p.ValidatedFullSequence),
                ("precursor_charge", p => p.PrecursorCharge),
                ("collisional_cross_section", p => p.PredictedCCS),
                ("warning", p => p.Warning?.Message),
    ];

    /// <summary>The detectability view's columns, named so a test can exercise the
    /// same projection the verb uses rather than a copy of it.</summary>
    internal static readonly (string Name, Func<PeptideDetectabilityPrediction, object?> Read)[] DetectabilityColumns =
    [
                ("sequence", p => p.FullSequence),
                ("validated_sequence", p => p.ValidatedFullSequence),
                ("not_detectable", p => p.DetectabilityProbabilities?.NotDetectable),
                ("low_detectability", p => p.DetectabilityProbabilities?.LowDetectability),
                ("intermediate_detectability", p => p.DetectabilityProbabilities?.IntermediateDetectability),
                ("high_detectability", p => p.DetectabilityProbabilities?.HighDetectability),
                ("warning", p => p.Warning?.Message),
    ];

    /// <summary>The crosslink view's columns, named so a test can exercise the
    /// same projection the verb uses rather than a copy of it.</summary>
    internal static readonly (string Name, Func<CrosslinkFragmentIntensityPrediction, object?> Read)[] CrosslinkColumns =
    [
                ("alpha_sequence", p => p.AlphaSequence),
                ("beta_sequence", p => p.BetaSequence),
                ("validated_alpha_sequence", p => p.ValidatedAlphaSequence),
                ("validated_beta_sequence", p => p.ValidatedBetaSequence),
                ("precursor_charge", p => p.PrecursorCharge),
                ("fragment_annotations", p => p.FragmentAnnotations),
                ("fragment_mz", p => p.FragmentMZs),
                ("fragment_intensity", p => p.FragmentIntensities),
                ("warning", p => p.Warning?.Message),
    ];

    // ---------------------------------------------------------------------------------------
    // The model catalogue
    // ---------------------------------------------------------------------------------------

    /// <summary>The verb each model family is called through.</summary>
    private static readonly Dictionary<string, string> FamilyVerbs = new()
    {
        ["retention_time"] = "predict retention-time",
        ["fragment_intensity"] = "predict fragments",
        ["collisional_cross_section"] = "predict ccs",
        ["detectability"] = "predict detectability",
        ["crosslink_intensity"] = "predict crosslink-fragments",
    };

    private static IEnumerable<Type> ConcreteModelTypes() =>
        typeof(FragmentIntensityModel).Assembly
            .GetTypes()
            .Where(type => type is { IsAbstract: false, IsPublic: true } && FamilyOf(type) != "unknown")
            .OrderBy(type => FamilyOf(type), StringComparer.Ordinal)
            .ThenBy(type => type.Name, StringComparer.Ordinal);

    private static string FamilyOf(Type type)
    {
        for (Type? current = type; current is not null; current = current.BaseType)
        {
            if (current == typeof(RetentionTimeModel)) return "retention_time";
            // Crosslink is checked before fragment intensity because it is NOT a subclass of it —
            // they are siblings — but the names invite the assumption, and getting it wrong would
            // route crosslink models to a verb that cannot take their beta sequence.
            if (current == typeof(CrosslinkFragmentIntensityModel)) return "crosslink_intensity";
            if (current == typeof(FragmentIntensityModel)) return "fragment_intensity";
            if (current == typeof(CollisionalCrossSectionModel)) return "collisional_cross_section";
            if (current == typeof(DetectabilityModel)) return "detectability";
        }

        return "unknown";
    }

    /// <summary>
    /// Builds a model with the arguments its own constructor declares as defaults.
    /// </summary>
    /// <remarks>
    /// <b>Not <see cref="Activator.CreateInstance(Type)"/>.</b> Every Koina model has a constructor
    /// whose parameters are all optional — <c>Prosit2019iRT(SequenceConversionHandlingMode = …, int
    /// = 500, int = 100)</c> — which reads as "constructible with no arguments" in C# and is not
    /// one to reflection: a default-valued parameter is still a parameter, so the parameterless
    /// overload the activator looks for does not exist and <i>every</i> model failed to construct.
    /// The declared defaults are read off the constructor and passed explicitly, so the models are
    /// configured exactly as a C# caller writing <c>new Prosit2019iRT()</c> would get them, rather
    /// than as this file's guess at sensible values.
    /// </remarks>
    private static object? TryConstruct(Type type)
    {
        try
        {
            ConstructorInfo? constructor = type
                .GetConstructors()
                .Where(candidate => candidate.GetParameters().All(parameter => parameter.HasDefaultValue))
                // Fewest parameters first, so a genuinely parameterless constructor wins over one
                // that merely defaults everything.
                .OrderBy(candidate => candidate.GetParameters().Length)
                .FirstOrDefault();

            return constructor?.Invoke(
                constructor.GetParameters().Select(parameter => parameter.DefaultValue).ToArray());
        }
        catch (Exception)
        {
            return null;
        }
    }

    /// <summary>One model's published constraints, read off the instance rather than transcribed.</summary>
    private static object Describe(object model, Type type, string family)
    {
        Func<string, object?> read = name => type.GetProperty(name)?.GetValue(model);

        return new
        {
            model = read("ModelName"),
            family,
            verb = FamilyVerbs.TryGetValue(family, out string? verb) ? verb : null,
            type = type.Name,
            min_peptide_length = read("MinPeptideLength"),
            max_peptide_length = read("MaxPeptideLength"),
            max_batch_size = read("MaxBatchSize"),
            // Tri-state, not a list — see Constraint. Emitting the raw collection would read
            // exactly backwards.
            precursor_charge = Constraint(read("AllowedPrecursorCharges")),
            collision_energy = Constraint(read("AllowedCollisionEnergies")),
            instrument_type = Constraint(read("AllowedInstrumentTypes")),
            fragmentation_type = Constraint(read("AllowedFragmentationTypes")),
            // Empty means the model accepts no modifications at all, which is a real answer.
            allowed_unimod_ids = Sorted(read("AllowedUnimodIds")),
            // Only meaningful for the retention-time family; null elsewhere.
            retention_time_unit = family == "retention_time"
                ? read("IsIndexedRetentionTimeModel") is true ? "indexed_retention_time" : "minutes"
                : null,
            number_of_predicted_fragment_ions = NullIfNegative(read("NumberOfPredictedFragmentIons")),
        };
    }

    /// <summary>
    /// What a model requires of an optional input parameter, as a value rather than as a
    /// collection whose emptiness has to be interpreted.
    /// </summary>
    /// <remarks>
    /// <para>
    /// mzLib encodes three states in one nullable set, and documents them in a comment on
    /// <c>FragmentIntensityModel.AllowedInstrumentTypes</c>:
    /// </para>
    /// <list type="bullet">
    /// <item><c>null</c> — the parameter is <b>not applicable</b>; validation is skipped.</item>
    /// <item><b>empty</b> — the parameter <b>is required</b>, and any value is accepted.</item>
    /// <item><b>populated</b> — required, and only these values.</item>
    /// </list>
    /// <para>
    /// Passing the raw collection through would be faithful and useless: a caller reads <c>null</c>
    /// as "no constraint, send anything" and <c>[]</c> as "nothing is allowed", which is the
    /// opposite of both meanings. <c>Prosit_2020_intensity_HCD</c> returns an empty set for
    /// collision energy and <i>requires</i> one; <c>Prosit_2020_intensity_CID</c> returns null and
    /// rejects the idea, being fixed at NCE 35. Encoding the requirement as a word is the whole
    /// point of a wire contract that claims to carry availability rather than raw fields.
    /// </para>
    /// </remarks>
    private static object Constraint(object? value) => value switch
    {
        null => new { requirement = "not_applicable", values = (object?)null },
        IEnumerable<int> numbers when !numbers.Any() => new { requirement = "any_value_required", values = (object?)null },
        IEnumerable<int> numbers => new { requirement = "one_of", values = (object?)numbers.Order().ToList() },
        IEnumerable<string> text when !text.Any() => new { requirement = "any_value_required", values = (object?)null },
        IEnumerable<string> text => new { requirement = "one_of", values = (object?)text.Order(StringComparer.Ordinal).ToList() },
        _ => new { requirement = "not_applicable", values = (object?)null },
    };

    private static List<int>? Sorted(object? value) =>
        value is IEnumerable<int> numbers ? numbers.Order().ToList() : null;

    /// <summary>mzLib's -1 "the count is dynamic" marker, as null.</summary>
    private static int? NullIfNegative(object? value) =>
        value is int number && number >= 0 ? number : null;

    // ---------------------------------------------------------------------------------------
    // Model construction
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// Builds the named model of a family, or explains which names that family has.
    /// </summary>
    /// <remarks>
    /// Matching is on the model's published Koina name, not its .NET class name, because that is
    /// the identifier a user finds in the Koina catalogue and in the literature.
    /// </remarks>
    internal static TModel Build<TModel>(Program.Arguments arguments) where TModel : class
    {
        string requested = arguments.Required("model");

        List<(Type Type, object Instance)> candidates = typeof(TModel).Assembly
            .GetTypes()
            .Where(type => type is { IsAbstract: false } && typeof(TModel).IsAssignableFrom(type))
            .Select(type => (Type: type, Instance: TryConstruct(type)))
            .Where(pair => pair.Instance is not null)
            .Select(pair => (pair.Type, Instance: pair.Instance!))
            .ToList();

        foreach ((Type type, object instance) in candidates)
        {
            if (string.Equals(NameOf(type, instance), requested, StringComparison.OrdinalIgnoreCase))
                return ApplyPoliteness((TModel)instance, arguments);
        }

        List<string> names = candidates
            .Select(pair => NameOf(pair.Type, pair.Instance))
            .Where(name => name is not null)
            .Select(name => name!)
            .Order(StringComparer.Ordinal)
            .ToList();

        throw new Program.UsageException(
            $"No model named '{requested}' for this verb. Available: {string.Join(", ", names)}. " +
            "The models listing describes every model with its constraints.");
    }

    private static string? NameOf(Type type, object instance) =>
        type.GetProperty("ModelName")?.GetValue(instance) as string;

    /// <summary>
    /// Applies the caller's throttling choices, if any.
    /// </summary>
    /// <remarks>
    /// <b>The defaults are mzLib's and are deliberately left alone.</b> Koina is a public,
    /// shared, community-run GPU server, not a service anyone here pays for. The knobs are exposed
    /// because a caller with a legitimate large job needs them, but a binding that quietly
    /// maximised throughput would be spending someone else's capacity by default.
    /// </remarks>
    private static TModel ApplyPoliteness<TModel>(TModel model, Program.Arguments arguments) where TModel : class
    {
        Set("max-batches", "MaxNumberOfBatchesPerRequest", minimum: 1);
        Set("throttle-ms", "ThrottlingDelayInMilliseconds", minimum: 0);
        return model;

        void Set(string option, string property, int minimum)
        {
            if (!arguments.WasProvided(option))
                return;
            if (string.IsNullOrWhiteSpace(arguments.Optional(option)))
                throw new Program.UsageException(
                    $"Option --{option} was given but has no value; omit it to use mzLib's default.");

            int value = arguments.OptionalInt(option, minimum);
            if (value < minimum)
                throw new Program.UsageException($"Option --{option} must be {minimum} or greater; got {value}.");

            // init-only on the model, so it is set through the backing property rather than the
            // initialiser; a model that stops exposing it should fail loudly rather than silently
            // ignore the caller's throttling choice.
            PropertyInfo info = model.GetType().GetProperty(property)
                ?? throw new Program.UsageException(
                    $"This model does not expose {property}, so --{option} cannot be honoured.");
            info.SetValue(model, value);
        }
    }

    // ---------------------------------------------------------------------------------------
    // Input
    // ---------------------------------------------------------------------------------------

    /// <summary>One row of the input table, by column name.</summary>
    private sealed class InputRow(Dictionary<string, string> values, int lineNumber)
    {
        public string Required(string column) =>
            values.TryGetValue(column, out string? value) && !string.IsNullOrWhiteSpace(value)
                ? value
                : throw new Program.UsageException(
                    $"Input line {lineNumber} has no value for the required column '{column}'.");

        public string? Optional(string column) =>
            values.TryGetValue(column, out string? value) && !string.IsNullOrWhiteSpace(value) ? value : null;

        public int RequiredInt(string column) => ParseInt(column, Required(column));

        public int? OptionalInt(string column)
        {
            string? raw = Optional(column);
            return raw is null ? null : ParseInt(column, raw);
        }

        private int ParseInt(string column, string raw) =>
            int.TryParse(raw, NumberStyles.Integer, CultureInfo.InvariantCulture, out int parsed)
                ? parsed
                : throw new Program.UsageException(
                    $"Input line {lineNumber}: column '{column}' must be a whole number; got '{raw}'.");
    }

    /// <summary>
    /// The input table, from <c>--sequence</c> or from a tab-separated table on stdin.
    /// </summary>
    /// <remarks>
    /// A header line is required when reading stdin, and an unknown column is an error rather than
    /// being ignored: a caller who wrote <c>charge</c> where the contract says
    /// <c>precursor_charge</c> has asked for something specific, and silently predicting at a
    /// default charge would give them a plausible wrong answer.
    /// </remarks>
    private static List<InputRow> ReadInput(Program.Arguments arguments, params string[] required)
    {
        string? single = arguments.Optional("sequence");
        if (arguments.WasProvided("sequence"))
        {
            if (string.IsNullOrWhiteSpace(single))
                throw new Program.UsageException("Option --sequence was given but is empty.");

            var values = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
            {
                [required[0]] = single,
            };
            foreach (string column in KnownColumns)
            {
                string? value = arguments.Optional(column.Replace('_', '-'));
                if (!string.IsNullOrWhiteSpace(value))
                    values[column] = value;
            }

            return [new InputRow(values, 1)];
        }

        var rows = new List<InputRow>();
        string? header = Console.In.ReadLine();
        if (header is null)
            throw new Program.UsageException(
                $"No input. Provide --sequence for one peptide, or a tab-separated table on stdin " +
                $"whose header names at least: {string.Join(", ", required)}.");

        string[] columns = header.Split('\t').Select(name => name.Trim()).ToArray();
        List<string> unknown = columns.Where(name => !KnownColumns.Contains(name, StringComparer.OrdinalIgnoreCase)).ToList();
        if (unknown.Count > 0)
            throw new Program.UsageException(
                $"Unknown input column(s): {string.Join(", ", unknown)}. Known columns: " +
                $"{string.Join(", ", KnownColumns)}.");

        foreach (string missing in required.Where(name => !columns.Contains(name, StringComparer.OrdinalIgnoreCase)))
            throw new Program.UsageException($"The input table has no '{missing}' column, which this verb requires.");

        int lineNumber = 1;
        string? line;
        while ((line = Console.In.ReadLine()) != null)
        {
            lineNumber++;
            if (string.IsNullOrWhiteSpace(line))
                continue;

            string[] cells = line.Split('\t');
            var values = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            for (int i = 0; i < columns.Length && i < cells.Length; i++)
                values[columns[i]] = cells[i].Trim();

            rows.Add(new InputRow(values, lineNumber));
        }

        if (rows.Count == 0)
            throw new Program.UsageException("The input table has a header but no rows.");

        return rows;
    }

    /// <summary>Every column any predict verb understands.</summary>
    private static readonly string[] KnownColumns =
    [
        "sequence", "alpha_sequence", "beta_sequence", "precursor_charge",
        "collision_energy", "instrument_type", "fragmentation_type",
    ];

    // ---------------------------------------------------------------------------------------
    // Output
    // ---------------------------------------------------------------------------------------

    /// <summary>
    /// The columnar payload every predict verb returns, in the shape the readers verbs use.
    /// </summary>
    /// <remarks>
    /// Deliberately the same shape as <c>readers read-records</c> — one array per column, names in
    /// order — so a caller who has learned one table shape has learned all of them, and a binding
    /// can reuse its table type rather than writing a second one.
    /// </remarks>
    internal static object Table<T>(
        Program.Arguments arguments,
        string modelName,
        int rowCount,
        (string Name, Func<T, object?> Read)[] columns,
        List<T> rows,
        object extra)
    {
        string? outputPath = arguments.Optional("out");
        if (arguments.WasProvided("out") && string.IsNullOrWhiteSpace(outputPath))
            throw new Program.UsageException("Option --out was given but has no value; omit it to return the table.");

        var built = new Dictionary<string, List<object?>>(columns.Length);
        foreach ((string name, Func<T, object?> read) in columns)
            built[name] = rows.Select(read).ToList();

        object? written = null;
        if (!string.IsNullOrWhiteSpace(outputPath))
            written = WriteTable(outputPath, columns, rows);

        // Rows whose prediction failed are still rows, with nulls and a warning — a caller must be
        // able to line predictions up against the peptides they sent.
        int failed = rows.Count(row => columns.Any(c => c.Name == "warning" && c.Read(row) is not null));

        return Merge(
            new
            {
                model = modelName,
                row_count = rowCount,
                // Not an error: an unpredictable peptide is a normal outcome (too long, an
                // unsupported modification), and the row survives so the alignment does.
                failed_row_count = failed,
                column_names = columns.Select(c => c.Name).ToList(),
                columns = written is null ? built : null,
                output = written,
            },
            extra);
    }

    internal static object WriteTable<T>(
        string outputPath, (string Name, Func<T, object?> Read)[] columns, List<T> rows)
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
            foreach ((string name, _) in columns)
                csv.WriteField(name);
            csv.NextRecord();

            foreach (T row in rows)
            {
                foreach ((_, Func<T, object?> read) in columns)
                    csv.WriteField(Render(read(row)));
                csv.NextRecord();
            }
        }

        return new
        {
            path = Path.GetFullPath(outputPath),
            format = "tsv",
            row_count = rows.Count,
        };
    }

    internal static string Render(object? value) => value switch
    {
        null => string.Empty,
        double number => number.ToString(CultureInfo.InvariantCulture),
        int number => number.ToString(CultureInfo.InvariantCulture),
        bool flag => flag ? "true" : "false",
        string text => text,
        System.Collections.IEnumerable sequence =>
            string.Join(";", sequence.Cast<object?>().Select(Render)),
        _ => value.ToString() ?? string.Empty,
    };

    /// <summary>Flattens the shared fields and a verb's own into one wire object.</summary>
    internal static Dictionary<string, object?> Merge(object shared, object extra)
    {
        var merged = new Dictionary<string, object?>();
        foreach (object source in new[] { shared, extra })
            foreach (PropertyInfo property in source.GetType().GetProperties())
                merged[SnakeCase(property.Name)] = property.GetValue(source);

        return merged;
    }

    private static string SnakeCase(string name) => name;

    internal static string[] RetentionTimeCaveats(RetentionTimeModel model)
    {
        var caveats = new List<string>
        {
            "A predicted retention time is a model's opinion, not a measurement, and is only " +
            "comparable with an observed gradient after alignment against shared peptides.",
        };

        if (model.IsIndexedRetentionTimeModel)
        {
            caveats.Add(
                "retention_time is an INDEXED retention time (iRT): a dimensionless scale anchored " +
                "to standard peptides, not minutes. Plotting it against a gradient without first " +
                "fitting the iRT-to-minutes line is the commonest way to misread this number.");
        }

        if (string.Equals(model.ModelName, "Chronologer_RT", StringComparison.Ordinal))
        {
            caveats.Add(
                "This is the NETWORK Chronologer, which returns absolute retention time. mzLib also " +
                "ships a local Chronologer that returns % acetonitrile from the same weights - the " +
                "two are not interchangeable, and the local one is not exposed here because it " +
                "needs x64-only native libraries.");
        }

        return caveats.ToArray();
    }
}
