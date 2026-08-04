using System.Diagnostics.CodeAnalysis;
using System.Text.Json;
using PredictionClients.Koina.AbstractClasses;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the <c>predict</c> verbs.
/// </summary>
/// <remarks>
/// <para>
/// Offline. Everything here exercises the catalogue and the input parsing, neither of which
/// touches the network: <c>predict models</c> reflects over the loaded assembly, and every input
/// validation happens before a request is composed. The live canaries belong in the Python and
/// Rust suites, where they can be skipped on an outage without failing this job.
/// </para>
/// <para>
/// The catalogue tests are the important ones. mzLib has no runtime model registry — the
/// equivalent reflection lives in its own <i>test</i> project — so this is the only place the
/// published constraints are checked against the classes that declare them.
/// </para>
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PredictionTests
{
    private string _tempDirectory = string.Empty;

    [SetUp]
    public void CreateTempDirectory()
    {
        _tempDirectory = Path.Combine(Path.GetTempPath(), $"pymzlib-prediction-tests-{Guid.NewGuid():N}");
        Directory.CreateDirectory(_tempDirectory);
    }

    [TearDown]
    public void RemoveTempDirectory()
    {
        if (Directory.Exists(_tempDirectory))
            Directory.Delete(_tempDirectory, recursive: true);
    }

    private static JsonElement Models(params string[] extra)
    {
        string[] args = ["predict", "models", .. extra];
        return Invoke(args);
    }

    private static IEnumerable<JsonElement> ModelList(JsonElement catalogue) =>
        catalogue.GetProperty("models").EnumerateArray();

    private static JsonElement ModelNamed(string name) =>
        ModelList(Models()).First(model =>
            model.TryGetProperty("model", out JsonElement value) && value.GetString() == name);

    [Test]
    public void EveryKoinaModelIsConstructibleAndDescribed()
    {
        // Enumerated by reflection rather than transcribed, so a model added to mzLib appears here
        // without anyone editing this file — and one whose constructor changes shape shows up as a
        // described failure rather than as a silently missing entry.
        JsonElement catalogue = Models();

        List<JsonElement> models = ModelList(catalogue).ToList();
        List<string> broken = models
            .Where(model => model.TryGetProperty("error", out _))
            .Select(model => model.GetProperty("type").GetString()!)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(models, Is.Not.Empty);
            Assert.That(broken, Is.Empty, $"models that could not be constructed: {string.Join(", ", broken)}");
            Assert.That(catalogue.GetProperty("model_count").GetInt32(), Is.EqualTo(models.Count));
        });
    }

    [Test]
    public void EveryModelBelongsToAKnownFamilyAndNamesItsVerb()
    {
        string[] families =
        [
            "retention_time", "fragment_intensity", "collisional_cross_section",
            "detectability", "crosslink_intensity",
        ];

        Assert.Multiple(() =>
        {
            foreach (JsonElement model in ModelList(Models()))
            {
                Assert.That(model.GetProperty("family").GetString(), Is.AnyOf(families));
                Assert.That(model.GetProperty("verb").GetString(), Is.Not.Null.And.Not.Empty);
            }
        });
    }

    [Test]
    public void CrosslinkModelsAreNotFiledUnderFragmentIntensity()
    {
        // They are siblings, not subclasses, but the names invite the assumption — and getting it
        // wrong would route them to a verb that cannot take their beta sequence.
        foreach (JsonElement model in ModelList(Models("--family", "crosslink_intensity")))
            Assert.That(model.GetProperty("verb").GetString(), Is.EqualTo("predict crosslink-fragments"));
    }

    [Test]
    public void TheRetentionTimeUnitIsPerModel_NotPerFamily()
    {
        // The distinction that makes mzLib's bare IsIndexed boolean insufficient on the wire: one
        // model in the family returns absolute minutes and the rest return a dimensionless index.
        Assert.Multiple(() =>
        {
            Assert.That(ModelNamed("Prosit_2019_irt").GetProperty("retention_time_unit").GetString(),
                Is.EqualTo("indexed_retention_time"));
            Assert.That(ModelNamed("Chronologer_RT").GetProperty("retention_time_unit").GetString(),
                Is.EqualTo("minutes"));
            // ...and it is meaningless outside that family.
            Assert.That(ModelNamed("IM2Deep").GetProperty("retention_time_unit").ValueKind,
                Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void AConstraintDistinguishesNotApplicableFromRequiredAnyValue()
    {
        // The trap. mzLib expresses both as a nullable set: null means "this model has no such
        // input", empty means "required, any value". Passing the raw collection through would make
        // CID look permissive and HCD look impossible — backwards for both.
        JsonElement hcd = ModelNamed("Prosit_2020_intensity_HCD").GetProperty("collision_energy");
        JsonElement cid = ModelNamed("Prosit_2020_intensity_CID").GetProperty("collision_energy");

        Assert.Multiple(() =>
        {
            Assert.That(hcd.GetProperty("requirement").GetString(), Is.EqualTo("any_value_required"));
            Assert.That(cid.GetProperty("requirement").GetString(), Is.EqualTo("not_applicable"));
        });
    }

    [Test]
    public void AConstraintListsItsValuesWhenItHasThem()
    {
        JsonElement altimeter = ModelNamed("Altimeter_2024_intensities").GetProperty("collision_energy");
        int[] energies = altimeter.GetProperty("values").EnumerateArray().Select(v => v.GetInt32()).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(altimeter.GetProperty("requirement").GetString(), Is.EqualTo("one_of"));
            Assert.That(energies, Does.Contain(20).And.Contain(40));
            Assert.That(energies, Does.Not.Contain(45));
            // Sorted, so a caller can present them without re-sorting and two runs agree.
            Assert.That(energies, Is.Ordered);
        });
    }

    [Test]
    public void AModelThatAcceptsNoModificationsSaysSoWithAnEmptyList()
    {
        // An empty allowed-UNIMOD list is a real answer, not a missing one.
        JsonElement pfly = ModelNamed("pfly_2024_fine_tuned");
        Assert.That(pfly.GetProperty("allowed_unimod_ids").GetArrayLength(), Is.Zero);
        Assert.That(ModelNamed("Prosit_2019_irt").GetProperty("allowed_unimod_ids").GetArrayLength(),
            Is.GreaterThan(0));
    }

    [Test]
    public void ADynamicFragmentIonCountIsNullRatherThanMinusOne()
    {
        // mzLib uses -1 to mean "the count depends on the peptide". A -1 on the wire would be a
        // number someone plots.
        JsonElement unispec = ModelNamed("UniSpec").GetProperty("number_of_predicted_fragment_ions");
        JsonElement prosit = ModelNamed("Prosit_2020_intensity_HCD").GetProperty("number_of_predicted_fragment_ions");

        Assert.Multiple(() =>
        {
            Assert.That(unispec.ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(prosit.GetInt32(), Is.EqualTo(174));
        });
    }

    [Test]
    public void TheCatalogueNamesTheServiceItWouldCall()
    {
        // Koina is a public shared server, and a caller deciding how hard to hit it should be able
        // to see who they are hitting.
        Assert.That(Models().GetProperty("service").GetString(), Does.Contain("koina"));
    }

    [Test]
    public void AnUnknownFamilyReturnsNothingRatherThanEverything()
    {
        // A filter that matched nothing must never silently widen to the whole catalogue.
        Assert.That(Models("--family", "not_a_family").GetProperty("model_count").GetInt32(), Is.Zero);
    }

    [Test]
    public void AnEmptyFamilyOptionIsAUsageErrorNotAWiderSearch()
    {
        JsonElement error = InvokeExpectingError("predict", "models", "--family");
        Assert.That(error.GetProperty("message").GetString(), Does.Contain("empty"));
    }

    // ---- input handling ------------------------------------------------------------------------

    [TestCase("retention-time")]
    [TestCase("fragments")]
    [TestCase("ccs")]
    [TestCase("detectability")]
    [TestCase("crosslink-fragments")]
    public void AnUnknownModelNamesTheOnesThisVerbHas(string verb)
    {
        JsonElement error = InvokeExpectingError("predict", verb, "--model", "not_a_model", "--sequence", "PEPTIDEK");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
            Assert.That(error.GetProperty("message").GetString(), Does.Contain("Available:"),
                "a caller who guessed a name has no other way to find the right one");
        });
    }

    [Test]
    public void AModelFromTheWrongFamilyIsRejectedByName()
    {
        // Prosit_2019_irt is a real model, and not a fragment-intensity one. Accepting it would
        // fail much later, inside a request, with an error about the payload rather than the name.
        JsonElement error = InvokeExpectingError(
            "predict", "fragments", "--model", "Prosit_2019_irt", "--sequence", "PEPTIDEK");

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("No model named"));
    }

    [Test]
    public void AnUnknownInputColumnIsRefusedRatherThanIgnored()
    {
        // A caller who wrote 'charge' where the contract says 'precursor_charge' asked for
        // something specific; predicting at a default charge would hand back a plausible wrong
        // answer.
        JsonElement error = InvokeWithStdin(
            "charge\n2\n", "predict", "ccs", "--model", "IM2Deep");

        Assert.Multiple(() =>
        {
            Assert.That(error.GetProperty("ok").GetBoolean(), Is.False);
            Assert.That(error.GetProperty("error").GetProperty("message").GetString(),
                Does.Contain("Unknown input column"));
        });
    }

    [Test]
    public void AMissingRequiredColumnNamesIt()
    {
        JsonElement error = InvokeWithStdin("sequence\nPEPTIDEK\n", "predict", "ccs", "--model", "IM2Deep");

        Assert.That(error.GetProperty("error").GetProperty("message").GetString(),
            Does.Contain("precursor_charge"));
    }

    [Test]
    public void AHeaderWithNoRowsIsAUsageErrorNotAnEmptySuccess()
    {
        JsonElement error = InvokeWithStdin("sequence\n", "predict", "retention-time", "--model", "Prosit_2019_irt");

        Assert.That(error.GetProperty("error").GetProperty("message").GetString(),
            Does.Contain("no rows"));
    }

    [Test]
    public void NoInputAtAllExplainsBothWaysToProvideIt()
    {
        JsonElement error = InvokeWithStdin("", "predict", "retention-time", "--model", "Prosit_2019_irt");

        string message = error.GetProperty("error").GetProperty("message").GetString()!;
        Assert.Multiple(() =>
        {
            Assert.That(message, Does.Contain("--sequence"));
            Assert.That(message, Does.Contain("stdin"));
        });
    }

    [Test]
    public void ANonNumericChargeNamesTheLineAndTheColumn()
    {
        JsonElement error = InvokeWithStdin(
            "sequence\tprecursor_charge\nPEPTIDEK\ttwo\n", "predict", "ccs", "--model", "IM2Deep");

        string message = error.GetProperty("error").GetProperty("message").GetString()!;
        Assert.Multiple(() =>
        {
            Assert.That(message, Does.Contain("line 2"));
            Assert.That(message, Does.Contain("precursor_charge"));
        });
    }

    [TestCase("max-batches", "0")]
    [TestCase("throttle-ms", "-1")]
    public void APolitenessOptionBelowItsFloorIsRefused(string option, string value)
    {
        // The knobs exist for a large legitimate job. Letting one be set to nonsense would let a
        // caller hammer a shared community server by accident.
        JsonElement error = InvokeExpectingError(
            "predict", "retention-time", "--model", "Prosit_2019_irt", "--sequence", "PEPTIDEK",
            $"--{option}", value);

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("or greater"));
    }

    [TestCase("max-batches")]
    [TestCase("throttle-ms")]
    public void APolitenessOptionWithNoValueIsRefusedRatherThanDefaulted(string option)
    {
        JsonElement error = InvokeExpectingError(
            "predict", "retention-time", "--model", "Prosit_2019_irt", "--sequence", "PEPTIDEK",
            $"--{option}");

        Assert.That(error.GetProperty("message").GetString(), Does.Contain("no value"));
    }

    // ---- the columnar projection, without a network ---------------------------------------------
    //
    // The verbs' bodies past model resolution are pure: they read the same column sets exercised
    // here off a list of mzLib prediction records. Testing them directly rather than through a
    // request keeps the offline suite able to cover the shaping — which is the part this bridge
    // actually owns, the prediction itself being mzLib's.

    private static List<PeptideRTPrediction> TwoRetentionTimes() =>
    [
        new("PEPTIDEK", "PEPTIDEK", 5.5, IsIndexed: true),
        new("TOOLONG", "TOOLONG", null, IsIndexed: true,
            Warning: new System.ComponentModel.WarningException("Sequence is too long for this model.")),
    ];

    private static JsonElement Shaped<T>(
        (string Name, Func<T, object?> Read)[] columns, List<T> rows, params string[] extraArgs)
    {
        var arguments = new Program.Arguments(["predict", "retention-time", .. extraArgs]);
        object table = Prediction.Table(arguments, "a_model", rows.Count, columns, rows,
            extra: new { caveats = new[] { "a caveat" } });
        return JsonSerializer.SerializeToElement(table, Program.JsonOptions);
    }

    [Test]
    public void TheColumnarPayloadHasOneArrayPerFieldAndKeepsEveryRow()
    {
        JsonElement data = Shaped(Prediction.RetentionTimeColumns, TwoRetentionTimes());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("model").GetString(), Is.EqualTo("a_model"));
            Assert.That(data.GetProperty("row_count").GetInt32(), Is.EqualTo(2));
            foreach (JsonProperty column in data.GetProperty("columns").EnumerateObject())
                Assert.That(column.Value.GetArrayLength(), Is.EqualTo(2),
                    $"column '{column.Name}' must carry one value per peptide sent");
        });
    }

    [Test]
    public void APeptideThatCouldNotBePredictedIsCountedAndKeepsItsRow()
    {
        // The row must survive, or predictions stop lining up with the peptides that were sent.
        JsonElement data = Shaped(Prediction.RetentionTimeColumns, TwoRetentionTimes());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("failed_row_count").GetInt32(), Is.EqualTo(1));
            Assert.That(data.GetProperty("columns").GetProperty("retention_time")[1].ValueKind,
                Is.EqualTo(JsonValueKind.Null));
            Assert.That(data.GetProperty("columns").GetProperty("warning")[1].GetString(),
                Does.Contain("too long"));
            // ...and the row that DID predict is untouched.
            Assert.That(data.GetProperty("columns").GetProperty("retention_time")[0].GetDouble(),
                Is.EqualTo(5.5));
        });
    }

    [Test]
    public void TheVerbsOwnFieldsAreMergedWithTheSharedOnes()
    {
        JsonElement data = Shaped(Prediction.RetentionTimeColumns, TwoRetentionTimes());

        Assert.Multiple(() =>
        {
            Assert.That(data.TryGetProperty("column_names", out _), Is.True);
            Assert.That(data.GetProperty("caveats")[0].GetString(), Is.EqualTo("a caveat"));
        });
    }

    [Test]
    public void RaggedFragmentArraysCrossAsPerRowArrays()
    {
        // Koina returns a fixed-width grid with -1 for impossible ions and mzLib drops those, so
        // two peptides in one call genuinely differ in length. A rectangle would be a lie.
        List<PeptideFragmentIntensityPrediction> rows =
        [
            new("PEPTIDEK", "PEPTIDEK", 2, ["y1+1", "b1+1"], [147.1, 98.1], [0.9, 0.2]),
            new("ELVISLIVESK", "ELVISLIVESK", 2, ["y1+1"], [147.1], [1.0]),
        ];

        JsonElement data = Shaped(Prediction.FragmentColumns, rows);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").GetProperty("fragment_mz")[0].GetArrayLength(), Is.EqualTo(2));
            Assert.That(data.GetProperty("columns").GetProperty("fragment_mz")[1].GetArrayLength(), Is.EqualTo(1));
            // Index-aligned within a row, which is the only alignment that holds.
            Assert.That(data.GetProperty("columns").GetProperty("fragment_intensity")[0].GetArrayLength(),
                Is.EqualTo(data.GetProperty("columns").GetProperty("fragment_annotations")[0].GetArrayLength()));
        });
    }

    [Test]
    public void WritingToDiskOmitsTheInlinePayloadAndAgreesOnItsHeader()
    {
        string output = Path.Combine(_tempDirectory, "predictions.tsv");
        JsonElement data = Shaped(Prediction.RetentionTimeColumns, TwoRetentionTimes(), "--out", output);

        Assert.That(File.Exists(output), Is.True);
        string[] lines = File.ReadAllLines(output);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").ValueKind, Is.EqualTo(JsonValueKind.Null),
                "materialising both would defeat the point of --out");
            Assert.That(lines[0].Split('\t'),
                Is.EqualTo(data.GetProperty("column_names").EnumerateArray().Select(c => c.GetString()).ToArray()));
            Assert.That(lines, Has.Length.EqualTo(3), "a header and two rows");
        });
    }

    [Test]
    public void AWrittenArrayCellIsSemicolonJoined()
    {
        string output = Path.Combine(_tempDirectory, "fragments.tsv");
        List<PeptideFragmentIntensityPrediction> rows =
        [
            new("PEPTIDEK", "PEPTIDEK", 2, ["y1+1", "b1+1"], [147.1, 98.1], [0.9, 0.2]),
        ];

        Shaped(Prediction.FragmentColumns, rows, "--out", output);

        // ToString() on a List would write the type's name into the cell.
        Assert.That(File.ReadAllText(output), Does.Contain("147.1;98.1"));
    }

    [Test]
    public void AnEmptyOutIsRefusedRatherThanTreatedAsAbsent()
    {
        var arguments = new Program.Arguments(["predict", "retention-time", "--out"]);
        Assert.That(
            () => Prediction.Table(arguments, "m", 0, Prediction.RetentionTimeColumns,
                new List<PeptideRTPrediction>(), extra: new { }),
            Throws.TypeOf<Program.UsageException>());
    }

    [TestCase(null, "")]
    [TestCase(true, "true")]
    [TestCase(false, "false")]
    public void RenderWritesEachScalarAsInvariantText(object? value, string expected)
    {
        Assert.That(Prediction.Render(value), Is.EqualTo(expected));
    }

    [Test]
    public void RenderUsesTheInvariantCultureForNumbers()
    {
        // A comma-decimal locale would otherwise write "1,5" into a tab-separated file, which
        // reads back as a different number in every locale but its own.
        Assert.That(Prediction.Render(1.5), Is.EqualTo("1.5"));
    }

    [Test]
    public void EveryColumnSetNamesAWarningColumn()
    {
        // failed_row_count is derived from it, so a family that dropped the column would silently
        // report every row as successful.
        Assert.Multiple(() =>
        {
            Assert.That(Prediction.RetentionTimeColumns.Select(c => c.Name), Does.Contain("warning"));
            Assert.That(Prediction.FragmentColumns.Select(c => c.Name), Does.Contain("warning"));
            Assert.That(Prediction.CcsColumns.Select(c => c.Name), Does.Contain("warning"));
            Assert.That(Prediction.DetectabilityColumns.Select(c => c.Name), Does.Contain("warning"));
            Assert.That(Prediction.CrosslinkColumns.Select(c => c.Name), Does.Contain("warning"));
        });
    }

    [Test]
    public void TheDetectabilityColumnsProjectAllFourClasses()
    {
        List<PeptideDetectabilityPrediction> rows =
        [
            new("PEPTIDEK", "PEPTIDEK", (0.6, 0.3, 0.08, 0.02)),
            new("BADONE", "BADONE", null,
                new System.ComponentModel.WarningException("unsupported")),
        ];

        JsonElement data = Shaped(Prediction.DetectabilityColumns, rows);

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("columns").GetProperty("not_detectable")[0].GetDouble(), Is.EqualTo(0.6));
            Assert.That(data.GetProperty("columns").GetProperty("high_detectability")[0].GetDouble(), Is.EqualTo(0.02));
            // A peptide with no probabilities at all gets nulls across the four, not zeros.
            Assert.That(data.GetProperty("columns").GetProperty("not_detectable")[1].ValueKind,
                Is.EqualTo(JsonValueKind.Null));
        });
    }

    [Test]
    public void TheCcsAndCrosslinkColumnsProjectTheirOwnFields()
    {
        JsonElement ccs = Shaped(Prediction.CcsColumns,
            new List<PeptideCCSPrediction> { new("PEPTIDEK", "PEPTIDEK", 2, 327.5) });
        JsonElement crosslink = Shaped(Prediction.CrosslinkColumns,
            new List<CrosslinkFragmentIntensityPrediction>
            {
                new("ALPHA", "BETA", "ALPHA", "BETA", 3, ["y1+1"], [147.1], [1.0]),
            });

        Assert.Multiple(() =>
        {
            Assert.That(ccs.GetProperty("columns").GetProperty("collisional_cross_section")[0].GetDouble(),
                Is.EqualTo(327.5));
            Assert.That(crosslink.GetProperty("columns").GetProperty("beta_sequence")[0].GetString(),
                Is.EqualTo("BETA"));
        });
    }

    [Test]
    public void AModelIsResolvedByItsPublishedKoinaNameNotItsClassName()
    {
        // The published name is the identifier a user finds in the Koina catalogue and in the
        // literature; the class name (Prosit2019iRT) is an mzLib detail.
        var arguments = new Program.Arguments(["predict", "retention-time", "--model", "Prosit_2019_irt"]);
        var model = Prediction.Build<RetentionTimeModel>(arguments);

        Assert.That(model.ModelName, Is.EqualTo("Prosit_2019_irt"));
    }

    [Test]
    public void ModelNamesAreMatchedCaseInsensitively()
    {
        var arguments = new Program.Arguments(["predict", "retention-time", "--model", "prosit_2019_IRT"]);
        Assert.That(Prediction.Build<RetentionTimeModel>(arguments).ModelName, Is.EqualTo("Prosit_2019_irt"));
    }

    [Test]
    public void ThePolitenessKnobsReachTheModelWhenGiven()
    {
        var arguments = new Program.Arguments(
            ["predict", "retention-time", "--model", "Prosit_2019_irt",
             "--max-batches", "42", "--throttle-ms", "250"]);
        var model = Prediction.Build<RetentionTimeModel>(arguments);

        Assert.Multiple(() =>
        {
            Assert.That(model.MaxNumberOfBatchesPerRequest, Is.EqualTo(42));
            Assert.That(model.ThrottlingDelayInMilliseconds, Is.EqualTo(250));
        });
    }

    [Test]
    public void MzLibsOwnDefaultsSurviveWhenThePolitenessKnobsAreOmitted()
    {
        // Koina is a shared community GPU. The defaults are mzLib's and this must not raise them.
        var arguments = new Program.Arguments(["predict", "retention-time", "--model", "Prosit_2019_irt"]);
        var configured = Prediction.Build<RetentionTimeModel>(arguments);
        var untouched = new PredictionClients.Koina.SupportedModels.RetentionTimeModels.Prosit2019iRT();

        Assert.Multiple(() =>
        {
            Assert.That(configured.MaxNumberOfBatchesPerRequest,
                Is.EqualTo(untouched.MaxNumberOfBatchesPerRequest));
            Assert.That(configured.ThrottlingDelayInMilliseconds,
                Is.EqualTo(untouched.ThrottlingDelayInMilliseconds));
        });
    }

    [Test]
    public void AnIndexedModelSaysSoInItsCaveats()
    {
        var arguments = new Program.Arguments(["predict", "retention-time", "--model", "Prosit_2019_irt"]);
        string[] caveats = Prediction.RetentionTimeCaveats(Prediction.Build<RetentionTimeModel>(arguments));

        Assert.That(caveats.Any(caveat => caveat.Contains("INDEXED")), Is.True,
            "an iRT value plotted against a gradient is the mistake this caveat exists to prevent");
    }

    [Test]
    public void TheNetworkChronologerWarnsThatTheLocalOneDiffers()
    {
        // Same weights, different units: absolute retention time over the network, % acetonitrile
        // locally. Naming one "chronologer" without saying so would be ambiguous.
        var arguments = new Program.Arguments(["predict", "retention-time", "--model", "Chronologer_RT"]);
        string[] caveats = Prediction.RetentionTimeCaveats(Prediction.Build<RetentionTimeModel>(arguments));

        Assert.Multiple(() =>
        {
            Assert.That(caveats.Any(caveat => caveat.Contains("acetonitrile")), Is.True);
            Assert.That(caveats.Any(caveat => caveat.Contains("INDEXED")), Is.False,
                "Chronologer_RT is the one model in the family that returns real minutes");
        });
    }

    // ---- harness -------------------------------------------------------------------------------

    private static JsonElement Invoke(params string[] args)
    {
        JsonElement envelope = Envelope(args, stdin: null);
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.True, $"Expected success, got: {envelope}");
        return envelope.GetProperty("data");
    }

    private static JsonElement InvokeExpectingError(params string[] args)
    {
        JsonElement envelope = Envelope(args, stdin: null);
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.False, $"Expected a failure, got: {envelope}");
        return envelope.GetProperty("error");
    }

    /// <summary>Runs a verb with a table on stdin, returning the whole envelope.</summary>
    private static JsonElement InvokeWithStdin(string stdin, params string[] args) =>
        Envelope(args, stdin);

    private static JsonElement Envelope(string[] args, string? stdin)
    {
        TextReader previousIn = Console.In;
        TextWriter previousOut = Console.Out;
        if (stdin is not null)
            Console.SetIn(new StringReader(stdin));
        Console.SetOut(new StringWriter());

        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.SerializeToElement(new { ok = true, data }, Program.JsonOptions);
        }
        catch (Program.UsageException usage)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = "usage", message = usage.Message } },
                Program.JsonOptions);
        }
        catch (Exception exception)
        {
            return JsonSerializer.SerializeToElement(
                new
                {
                    ok = false,
                    error = new
                    {
                        type = Program.ClassifyError(exception),
                        message = Program.Unwrap(exception).Message,
                    },
                },
                Program.JsonOptions);
        }
        finally
        {
            Console.SetOut(previousOut);
            if (stdin is not null)
                Console.SetIn(previousIn);
        }
    }
}
