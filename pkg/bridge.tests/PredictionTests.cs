using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

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
