using System.Diagnostics.CodeAnalysis;
using System.Text.Json;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for <c>isobaric kits</c>, the projection of mzLib's <c>IsobaricMassTag</c>.
/// </summary>
/// <remarks>
/// The m/z values are mzLib's (derived from its embedded TMT.txt); these tests prove the table is
/// built from them faithfully, in order and paired with the right labels, and pin one value per
/// kit family against the published reporter-ion mass so a change to TMT.txt is noticed here.
/// </remarks>
[TestFixture]
[ExcludeFromCodeCoverage]
public class IsobaricTests
{
    [Test]
    public void Kits_ListsEveryKitMzLibKnows()
    {
        JsonElement data = Invoke("isobaric", "kits");

        string[] kits = data.GetProperty("kits").EnumerateArray().Select(k => k.GetProperty("kit").GetString()!).ToArray();
        int channels = data.GetProperty("kits").EnumerateArray().Sum(k => k.GetProperty("channel_count").GetInt32());

        Assert.Multiple(() =>
        {
            Assert.That(data.GetProperty("kit").ValueKind, Is.EqualTo(JsonValueKind.Null));
            Assert.That(kits, Is.EqualTo(new[] { "TMT6", "TMT10", "TMT11", "TMT16", "TMT18", "iTRAQ4", "iTRAQ8", "diLeu4", "diLeu12" }));
            Assert.That(data.GetProperty("kit_count").GetInt32(), Is.EqualTo(kits.Length));
            Assert.That(data.GetProperty("record_count").GetInt32(), Is.EqualTo(channels));
            Assert.That(data.GetProperty("absolute_tolerance").GetDouble(), Is.EqualTo(0.003));
            Assert.That(Strings(data.GetProperty("column_names")), Is.EqualTo(new[]
                { "kit", "channel_index", "channel_label", "reporter_ion_mz", "mz_min", "mz_max" }));
            Assert.That(data.GetProperty("caveats").GetArrayLength(), Is.EqualTo(5));
        });
    }

    [Test]
    public void Kits_Tmt18HasEighteenChannelsInAscendingMz()
    {
        JsonElement data = Invoke("isobaric", "kits", "--kit", "TMT18");
        JsonElement columns = data.GetProperty("columns");
        double[] mz = Doubles(columns.GetProperty("reporter_ion_mz"));
        double[] low = Doubles(columns.GetProperty("mz_min"));
        double[] high = Doubles(columns.GetProperty("mz_max"));

        Assert.That(data.GetProperty("kit").GetString(), Is.EqualTo("TMT18"));
        Assert.That(Strings(columns.GetProperty("channel_label")), Is.EqualTo(new[]
            { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
              "132N", "132C", "133N", "133C", "134N", "134C", "135N" }));
        Assert.That(Ints(columns.GetProperty("channel_index")), Is.EqualTo(Enumerable.Range(0, 18)));
        Assert.That(mz, Is.Ordered.Ascending);
        // The TMT 126 reporter ion's widely quoted m/z; a change to mzLib's TMT.txt shows up here.
        Assert.That(mz[0], Is.EqualTo(126.127726).Within(1e-5));
        for (int i = 0; i < mz.Length; i++)
        {
            Assert.That(low[i], Is.EqualTo(mz[i] - 0.003).Within(1e-9));
            Assert.That(high[i], Is.EqualTo(mz[i] + 0.003).Within(1e-9));
        }
    }

    [Test]
    public void Kits_Tmt16IsTheLowestSixteenOfTmt18()
    {
        double[] tmt16 = Doubles(Invoke("isobaric", "kits", "--kit", "TMT16").GetProperty("columns").GetProperty("reporter_ion_mz"));
        double[] tmt18 = Doubles(Invoke("isobaric", "kits", "--kit", "TMT18").GetProperty("columns").GetProperty("reporter_ion_mz"));

        Assert.That(tmt16, Is.EqualTo(tmt18.Take(16)));
    }

    [Test]
    public void Kits_Itraq8HasNo120Channel()
    {
        string[] labels = Strings(Invoke("isobaric", "kits", "--kit", "iTRAQ8").GetProperty("columns").GetProperty("channel_label"));

        Assert.That(labels, Is.EqualTo(new[] { "113", "114", "115", "116", "117", "118", "119", "121" }));
    }

    [TestCase("iTRAQ-4plex on K", "iTRAQ4")]
    [TestCase("tmt10", "TMT10")]
    [TestCase("TMT6-plex", "TMT6")]
    [TestCase("DiLeu-12plex", "diLeu12")]
    public void Kits_AcceptsMetaMorpheusModificationNames(string name, string kit)
    {
        JsonElement data = Invoke("isobaric", "kits", "--kit", name);

        Assert.That(data.GetProperty("kit").GetString(), Is.EqualTo(name), "the name as given is echoed");
        Assert.That(data.GetProperty("kits")[0].GetProperty("kit").GetString(), Is.EqualTo(kit));
    }

    [TestCase("TMT10plex")]
    [TestCase("TMT")]
    [TestCase("SILAC")]
    public void Kits_AnUnknownOrPartialNameIsAUsageError(string name)
    {
        JsonElement error = InvokeExpectingError("isobaric", "kits", "--kit", name);

        Assert.That(error.GetProperty("type").GetString(), Is.EqualTo("usage"));
        Assert.That(error.GetProperty("message").GetString(), Does.Contain(name).And.Contain("TMT18"));
    }

    [Test]
    public void Kits_ABlankKitIsAUsageError()
    {
        Assert.That(InvokeExpectingError("isobaric", "kits", "--kit").GetProperty("type").GetString(), Is.EqualTo("usage"));
    }

    // ---- helpers ------------------------------------------------------------------------------

    private static string[] Strings(JsonElement array) => array.EnumerateArray().Select(x => x.GetString()!).ToArray();

    private static int[] Ints(JsonElement array) => array.EnumerateArray().Select(x => x.GetInt32()).ToArray();

    private static double[] Doubles(JsonElement array) => array.EnumerateArray().Select(x => x.GetDouble()).ToArray();

    private static JsonElement Invoke(params string[] args) => Unwrap(Envelope(args), expectOk: true);

    private static JsonElement InvokeExpectingError(params string[] args) => Unwrap(Envelope(args), expectOk: false);

    private static JsonElement Unwrap(JsonElement envelope, bool expectOk)
    {
        Assert.That(envelope.GetProperty("ok").GetBoolean(), Is.EqualTo(expectOk), $"Unexpected envelope: {envelope}");
        return envelope.GetProperty(expectOk ? "data" : "error");
    }

    private static JsonElement Envelope(string[] args)
    {
        try
        {
            object data = Program.DispatchAsync(args).GetAwaiter().GetResult();
            return JsonSerializer.SerializeToElement(new { ok = true, data }, Program.JsonOptions);
        }
        catch (Program.UsageException usage)
        {
            return JsonSerializer.SerializeToElement(
                new { ok = false, error = new { type = "usage", message = usage.Message } }, Program.JsonOptions);
        }
    }
}
