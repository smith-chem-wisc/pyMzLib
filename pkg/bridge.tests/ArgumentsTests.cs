using MzLibBridge;

namespace MzLibBridge.Tests;

/// <summary>
/// Tests for the command-line parser. It is hand-rolled — the bridge takes no dependency it does
/// not need, because every dependency is weight inside the published wheel — so the usual
/// argument-parsing edge cases are ours to get right rather than a library's.
/// </summary>
[TestFixture]
public class ArgumentsTests
{
    [Test]
    public void LeadingWordsBecomeTheVerb()
    {
        var arguments = new Program.Arguments(["pride", "files", "--accession", "PXD000001"]);
        Assert.That(arguments.Verb, Is.EqualTo("pride files"));
    }

    [Test]
    public void ASingleWordIsAValidVerb()
    {
        Assert.That(new Program.Arguments(["version"]).Verb, Is.EqualTo("version"));
    }

    [Test]
    public void NoArgumentsYieldsAnEmptyVerbRatherThanThrowing()
    {
        // Dispatch turns this into a usage error with the list of known commands, which is more
        // helpful than a parser exception with no context.
        Assert.That(new Program.Arguments([]).Verb, Is.Empty);
    }

    [Test]
    public void NamedOptionsAreRead()
    {
        var arguments = new Program.Arguments(["pride", "files", "--accession", "PXD000001"]);
        Assert.That(arguments.Required("accession"), Is.EqualTo("PXD000001"));
    }

    [Test]
    public void OptionNamesAreCaseInsensitive()
    {
        var arguments = new Program.Arguments(["v", "--Accession", "PXD1"]);
        Assert.That(arguments.Required("accession"), Is.EqualTo("PXD1"));
    }

    [Test]
    public void AnOptionWithNoValueIsAFlag()
    {
        var arguments = new Program.Arguments(["pride", "download", "--no-overwrite"]);
        Assert.Multiple(() =>
        {
            Assert.That(arguments.Flag("no-overwrite"), Is.True);
            Assert.That(arguments.Flag("overwrite"), Is.False);
        });
    }

    [Test]
    public void AFlagFollowedByAnotherOptionDoesNotSwallowIt()
    {
        // The bug this guards: treating "--dest" as the value of "--no-overwrite" would both
        // lose the flag and lose the destination.
        var arguments = new Program.Arguments(["pride", "download", "--no-overwrite", "--dest", "out"]);
        Assert.Multiple(() =>
        {
            Assert.That(arguments.Flag("no-overwrite"), Is.True);
            Assert.That(arguments.Required("dest"), Is.EqualTo("out"));
        });
    }

    [Test]
    public void ValuesMayContainSpacesAndPathSeparators()
    {
        var arguments = new Program.Arguments(["pride", "download", "--dest", @"C:\my data\out"]);
        Assert.That(arguments.Required("dest"), Is.EqualTo(@"C:\my data\out"));
    }

    [Test]
    public void MissingRequiredOptionThrowsUsageNotSomethingElse()
    {
        // Must be UsageException specifically: it is what produces exit code 2 and surfaces in
        // Python as UsageError rather than as a runtime failure.
        var arguments = new Program.Arguments(["pride", "files"]);
        var ex = Assert.Throws<Program.UsageException>(() => arguments.Required("accession"));
        Assert.That(ex!.Message, Does.Contain("--accession"));
    }

    [Test]
    public void AnEmptyRequiredValueCountsAsMissing()
    {
        var arguments = new Program.Arguments(["pride", "files", "--accession", "   "]);
        Assert.Throws<Program.UsageException>(() => arguments.Required("accession"));
    }

    [Test]
    public void OptionalReturnsNullWhenAbsent()
    {
        Assert.That(new Program.Arguments(["v"]).Optional("category"), Is.Null);
    }

    [Test]
    public void OptionalIntFallsBackWhenAbsent()
    {
        Assert.That(new Program.Arguments(["v"]).OptionalInt("page-size", 100), Is.EqualTo(100));
    }

    [Test]
    public void OptionalIntParsesWhenPresent()
    {
        var arguments = new Program.Arguments(["v", "--page-size", "500"]);
        Assert.That(arguments.OptionalInt("page-size", 100), Is.EqualTo(500));
    }

    [Test]
    public void OptionalIntRejectsNonNumbersWithAUsableMessage()
    {
        var arguments = new Program.Arguments(["v", "--page-size", "lots"]);
        var ex = Assert.Throws<Program.UsageException>(() => arguments.OptionalInt("page-size", 100));
        Assert.That(ex!.Message, Does.Contain("lots"), "the bad value should appear in the message");
    }

    [Test]
    public void SingleDashIsAcceptedAsWellAsDouble()
    {
        var arguments = new Program.Arguments(["v", "--accession", "PXD1"]);
        Assert.That(arguments.Optional("accession"), Is.EqualTo("PXD1"));
    }
}
