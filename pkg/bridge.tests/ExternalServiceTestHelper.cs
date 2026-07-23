using System.Diagnostics.CodeAnalysis;

namespace MzLibBridge.Tests;

/// <summary>
/// Support for tests that exercise a live external web service. Adapted from mzLib's helper of the
/// same name so both repositories behave identically when EBI has a bad day.
/// </summary>
/// <remarks>
/// <para>
/// The hard part is telling two failures apart:
/// </para>
/// <list type="bullet">
///   <item>the service is unavailable (down, rate-limited, 5xx, timeout) — <b>not</b> our bug, so
///     the test should be <b>skipped</b>: "we tried, the service is down, don't worry"; versus</item>
///   <item>the service answered but the contract is broken (wrong URL, response no longer parses,
///     an expected value missing) — a real regression that must <b>fail</b>.</item>
/// </list>
/// <para>
/// Here the classification has already happened one layer down: the bridge labels availability
/// failures <c>ServiceUnavailable</c> in its error envelope, so every consumer gets the
/// distinction rather than just this test suite. This helper turns that label into a skip.
/// </para>
/// </remarks>
[ExcludeFromCodeCoverage]
public static class ExternalServiceTestHelper
{
    /// <summary>
    /// Runs <paramref name="testBody"/>; if the external service proves unavailable, marks the test
    /// Skipped rather than Failed. Genuine assertion failures propagate.
    /// </summary>
    public static async Task RunAsync(string serviceName, Func<Task> testBody)
    {
        try
        {
            await testBody();
        }
        catch (ExternalServiceUnavailableException e)
        {
            Skip(serviceName, e.Message);
        }
        catch (HttpRequestException e)
        {
            Skip(serviceName, $"unavailable ({e.Message})");
        }
        catch (TaskCanceledException)
        {
            Skip(serviceName, "timed out");
        }
    }

    /// <summary>
    /// Writes the reason to the CI log — where it is visible immediately, regardless of verbosity —
    /// and then skips. This is what surfaces "we tried, the service is down" instead of leaving
    /// someone to wonder whether they broke something.
    /// </summary>
    private static void Skip(string serviceName, string reason)
    {
        string message = $"Skipping external-service test: {serviceName} {reason}. " +
                         "This is a third-party availability problem, not a code failure.";
        TestContext.Progress.WriteLine(message);
        Assert.Ignore(message);
    }
}

/// <summary>
/// Marker meaning "the external service is unavailable", as opposed to a real bug.
/// <see cref="ExternalServiceTestHelper.RunAsync"/> turns this into a skipped test.
/// </summary>
[ExcludeFromCodeCoverage]
public class ExternalServiceUnavailableException(string message) : Exception(message);
