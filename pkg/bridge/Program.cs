using System.Net.Http;
using System.Net.Sockets;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Text.RegularExpressions;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace MzLibBridge;

/// <summary>
/// A language-neutral command-line bridge over mzLib. Each invocation runs one verb and writes a
/// single JSON envelope to stdout: <c>{"ok":true,"data":…}</c> or <c>{"ok":false,"error":{…}}</c>.
/// </summary>
/// <remarks>
/// <para>
/// This executable is published self-contained and single-file, then staged inside the pyMzLib
/// Python wheel, so a consumer needs no .NET installation (decision D2).
/// </para>
/// <para>
/// The contract is deliberately <i>not</i> Python-shaped: stdin/stdout JSON with a stable envelope
/// is equally consumable from Rust, R, or a shell. Nothing about the wire format should assume its
/// caller is Python.
/// </para>
/// <para>
/// Diagnostics go to stderr. stdout carries the JSON envelope and nothing else, so a caller can
/// parse it without stripping log lines.
/// </para>
/// </remarks>
public static class Program
{
    /// <summary>The wire-format version. Bumped when the JSON envelope's shape changes incompatibly.</summary>
    private const int ProtocolVersion = 1;

    internal static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower,
        DefaultIgnoreCondition = JsonIgnoreCondition.Never,
        WriteIndented = false,
    };

    /// <summary>
    /// How a <see cref="PrideArchiveClient"/> is obtained. Tests replace this with one built over a
    /// stub <see cref="HttpMessageHandler"/> so the verb handlers can be exercised without touching
    /// the network — the same approach mzLib's own PrideArchiveClient tests use.
    /// </summary>
    internal static Func<PrideArchiveClient> PrideClientFactory { get; set; } = () => new PrideArchiveClient();

    /// <summary>
    /// The error type reported when an external service is unavailable, as opposed to something
    /// here being wrong. Callers key off this to decide whether retrying later is sensible; the
    /// test suites use it to skip rather than fail when EBI is down.
    /// </summary>
    internal const string ServiceUnavailableType = "ServiceUnavailable";

    /// <summary>
    /// Names the error type a failure should cross the boundary as.
    /// </summary>
    /// <remarks>
    /// The distinction worth drawing is availability versus correctness. A timeout, a dropped
    /// socket, or an HTTP 408/429/5xx means the service is having a bad day and the caller should
    /// try later. Anything else — a 404 from a wrong URL, a response that will not parse — means
    /// something is genuinely broken and silence would be misleading.
    /// <para>
    /// The status code is read out of the message first and from
    /// <see cref="HttpRequestException.StatusCode"/> only as a fallback, because mzLib's client
    /// composes its exception message by hand and leaves the property null.
    /// </para>
    /// </remarks>
    internal static string ClassifyError(Exception exception) => exception switch
    {
        TaskCanceledException or OperationCanceledException or TimeoutException => ServiceUnavailableType,
        SocketException => ServiceUnavailableType,
        HttpRequestException http => ClassifyHttpFailure(http),
        _ => exception.GetType().Name,
    };

    /// <summary>
    /// Classifies an HTTP failure as availability or correctness.
    /// </summary>
    /// <remarks>
    /// The case that is easy to get wrong, and was: a connection refused, a DNS failure, or a TLS
    /// handshake failure produces an <see cref="HttpRequestException"/> with <b>no status code at
    /// all</b>, because no response ever arrived. HttpClient wraps the underlying
    /// <see cref="SocketException"/>, so the socket arm above never sees it. Those are the most
    /// common outage modes there are — treating them as correctness failures would turn every
    /// real outage into a red build, which is exactly backwards.
    /// </remarks>
    private static string ClassifyHttpFailure(HttpRequestException exception)
    {
        int? status = StatusCodeOf(exception);
        if (status is not null)
            return IsAvailabilityStatus(status.Value) ? ServiceUnavailableType : nameof(HttpRequestException);

        // No status at all. Two very different things land here, and the inner exception tells
        // them apart: HttpClient wraps the real cause when a request never completed (a refused
        // connection, a DNS failure, a TLS handshake failure), whereas mzLib composes an
        // HttpRequestException by hand — the paging guard, for instance — to signal a condition
        // that is our problem and must not be excused as an outage.
        return exception.InnerException is not null ? ServiceUnavailableType : nameof(HttpRequestException);
    }

    /// <summary>An HTTP status that means "not now" rather than "not ever".</summary>
    private static bool IsAvailabilityStatus(int status) =>
        status is 408 or 429 || status >= 500;

    /// <summary>
    /// Recovers the status code from the exception's message, or failing that from its property.
    /// </summary>
    /// <remarks>
    /// mzLib composes its message by hand and leaves <see cref="HttpRequestException.StatusCode"/> null, so the message is authoritative here and is read FIRST; the property is the fallback for exceptions raised by HttpClient itself. The pattern is anchored to mzLib's own wording rather
    /// than matching "status NNN" anywhere in the text: an unanchored pattern can be steered by
    /// any caller-supplied string that reaches the message — an accession of
    /// <c>"x status 503 x"</c> would classify itself as an outage.
    /// </remarks>
    private static int? StatusCodeOf(HttpRequestException exception)
    {
        Match match = Regex.Match(
            exception.Message, @"failed with status (\d{3})\b", RegexOptions.IgnoreCase);
        if (match.Success)
            return int.Parse(match.Groups[1].Value);

        return exception.StatusCode.HasValue ? (int)exception.StatusCode.Value : null;
    }

    /// <summary>Runs one verb. Exit code 0 on success, 1 on a handled failure, 2 on bad usage.</summary>
    private static async Task<int> Main(string[] args)
    {
        try
        {
            object data = await DispatchAsync(args).ConfigureAwait(false);
            WriteJson(new Envelope { Ok = true, Data = data });
            return 0;
        }
        catch (UsageException ex)
        {
            WriteError("usage", ex.Message);
            return 2;
        }
        catch (Exception ex)
        {
            // Every failure crosses the boundary as structured data, never as a stack trace on
            // stdout: the caller's language has no way to interpret a .NET exception dump.
            //
            // Availability failures are labelled as such rather than left for the caller to
            // guess at. The distinction that matters to anyone consuming this is "the service is
            // down, try later" versus "something here is actually broken", and the bridge is the
            // only layer with enough information to tell them apart. Classifying here rather than
            // in a Python test helper means every consumer gets it — including a future binding
            // in another language.
            WriteError(ClassifyError(ex), ex.Message);
            return 1;
        }
    }

    /// <summary>Routes to the verb named by the leading positional arguments.</summary>
    internal static async Task<object> DispatchAsync(string[] args)
    {
        var arguments = new Arguments(args);

        return arguments.Verb switch
        {
            "version" => VersionInfo(),
            "pride files" => await PrideFilesAsync(arguments).ConfigureAwait(false),
            "pride download" => await PrideDownloadAsync(arguments).ConfigureAwait(false),
            _ => throw new UsageException(
                $"Unknown command '{arguments.Verb}'. Known commands: version, pride files, pride download."),
        };
    }

    /// <summary>Reports the bridge and protocol versions so a caller can check compatibility.</summary>
    private static object VersionInfo() => new
    {
        bridge = typeof(Program).Assembly.GetName().Version?.ToString() ?? "0.0.0",
        protocol = ProtocolVersion,
        runtime = Environment.Version.ToString(),
    };

    /// <summary>
    /// <c>pride files --accession PXD000001 [--page-size 100]</c> — the full file manifest of a
    /// PRIDE Archive project, with paging already resolved by mzLib.
    /// </summary>
    private static async Task<object> PrideFilesAsync(Arguments arguments)
    {
        string accession = arguments.Required("accession");
        int pageSize = arguments.OptionalInt("page-size", 100);

        using PrideArchiveClient client = PrideClientFactory();
        List<PrideArchiveFile> files = await client
            .GetProjectFilesAsync(accession, pageSize, CancellationToken.None)
            .ConfigureAwait(false);

        return new
        {
            accession,
            file_count = files.Count,
            total_size_bytes = files.TotalSizeBytes(),
            files = files.Select(ToWireFile).ToList(),
        };
    }

    /// <summary>
    /// <c>pride download --accession PXD000001 --dest DIR [--category RAW] [--ext .raw,.mzML]
    /// [--no-overwrite]</c> — downloads the selected files and reports where they landed.
    /// </summary>
    private static async Task<object> PrideDownloadAsync(Arguments arguments)
    {
        string accession = arguments.Required("accession");
        string destination = arguments.Required("dest");
        bool overwrite = !arguments.Flag("no-overwrite");

        // A filter that was ASKED FOR but degenerates to nothing must fail, never silently widen.
        // Previously a blank --category, or --ext "," , left the filter null and downloaded the
        // entire project — hundreds of megabytes instead of an error. A selection that cannot be
        // honoured is a usage error; quietly downloading everything is the worst possible reading
        // of the caller's intent.
        bool categoryRequested = arguments.WasProvided("category");
        bool extensionsRequested = arguments.WasProvided("ext");

        string? category = arguments.Optional("category");
        if (categoryRequested && string.IsNullOrWhiteSpace(category))
            throw new UsageException("Option --category was given but is empty; omit it to download every category.");

        string[] extensions = (arguments.Optional("ext") ?? string.Empty)
            .Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
        if (extensionsRequested && extensions.Length == 0)
            throw new UsageException("Option --ext was given but names no extensions; omit it to download every file type.");

        // An explicit selection arrives on stdin, one file name per line, rather than on argv:
        // argv has a hard ceiling of roughly 32 KB and a few thousand names exceed it, and any
        // separator we might choose could legitimately occur in a file name.
        // An explicit selection and a filter are contradictory instructions. Honouring one and
        // discarding the other silently would give the caller fewer files than either request
        // implies, with no indication which rule won.
        if (arguments.Flag("names-from-stdin") && (categoryRequested || extensionsRequested))
            throw new UsageException(
                "--names-from-stdin selects specific files, so --category and --ext cannot also be given.");

        HashSet<string>? selectedNames = null;
        if (arguments.Flag("names-from-stdin"))
        {
            selectedNames = new HashSet<string>(StringComparer.Ordinal);
            string? line;
            while ((line = Console.In.ReadLine()) != null)
            {
                if (!string.IsNullOrWhiteSpace(line))
                    selectedNames.Add(line.Trim());
            }

            if (selectedNames.Count == 0)
                throw new UsageException("--names-from-stdin was given but no file names were supplied.");
        }

        // Compose the filter from the mzLib extension methods rather than re-implementing the
        // selection logic here; the bridge stays a translation layer, not a second implementation.
        Func<PrideArchiveFile, bool>? filter = null;
        // Records which requested names the manifest actually contained, so a name that matched
        // nothing can be reported afterwards rather than silently reducing the download. A typo in
        // a selection previously produced fewer files and a success.
        HashSet<string>? matchedNames = null;

        if (selectedNames is not null)
        {
            HashSet<string> names = selectedNames;
            matchedNames = new HashSet<string>(StringComparer.Ordinal);
            filter = file =>
            {
                if (!names.Contains(file.FileName))
                    return false;
                matchedNames.Add(file.FileName);
                return true;
            };
        }
        else if (!string.IsNullOrWhiteSpace(category) || extensions.Length > 0)
        {
            filter = file =>
            {
                IEnumerable<PrideArchiveFile> single = new[] { file };
                if (!string.IsNullOrWhiteSpace(category) && !single.WhereCategory(category).Any())
                    return false;
                if (extensions.Length > 0 && !single.WhereExtension(extensions).Any())
                    return false;
                return true;
            };
        }

        using PrideArchiveClient client = PrideClientFactory();
        IReadOnlyList<string> paths = await client
            .DownloadProjectFilesAsync(accession, destination, filter, overwrite, CancellationToken.None)
            .ConfigureAwait(false);

        if (selectedNames is not null && matchedNames is not null && matchedNames.Count != selectedNames.Count)
        {
            List<string> missing = selectedNames.Except(matchedNames, StringComparer.Ordinal).Order().ToList();
            throw new UsageException(
                $"{missing.Count} of {selectedNames.Count} selected files are not in project '{accession}': " +
                $"{string.Join(", ", missing.Take(5))}{(missing.Count > 5 ? ", …" : string.Empty)}. " +
                "Selections must come from this project's own manifest.");
        }

        return new
        {
            accession,
            destination_directory = Path.GetFullPath(destination),
            downloaded_count = paths.Count,
            paths,
        };
    }

    /// <summary>
    /// Flattens a <see cref="PrideArchiveFile"/> into the wire shape. Controlled-vocabulary terms
    /// are reduced to their useful parts so a caller does not have to understand CV structure to
    /// use the manifest, while the full terms remain available for callers that do.
    /// </summary>
    internal static object ToWireFile(PrideArchiveFile file) => new
    {
        file_name = file.FileName,
        file_size_bytes = file.FileSizeBytes,
        checksum = file.Checksum,
        category = file.FileCategory.Value,
        category_accession = file.FileCategory.Accession,
        https_url = file.TryGetHttpsDownloadUrl(out string httpsUrl) ? httpsUrl : null,
        locations = file.PublicFileLocations.Select(ToWireCvParam).ToList(),
        submission_date = file.SubmissionDate,
        publication_date = file.PublicationDate,
        updated_date = file.UpdatedDate,
    };

    private static object ToWireCvParam(CvParam term) => new
    {
        accession = term.Accession,
        name = term.Name,
        value = term.Value,
    };

    private static void WriteJson(Envelope envelope) =>
        Console.Out.WriteLine(JsonSerializer.Serialize(envelope, JsonOptions));

    private static void WriteError(string type, string message) =>
        WriteJson(new Envelope { Ok = false, Error = new ErrorInfo { Type = type, Message = message } });

    /// <summary>The single object written to stdout by every invocation.</summary>
    internal sealed class Envelope
    {
        public bool Ok { get; init; }
        public object? Data { get; init; }
        public ErrorInfo? Error { get; init; }
    }

    /// <summary>A failure, described in terms a non-.NET caller can act on.</summary>
    internal sealed class ErrorInfo
    {
        public string Type { get; init; } = string.Empty;
        public string Message { get; init; } = string.Empty;
    }

    /// <summary>Signals malformed or missing command-line input (exit code 2), not a runtime failure.</summary>
    internal sealed class UsageException(string message) : Exception(message);

    /// <summary>
    /// A minimal parser for <c>&lt;verb…&gt; --name value --flag</c>. Deliberately hand-rolled: the
    /// bridge takes no third-party dependency it does not need, because every dependency is weight
    /// inside the published wheel.
    /// </summary>
    internal sealed class Arguments
    {
        private readonly Dictionary<string, string> _named = new(StringComparer.OrdinalIgnoreCase);
        private readonly HashSet<string> _flags = new(StringComparer.OrdinalIgnoreCase);

        /// <summary>The leading positional words joined by spaces, e.g. "pride files".</summary>
        public string Verb { get; }

        public Arguments(string[] args)
        {
            var positional = new List<string>();
            int i = 0;

            for (; i < args.Length && !args[i].StartsWith("--", StringComparison.Ordinal); i++)
                positional.Add(args[i]);

            for (; i < args.Length; i++)
            {
                string name = args[i].TrimStart('-');
                if (i + 1 < args.Length && !args[i + 1].StartsWith("--", StringComparison.Ordinal))
                    _named[name] = args[++i];
                else
                    _flags.Add(name);
            }

            Verb = string.Join(' ', positional);
        }

        public string Required(string name) =>
            _named.TryGetValue(name, out string? value) && !string.IsNullOrWhiteSpace(value)
                ? value
                : throw new UsageException($"Missing required option --{name}.");

        public string? Optional(string name) => _named.TryGetValue(name, out string? value) ? value : null;

        public bool Flag(string name) => _flags.Contains(name);

        /// <summary>
        /// Whether the option appeared at all, with or without a value.
        /// </summary>
        /// <remarks>
        /// Distinguishing "not asked for" from "asked for, but empty" is what stops a degenerate
        /// filter from silently widening to everything. An option written without a value —
        /// <c>--category --ext .raw</c> — lands in the flag set rather than the named set, so both
        /// collections have to be consulted or the option looks absent and the caller's intent is
        /// discarded without complaint.
        /// </remarks>
        public bool WasProvided(string name) => _named.ContainsKey(name) || _flags.Contains(name);

        public int OptionalInt(string name, int fallback)
        {
            string? raw = Optional(name);
            if (raw == null)
                return fallback;
            if (!int.TryParse(raw, out int parsed))
                throw new UsageException($"Option --{name} must be an integer; got '{raw}'.");
            return parsed;
        }
    }
}
