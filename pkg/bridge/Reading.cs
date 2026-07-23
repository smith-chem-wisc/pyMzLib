using MassSpectrometry;
using Readers;

namespace MzLibBridge;

/// <summary>
/// Reading proteomics result files: what a file <i>is</i>, and what can be done with it.
/// </summary>
/// <remarks>
/// <para>
/// mzLib's <c>Readers</c> recognises 29 file types written by a dozen different search and
/// deconvolution tools, and dispatches each to a parser it maintains. That dispatch is the whole
/// value here — the bridge adds no parsing of its own; it asks mzLib what a path is
/// (<see cref="SupportedFileTypeExtensions.ParseFileType"/>) and reports the answer.
/// </para>
/// <para>
/// The important honesty in this verb is <c>views</c>. It is tempting to describe mzLib as reading
/// 29 formats into one uniform shape; it does not. Only three types implement
/// <see cref="IQuantifiableResultFile"/> (MetaMorpheus <c>.psmtsv</c> and <c>.osmtsv</c>, and
/// MSFragger <c>psm.tsv</c>); deconvolution features have their own interface
/// (<see cref="IMs1FeatureFile"/>); spectra files are an <see cref="MsDataFile"/>; and several
/// formats belong to no common interface at all. Rather than discovering that as a cast failure
/// deep in a later call, a caller is told up front which views a given file actually supports.
/// </para>
/// </remarks>
internal static class Reading
{
    /// <summary>
    /// <c>readers formats</c> — every file type mzLib can recognise, with its extension and the
    /// views it supports.
    /// </summary>
    /// <remarks>
    /// Enumerated from <see cref="SupportedFileType"/> rather than transcribed, so the published
    /// table of supported formats cannot drift from what mzLib actually dispatches. This mirrors
    /// how mzLib's own test suite enumerates the enum instead of maintaining a list.
    /// </remarks>
    public static object Formats(Program.Arguments arguments)
    {
        var formats = new List<object>();

        foreach (SupportedFileType fileType in Enum.GetValues<SupportedFileType>())
        {
            // A member missing its extension or reader mapping is a broken mzLib build, not a
            // reason to fail the whole listing: report what is known and keep going, so this verb
            // stays usable for diagnosing exactly that.
            string? extension = Try(() => fileType.GetFileExtension());
            Type? readerType = Try(() => fileType.GetResultFileType());

            formats.Add(new
            {
                file_type = fileType.ToString(),
                extension,
                reader = readerType?.Name,
                views = readerType is null ? new List<string>() : ViewsOf(readerType),
            });
        }

        return new
        {
            format_count = formats.Count,
            formats,
        };
    }

    /// <summary>
    /// <c>readers identify --path FILE</c> — what kind of result file this is, and what can be
    /// done with it, without parsing its contents.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Deliberately cheap: <see cref="FileReader.ReadResultFile"/> is lazy — it resolves the type
    /// and sets the path but does not call <c>LoadResults()</c> — so identifying a million-row file
    /// costs no more than identifying an empty one.
    /// </para>
    /// <para>
    /// "Without parsing" is not quite "without reading". mzLib's own type detection opens the file
    /// for three families: a bare <c>.tsv</c> is disambiguated by its first line, a <c>.mztab</c>
    /// by its first five, and a <c>.d</c> by which analysis file the directory contains. So this
    /// verb can fail on an unreadable file, which is why the path is checked first.
    /// </para>
    /// <para>
    /// <b>No <c>software</c> field, deliberately.</b> <see cref="IResultFile.Software"/> looks like
    /// the answer to "which tool wrote this" and is not. Every reader carries its software constant
    /// on its <i>path-taking</i> constructor (<c>MsFraggerPsmFile(string) : base(filePath,
    /// Software.MsFragger)</c>), but <see cref="FileReader"/> builds readers through the
    /// parameterless one, which sets <see cref="Software.Unspecified"/> — so the property is
    /// <c>Unspecified</c> for everything mzLib's own factory returns. Nor is it reliable through
    /// the other constructor: <c>MsFraggerPeptideFile</c> passes <c>Unspecified</c> outright and the
    /// psmtsv readers route to a base overload that never sets it. Reconstructing the value here —
    /// by probing the other constructor, or by mapping file type to tool — would mean the bridge
    /// answering a question mzLib cannot, which is precisely the drift the translation-layer rule
    /// exists to prevent. <c>file_type</c> and <c>reader</c> already name the tool
    /// (<c>MsFraggerPsm</c>, <c>ToppicPrsm</c>), so nothing is lost but a false certainty.
    /// </para>
    /// </remarks>
    public static object Identify(Program.Arguments arguments)
    {
        string path = arguments.Required("path");

        // Checked here rather than left to mzLib because FileReader throws a bare
        // FileNotFoundException carrying neither a message nor a file name, which would reach the
        // caller as an empty error. A missing path is the caller's mistake, so it is a usage
        // failure (exit 2), not a runtime one.
        // Directory.Exists as well as File.Exists: Bruker's .d "file" is a directory.
        if (!File.Exists(path) && !Directory.Exists(path))
            throw new Program.UsageException($"File not found: '{path}'.");

        IResultFile resultFile;
        try
        {
            resultFile = FileReader.ReadResultFile(path);
        }
        catch (MzLibUtil.MzLibException exception)
        {
            // mzLib recognises the file's extension or it does not; there is no sentinel and no
            // "unknown" enum member to return. An unsupported file is a bad argument rather than a
            // fault, so it surfaces as a usage error pointing at the verb that lists what IS
            // supported — not as a BridgeError, which would read as "pyMzLib is broken".
            throw new Program.UsageException(
                $"{exception.Message}: '{path}'. Run 'readers formats' for the file types mzLib recognises.");
        }

        Type readerType = resultFile.GetType();

        return new
        {
            path = Path.GetFullPath(path),
            file_type = resultFile.FileType.ToString(),
            extension = Try(() => resultFile.FileType.GetFileExtension()),
            reader = readerType.Name,
            views = ViewsOf(readerType),
        };
    }

    /// <summary>
    /// The common interfaces a reader implements, named as capabilities rather than as .NET types.
    /// </summary>
    /// <remarks>
    /// An empty list is a real and common answer — TopPIC, MsPathFinderT, Crux and Casanovo readers
    /// each parse into their own record type and share no cross-format interface. Saying so is the
    /// point of the verb.
    /// </remarks>
    private static List<string> ViewsOf(Type readerType)
    {
        var views = new List<string>();

        if (typeof(IQuantifiableResultFile).IsAssignableFrom(readerType))
            views.Add("quantifiable");
        if (typeof(IMs1FeatureFile).IsAssignableFrom(readerType))
            views.Add("ms1_features");
        if (typeof(MsDataFile).IsAssignableFrom(readerType))
            views.Add("spectra");

        // ISpectralMatch is implemented by the RECORD, not the file, so it cannot be found on the
        // reader itself — it is the type argument of the ResultFile<T> the reader derives from.
        Type? recordType = RecordTypeOf(readerType);
        if (recordType is not null && typeof(ISpectralMatch).IsAssignableFrom(recordType))
            views.Add("spectral_match");

        return views;
    }

    /// <summary>The <c>T</c> of the <c>ResultFile&lt;T&gt;</c> a reader derives from, if it does.</summary>
    private static Type? RecordTypeOf(Type readerType)
    {
        for (Type? type = readerType; type is not null; type = type.BaseType)
        {
            if (type.IsGenericType && type.GetGenericTypeDefinition() == typeof(ResultFile<>))
                return type.GetGenericArguments()[0];
        }

        return null;
    }

    /// <summary>
    /// Evaluates a lookup that throws rather than returning null for an unmapped value.
    /// </summary>
    /// <remarks>
    /// Used only by <see cref="Formats"/>, where one unmapped enum member must not take the whole
    /// listing down with it. Not a general-purpose swallow: every other path here lets failures
    /// propagate so the caller gets the real exception type.
    /// </remarks>
    private static T? Try<T>(Func<T> lookup) where T : class
    {
        try
        {
            return lookup();
        }
        catch (MzLibUtil.MzLibException)
        {
            return null;
        }
    }
}
