using Readers;

namespace MzLibBridge;

/// <summary>
/// Three questions about SDRF documents, each answered by the mzLib type that owns it:
/// <c>sdrf validate</c> (is one document well-formed — <see cref="SdrfValidator"/>),
/// <c>sdrf lint</c> (do several documents write the same thing the same way —
/// <see cref="SdrfDriftLint"/>) and <c>sdrf assess</c> (does one document's sample half say
/// anything — <see cref="SdrfSampleInformativeness"/>).
/// </summary>
/// <remarks>
/// <para>
/// <b>They are three instruments because each is blind where another sees.</b> A document of
/// reserved words passes the validator (reserved words are the specification's correct way to say
/// "no value") and the drift lint (a reserved word cannot drift) and still answers no biological
/// question; only <c>assess</c> sees that. Twenty individually valid documents can still be
/// unpoolable; only <c>lint</c> looks across files. And <c>assess</c> never looks at structure at
/// all. The caveats on each verb say which blind spot the caller is standing in.
/// </para>
/// <para>
/// <b>No rule is re-derived here.</b> Every finding, severity, verdict and count is mzLib's own,
/// shaped into a table and nothing more. A second implementation in the bridge — and then again in
/// Rust and R — is exactly the per-binding repair this bridge exists to prevent
/// (<c>bridge/design/PRINCIPLE.md</c>).
/// </para>
/// </remarks>
internal static partial class Sdrf
{
    private static readonly string[] ValidationColumns =
        { "severity", "rule", "message", "row_index", "line_number", "column_name" };

    private static readonly string[] AssessmentColumns =
        { "role", "column_name", "rows", "filled", "absent", "distinct_values", "fill_rate" };

    // ---- sdrf validate --------------------------------------------------------------------------

    /// <summary>
    /// <c>sdrf validate (--path FILE | --paths-stdin) [--threads N] [--on-error fail|skip]</c> —
    /// the structural findings mzLib's <see cref="SdrfValidator"/> reports for each document, one
    /// row per finding.
    /// </summary>
    public static object Validate(Program.Arguments arguments)
    {
        BulkInput input = BulkFrom(arguments);
        List<BulkOutcome<ValidatedDocument>> outcomes = RunEach(input, ValidateOne);

        if (!input.IsBulk)
        {
            ValidatedDocument only = outcomes[0].Value!;
            var table = new ColumnTable(ValidationColumns);
            foreach (SdrfValidationMessage m in only.Result.Messages)
                table.Add(MessageCells(m));

            return new
            {
                path = only.Path,
                is_valid = only.Result.IsValid,
                error_count = only.Result.Errors.Count(),
                warning_count = only.Result.Warnings.Count(),
                message_count = only.Result.Messages.Count,
                row_count = only.RowCount,
                column_names = table.Names,
                columns = table.Columns(),
                caveats = ValidateCaveats(),
            };
        }

        var bulk = ColumnTable.WithSource(ValidationColumns);
        var files = new List<object>(outcomes.Count);
        for (int i = 0; i < outcomes.Count; i++)
        {
            BulkOutcome<ValidatedDocument> outcome = outcomes[i];
            ValidatedDocument? v = outcome.Value;
            if (v is not null)
                foreach (SdrfValidationMessage m in v.Result.Messages)
                    bulk.Add(new object?[] { i, outcome.Path }.Concat(MessageCells(m)).ToArray());

            files.Add(new
            {
                path = outcome.Path,
                is_valid = v?.Result.IsValid,
                error_count = v?.Result.Errors.Count(),
                warning_count = v?.Result.Warnings.Count(),
                message_count = v?.Result.Messages.Count,
                row_count = v?.RowCount,
                error = outcome.Error,
            });
        }

        var read = outcomes.Where(o => o.Value is not null).Select(o => o.Value!).ToList();
        return new
        {
            file_count = outcomes.Count,
            read_count = read.Count,
            failed_count = outcomes.Count - read.Count,
            valid_count = read.Count(v => v.Result.IsValid),
            message_count = bulk.RowCount,
            // BULK.md section 2: every bulk envelope carries record_count, the rows in columns.
            record_count = bulk.RowCount,
            files,
            column_names = bulk.Names,
            columns = bulk.Columns(),
            caveats = ValidateCaveats(),
        };
    }

    private sealed record ValidatedDocument(string Path, int RowCount, SdrfValidationResult Result);

    private static ValidatedDocument ValidateOne(string path)
    {
        SdrfDocument document = Load(path);
        return new ValidatedDocument(path, document.Results.Count, SdrfValidator.Validate(document));
    }

    private static object?[] MessageCells(SdrfValidationMessage m) => new object?[]
    {
        m.Severity.ToString(), m.Rule, m.Message, m.RowIndex, m.LineNumber, m.ColumnName,
    };

    private static List<string> ValidateCaveats() => new()
    {
        "STRUCTURE ONLY, against SDRF-Proteomics v1.1.0: column presence, order, casing, cell shape " +
        "and row-key uniqueness. Controlled-vocabulary terms are NOT resolved, so a well-formed but " +
        "wrong accession passes.",

        "is_valid means no Error. Warnings never make a document invalid: mzLib calibrated every " +
        "severity against the 1,236-file curated corpus, where a rule that fired on most curated " +
        "files was judged a wrong rule, not 1,236 wrong files.",

        "row_index is 0-based over data rows; line_number is the 1-based line in the file " +
        "(row_index + 2, because the header is line 1). Both are null for a document-level finding. " +
        "column_name is null when the finding is not about one column.",

        "A document full of reserved words (\"not available\") validates cleanly, because that is " +
        "the specification's correct way to say nothing. Whether a document SAYS anything is sdrf " +
        "assess; whether several documents agree is sdrf lint.",
    };

    // ---- sdrf assess ----------------------------------------------------------------------------

    /// <summary>
    /// <c>sdrf assess (--path FILE | --paths-stdin) [--threads N] [--on-error fail|skip]</c> —
    /// mzLib's verdict on whether each document's sample half describes an experimental design
    /// (Informative / Partial / Skeleton), with the per-column counts the verdict was decided from.
    /// </summary>
    public static object Assess(Program.Arguments arguments)
    {
        BulkInput input = BulkFrom(arguments);
        List<BulkOutcome<AssessedDocument>> outcomes = RunEach(input, AssessOne);

        if (!input.IsBulk)
        {
            AssessedDocument only = outcomes[0].Value!;
            var table = new ColumnTable(AssessmentColumns);
            foreach (object?[] row in EvidenceRows(only.Assessment))
                table.Add(row);

            return new
            {
                path = only.Path,
                verdict = only.Assessment.Verdict.ToString(),
                factor_value_varies = only.Assessment.FactorValueVaries,
                sample_is_described = only.Assessment.SampleIsDescribed,
                biological_replicate_varies = only.Assessment.BiologicalReplicateVaries,
                row_count = only.RowCount,
                column_names = table.Names,
                columns = table.Columns(),
                caveats = AssessCaveats(),
            };
        }

        var bulk = ColumnTable.WithSource(AssessmentColumns);
        var files = new List<object>(outcomes.Count);
        for (int i = 0; i < outcomes.Count; i++)
        {
            BulkOutcome<AssessedDocument> outcome = outcomes[i];
            SdrfSampleAssessment? a = outcome.Value?.Assessment;
            if (a is not null)
                foreach (object?[] row in EvidenceRows(a))
                    bulk.Add(new object?[] { i, outcome.Path }.Concat(row).ToArray());

            files.Add(new
            {
                path = outcome.Path,
                verdict = a?.Verdict.ToString(),
                factor_value_varies = a?.FactorValueVaries,
                sample_is_described = a?.SampleIsDescribed,
                biological_replicate_varies = a?.BiologicalReplicateVaries,
                row_count = outcome.Value?.RowCount,
                error = outcome.Error,
            });
        }

        var verdicts = outcomes.Where(o => o.Value is not null).Select(o => o.Value!.Assessment.Verdict).ToList();
        return new
        {
            file_count = outcomes.Count,
            read_count = verdicts.Count,
            failed_count = outcomes.Count - verdicts.Count,
            // BULK.md section 2: every bulk envelope carries record_count, the rows in columns.
            record_count = bulk.RowCount,
            verdict_counts = new
            {
                informative = verdicts.Count(v => v == SdrfSampleVerdict.Informative),
                partial = verdicts.Count(v => v == SdrfSampleVerdict.Partial),
                skeleton = verdicts.Count(v => v == SdrfSampleVerdict.Skeleton),
            },
            files,
            column_names = bulk.Names,
            columns = bulk.Columns(),
            caveats = AssessCaveats(),
        };
    }

    private sealed record AssessedDocument(string Path, int RowCount, SdrfSampleAssessment Assessment);

    private static AssessedDocument AssessOne(string path)
    {
        SdrfDocument document = Load(path);
        return new AssessedDocument(path, document.Results.Count, SdrfSampleInformativeness.Assess(document));
    }

    /// <summary>
    /// The evidence behind a verdict: one row per column each check read, in mzLib's order —
    /// factor values, then sample characteristics, then biological replicate.
    /// </summary>
    private static IEnumerable<object?[]> EvidenceRows(SdrfSampleAssessment a)
    {
        foreach (SdrfColumnCoverage c in a.FactorValueColumns)
            yield return CoverageCells("factor_value", c);
        foreach (SdrfColumnCoverage c in a.SampleCharacteristicColumns)
            yield return CoverageCells("sample_characteristic", c);
        if (a.BiologicalReplicate is not null)
            yield return CoverageCells("biological_replicate", a.BiologicalReplicate);
    }

    private static object?[] CoverageCells(string role, SdrfColumnCoverage c) => new object?[]
    {
        role, c.Column, c.Rows, c.Filled, c.Absent, c.DistinctValues, c.FillRate,
    };

    private static List<string> AssessCaveats() => new()
    {
        "Informative = all three checks pass (a factor value varies, a sample characteristic other " +
        "than organism and biological replicate is filled, the biological replicates differ); " +
        "Skeleton = none pass; Partial = anything between. On the 1,236-file bigbio curated corpus " +
        "mzLib measured 300 Informative, 742 Partial and 194 Skeleton.",

        "Partial is often legitimate: a single-condition study has no factor to vary. Whether Partial " +
        "is good enough is the caller's decision; the three booleans say which check failed.",

        "A reserved word (\"not available\", \"not applicable\", \"anonymized\", \"pooled\") and an " +
        "empty cell both count as ABSENT, so they neither fill a column nor add a distinct value. " +
        "Distinct values are compared case-insensitively after trimming: \"liver\" and \"Liver\" are one.",

        "Rows are counted once per row per column NAME, so a repeated column does not inflate them. " +
        "This verb looks only at the sample half; it says nothing about structure (sdrf validate) or " +
        "agreement across documents (sdrf lint).",
    };

    // ---- sdrf lint ------------------------------------------------------------------------------

    private static readonly string[] LintColumns =
        { "finding_index", "kind", "concept", "column_name", "variant_rank", "value", "occurrences", "documents" };

    /// <summary>
    /// <c>sdrf lint</c> — the concepts a set of SDRF documents wrote inconsistently, from mzLib's
    /// <see cref="SdrfDriftLint"/>. Documents arrive on stdin exactly as for <c>sdrf pool</c>:
    /// <c>path[\tlabel]</c>, one per line.
    /// </summary>
    /// <remarks>
    /// Long format: one row per finding × variant, so a finding with three spellings is three rows
    /// sharing a <c>finding_index</c>, the majority spelling first (<c>variant_rank</c> 0). The
    /// documents that used a spelling cross as a list of labels in one cell; a label is
    /// caller-chosen text, so joining them on any delimiter could not be undone.
    /// </remarks>
    public static object Lint(Program.Arguments arguments)
    {
        (List<string> paths, List<string> labels) = ReadLabelledDocuments();

        var collection = new SdrfCollection(
            paths.Select(p => new SdrfDocument(p)),
            labels.Count == 0 ? null : labels);
        IReadOnlyList<SdrfDriftFinding> findings = SdrfDriftLint.Analyze(collection);

        var table = new ColumnTable(LintColumns);
        for (int f = 0; f < findings.Count; f++)
        {
            SdrfDriftFinding finding = findings[f];
            for (int v = 0; v < finding.Variants.Count; v++)
            {
                SdrfDriftVariant variant = finding.Variants[v];
                table.Add(f, finding.Kind.ToString(), finding.Concept, finding.Column, v,
                    variant.Value, variant.Occurrences, variant.Documents.ToList());
            }
        }

        return new
        {
            document_count = collection.Count,
            paths,
            labels = collection.Labels.ToList(),
            finding_count = findings.Count,
            column_names = table.Names,
            columns = table.Columns(),
            caveats = LintCaveats(labelsWereSupplied: labels.Count != 0),
        };
    }

    private static List<string> LintCaveats(bool labelsWereSupplied)
    {
        var caveats = new List<string>
        {
            "variant_rank 0 is the MAJORITY spelling. That describes what these documents do; it is " +
            "not advice. In the curated corpus characteristics[organism] is free text in 963 " +
            "documents and a CV term in 275, so following the majority would throw the accession away.",

            "Columns whose values are unique per row by construction are never compared: source name, " +
            "assay name, comment[data file], comment[searched data file] (since mzLib 1.0.592, #1335), " +
            $"comment[file uri] and {SdrfCollection.SourceDocumentColumn}. File names that differ " +
            "between documents are data, not drift.",

            "Reserved words and empty cells are skipped: a reserved word cannot drift. A set of " +
            "documents that say nothing therefore lints clean; sdrf assess is the check for that.",

            "Findings are ordered by impact (occurrences outside the majority), then kind, column and " +
            "concept, so the same documents always produce the same table and two runs can be diffed.",

            "column_name is null for a ColumnNameVariant finding, which is about a column's NAME; there " +
            "concept is the name lower-cased with runs of whitespace collapsed.",
        };

        if (!labelsWereSupplied)
        {
            caveats.Add(
                "No labels were supplied, so documents are named by mzLib's default of " +
                "containing-folder/file-stem - which depends on WHERE THE FILES SIT. Supply a label " +
                "per document for anything you intend to keep or compare across machines.");
        }

        return caveats;
    }
}
