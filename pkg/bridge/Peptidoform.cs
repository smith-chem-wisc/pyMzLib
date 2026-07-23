using System.Net.Http;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace MzLibBridge;

/// <summary>
/// The peptidoform workflow: fetch an annotated protein, digest it, and fragment the peptides.
/// </summary>
/// <remarks>
/// <para>
/// This is deliberately ONE verb answering one question a mass spectrometrist actually asks —
/// "what fragments would I see for this protein's peptides?" — rather than three verbs the caller
/// must compose. The three mzLib stages involved (<see cref="ProteinDbLoader"/>,
/// <see cref="Protein.Digest"/>, <see cref="PeptideWithSetModifications.Fragment"/>) traffic in
/// rich objects that carry behaviour: a <see cref="PeptideWithSetModifications"/> knows how to
/// fragment itself. Serialising those objects out and back between steps would both cost more
/// than the work and lose the behaviour that makes them useful, so the composition happens here,
/// where the objects live.
/// </para>
/// <para>
/// The defaults are the lab's opinion — tryptic, two missed cleavages, ETD, both termini — so
/// that the common question needs no parameters at all. Every one of them is reachable, because
/// the point is to open the doors, not to hide them.
/// </para>
/// </remarks>
internal static class Peptidoform
{
    /// <summary>
    /// How a UniProtKB entry's XML is obtained, given an accession, returning a local file path.
    /// </summary>
    /// <remarks>
    /// Injectable for the same reason <see cref="Program.PrideClientFactory"/> is: the parts of
    /// this workflow worth testing — the annotation census, the digestion, the modification
    /// combinatorics, the isoform cap — are all downstream of the download, and none of them
    /// should require the network to exercise.
    /// </remarks>
    /// <remarks>
    /// Returns the XML path and whether the caller is responsible for deleting it. The real
    /// downloader creates a temp file and hands ownership over (delete it); a test fixture is
    /// owned by the test and must be left alone. Making ownership explicit is what stops the
    /// cleanup from deleting a file it did not create.
    /// </remarks>
    internal static Func<string, Task<(string Path, bool CallerDeletes)>> UniProtXmlSource { get; set; }
        = async accession => (await DownloadUniProtXmlAsync(accession).ConfigureAwait(false), true);

    /// <summary>Where a single UniProtKB entry's annotated XML comes from.</summary>
    private const string UniProtEntryUrlFormat = "https://rest.uniprot.org/uniprotkb/{0}.xml";

    /// <summary>
    /// UniProt's modification definitions, resolved once per process.
    /// </summary>
    /// <remarks>
    /// A UniProt entry's XML <i>names</i> its modifications ("Phosphoserine", "N-linked
    /// (GlcNAc) asparagine") but does not define them — no formula, no mass. Without the
    /// definitions from UniProt's own <c>ptmlist.txt</c>, <see cref="ProteinDbLoader"/> parses the
    /// entry perfectly, resolves nothing, and returns a protein with zero modification sites.
    /// <para>
    /// That failure is completely silent, which is the point worth remembering: the first version
    /// of this code omitted the list and reported 195 peptides and 0 modifications for serum
    /// albumin, a protein with dozens of annotated sites. Nothing errored. The only reason it was
    /// caught is that the count was printed and looked implausible.
    /// </para>
    /// </remarks>
    private static readonly Lazy<List<Modification>> KnownUniProtModifications = new(LoadUniProtModifications);

    /// <summary>Reads UniProt's PTM list, with formal charges from the PSI-MOD ontology.</summary>
    private static List<Modification> LoadUniProtModifications()
    {
        string resources = Path.Combine(AppContext.BaseDirectory, "Resources");
        string ptmList = Path.Combine(resources, "ptmlist.txt");
        string psiModObo = Path.Combine(resources, "PSI-MOD.obo.xml");

        if (!File.Exists(ptmList))
            throw new FileNotFoundException(
                $"UniProt's modification list is missing from the payload at '{ptmList}'. Without it, " +
                "annotated modifications resolve to nothing and every protein looks unmodified.",
                ptmList);

        // PSI-MOD carries the formal charges that make a permanently charged modification's mass
        // encode its charge - which is exactly what mz() reads back to avoid double-counting a
        // proton (see FormalChargeOf). Substituting an empty dictionary when the file is absent
        // was a silent hole: trimethyllysine and the like would load without their fixed charge,
        // and every downstream m/z would be wrong by a proton with nothing to indicate it. Fail as
        // loudly here as for a missing ptmlist, since correctness now depends on it.
        if (!File.Exists(psiModObo))
            throw new FileNotFoundException(
                $"PSI-MOD is missing from the payload at '{psiModObo}'. It supplies the formal " +
                "charges that let m/z avoid double-counting a proton on charged modifications; " +
                "without it those masses load uncharged and every m/z is silently wrong.",
                psiModObo);

        Dictionary<string, int> formalCharges =
            Loaders.GetFormalChargesDictionary(Loaders.ReadPsiModFile(psiModObo));

        return Loaders.LoadUniprot(ptmList, formalCharges).ToList();
    }

    /// <summary>
    /// <c>peptidoform fragments --accession P02768 [--protease trypsin] [--dissociation ETD]
    /// [--no-modifications] [--missed-cleavages 2] [--min-length 7] [--max-length 0]
    /// [--max-mods 2] [--terminus Both]</c>
    /// </summary>
    public static async Task<object> FragmentsAsync(Program.Arguments arguments)
    {
        string accession = arguments.Required("accession");

        // Default to the classic Keil rule — cleave after K/R except before proline — because that
        // is what a mass spectrometrist means by "trypsin".
        //
        // The naming here is a genuine trap and is worth stating plainly. In mzLib, "trypsin" is
        // K|,R| and cleaves before proline, while "trypsin|P" is K[P]|,R[P]| and does not: the
        // bracket denotes a PREVENTING residue. That is the opposite of the MaxQuant and Mascot
        // convention, where "Trypsin/P" means the proline rule is *ignored*. Someone reaching for
        // the familiar-looking name gets the opposite of their intent, in either direction, with
        // no indication anything happened.
        string protease = arguments.Optional("protease") ?? "trypsin|P";
        string dissociationName = arguments.Optional("dissociation") ?? "ETD";
        string terminusName = arguments.Optional("terminus") ?? "Both";
        bool applyModifications = !arguments.Flag("no-modifications");
        int missedCleavages = arguments.OptionalInt("missed-cleavages", 2);
        int minLength = arguments.OptionalInt("min-length", 7);
        int maxLength = arguments.OptionalInt("max-length", 0);
        int maxMods = arguments.OptionalInt("max-mods", 2);

        // DigestionParams caps modification isoforms per peptide at 1024 by default and truncates
        // silently when it binds. On a histone at three or more modifications that ceiling is
        // reachable, and the caller would receive a short list with no indication it was cut.
        // Exposed so it can be raised, and reported below so it can at least be noticed.
        int maxIsoforms = arguments.OptionalInt("max-isoforms", 1024);

        if (!ProteaseDictionary.Dictionary.ContainsKey(protease))
            throw new Program.UsageException(
                $"Unknown protease '{protease}'. Known proteases include: " +
                string.Join(", ", ProteaseDictionary.Dictionary.Keys.Take(12)) + ", …");

        if (!Enum.TryParse(dissociationName, ignoreCase: true, out DissociationType dissociation))
            throw new Program.UsageException(
                $"Unknown dissociation type '{dissociationName}'. Known types: " +
                string.Join(", ", Enum.GetNames<DissociationType>()) + ".");

        if (!Enum.TryParse(terminusName, ignoreCase: true, out FragmentationTerminus terminus))
            throw new Program.UsageException(
                $"Unknown terminus '{terminusName}'. Known values: " +
                string.Join(", ", Enum.GetNames<FragmentationTerminus>()) + ".");

        (Protein annotated, Dictionary<string, int> annotationCensus, List<string> unresolved) =
            await FetchAnnotatedProteinAsync(accession).ConfigureAwait(false);

        // "Without modifications" is the same sequence with the annotations discarded, not a
        // different protein. Keeping the accession and sequence identical is what makes the two
        // runs comparable — the only variable is whether UniProt's annotations were applied.
        Protein subject = applyModifications
            ? annotated
            : new Protein(annotated.BaseSequence, annotated.Accession, annotated.Organism);

        var digestionParams = new DigestionParams(
            protease: protease,
            maxMissedCleavages: missedCleavages,
            minPeptideLength: minLength,
            maxPeptideLength: maxLength > 0 ? maxLength : int.MaxValue,
            maxModsForPeptides: maxMods,
            maxModificationIsoforms: maxIsoforms,
            fragmentationTerminus: terminus);

        // The raw digest, BEFORE de-duplication. The isoform-cap check below must count against
        // this, not the deduped list: mzLib truncates a locus at maxIsoforms *distinct* forms,
        // and de-duplicating first can drop a truncated locus below the cap and report it as
        // untruncated — the exact silent-truncation failure the cap reporting exists to prevent.
        List<PeptideWithSetModifications> rawDigest = subject
            .Digest(digestionParams, new List<Modification>(), new List<Modification>())
            .ToList();

        int peptidesAtCap = rawDigest
            .GroupBy(p => (p.OneBasedStartResidue, p.OneBasedEndResidue))
            .Count(g => g.Count() >= maxIsoforms);

        // WORKAROUND for smith-chem-wisc/mzLib#1108: Digest returns the same peptidoform twice
        // where a chain annotation's start coincides with the initiator-methionine cleavage site
        // (the mature N-terminus). On histone H3.1 that is 19 exact duplicates - identical
        // sequence, position, mods and mass, differing only in PeptideDescription. It is invisible
        // in a search, where PSMs are de-duplicated downstream, but a verb that ENUMERATES
        // peptidoforms would report a count ~3.5% high. Group to distinct peptidoforms until the
        // duplicate emission is fixed upstream; then this collapses to a plain ToList().
        List<PeptideWithSetModifications> peptides = rawDigest
            .GroupBy(p => (p.OneBasedStartResidue, p.OneBasedEndResidue, p.FullSequence))
            .Select(g => g.First())
            .ToList();

        var products = new List<Product>();
        var wirePeptides = new List<object>(peptides.Count);
        foreach (PeptideWithSetModifications peptide in peptides)
        {
            products.Clear();
            peptide.Fragment(dissociation, terminus, products);
            wirePeptides.Add(ToWirePeptide(peptide, products));
        }

        return new
        {
            accession = annotated.Accession,
            name = annotated.Name,
            full_name = annotated.FullName,
            organism = annotated.Organism,
            sequence_length = annotated.BaseSequence.Length,
            modifications_applied = applyModifications,
            // Two different numbers, and conflating them is exactly the sort of quiet error this
            // whole verb exists to avoid. OneBasedPossibleLocalizedModifications is keyed by
            // POSITION, and a histone carries several alternatives at one residue - K9me1, K9me2,
            // K9me3, K9ac are four modifications at one site. Reporting the dictionary count as
            // "modifications" made H3.1 look as though 93 annotations had been dropped when they
            // had all been loaded and merely shared residues.
            annotated_modification_sites = annotated.OneBasedPossibleLocalizedModifications.Count,
            annotated_modifications_loaded = annotated.OneBasedPossibleLocalizedModifications
                .SelectMany(kv => kv.Value).Count(),
            uniprot_annotated_features = annotationCensus.Sum(kv => kv.Value),
            // The exclusions LoadProteinXML could not resolve individually - a name absent from
            // UniProt's ptmlist, for instance. Reporting only the excluded *types* accounted for
            // 3 of 10 exclusions on histone H3.1 and said nothing about the other 7, which is
            // worse than saying nothing at all: it creates the feeling of having been told.
            unresolved_modifications = unresolved,
            uniprot_features_by_type = annotationCensus.OrderBy(kv => kv.Key)
                .Select(kv => new { type = kv.Key, count = kv.Value, loaded = kv.Key is "modified residue" or "lipid moiety-binding region" })
                .ToList(),
            protease,
            dissociation = dissociation.ToString(),
            terminus = terminus.ToString(),
            max_modifications = maxMods,
            max_modification_isoforms = maxIsoforms,
            peptides_at_isoform_cap = peptides
                .GroupBy(p => (p.OneBasedStartResidue, p.OneBasedEndResidue))
                .Count(g => g.Count() >= maxIsoforms),
            peptide_count = wirePeptides.Count,
            peptides = wirePeptides,
        };
    }

    /// <summary>
    /// Downloads one UniProtKB entry's XML and loads it with its annotated modifications.
    /// </summary>
    /// <remarks>
    /// The XML carries its own PTM list, which <see cref="ProteinDbLoader"/> reads, so no external
    /// modification database is needed for annotations UniProt already states.
    /// </remarks>
    private static async Task<(Protein Protein, Dictionary<string, int> AnnotationCensus, List<string> Unresolved)> FetchAnnotatedProteinAsync(string accession)
    {
        // Ownership is explicit: only delete the file if the source handed it over. An injected
        // fixture is owned by the test and left alone; the real downloader's temp is ours to clean.
        (string xmlPath, bool callerDeletes) = await UniProtXmlSource(accession).ConfigureAwait(false);

        try
        {
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: KnownUniProtModifications.Value,
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out Dictionary<string, Modification> unresolved);

            if (proteins.Count == 0)
                throw new Program.UsageException(
                    $"UniProt returned an entry for '{accession}' that contained no protein sequence.");

            return (proteins[0], CensusAnnotatedFeatures(xmlPath), unresolved.Keys.OrderBy(k => k).ToList());
        }
        finally
        {
            if (callerDeletes && File.Exists(xmlPath))
                File.Delete(xmlPath);
        }
    }

    /// <summary>
    /// Turns a UniProt response status into the right kind of failure, or returns quietly.
    /// </summary>
    /// <remarks>
    /// The distinction is the same one the whole bridge rests on, applied to a second service.
    /// UniProt answers a malformed accession with <b>400</b> and an unknown-but-well-formed one
    /// with <b>404</b>: both are the caller's mistake and permanent, so they must surface as usage
    /// errors. Telling someone to retry later for a typo would mean retrying forever, and would
    /// make a genuine bad-accession regression skip in the live suites rather than fail.
    /// <para>
    /// Everything else non-success is phrased as "failed with status NNN" so that
    /// <see cref="Program.ClassifyError"/> can read the code and decide: 408/429/5xx are
    /// availability, the rest are correctness.
    /// </para>
    /// </remarks>
    internal static void ThrowIfUniProtRejected(
        System.Net.HttpStatusCode status, string? reasonPhrase, string accession, string url)
    {
        if (status is System.Net.HttpStatusCode.NotFound or System.Net.HttpStatusCode.BadRequest)
        {
            throw new Program.UsageException(
                $"UniProt has no entry '{accession}'. Accessions look like 'P02768' or 'A0A0B4J2D5'.");
        }

        if ((int)status is >= 200 and < 300)
            return;

        throw new HttpRequestException(
            $"UniProt request failed with status {(int)status} {reasonPhrase} for '{url}'.");
    }

    /// <summary>Downloads one UniProtKB entry's XML to a temporary file and returns its path.</summary>
    private static async Task<string> DownloadUniProtXmlAsync(string accession)
    {
        string url = string.Format(UniProtEntryUrlFormat, Uri.EscapeDataString(accession));
        string temp = Path.Combine(Path.GetTempPath(), $"pymzlib-{Guid.NewGuid():N}.xml");

        using var client = new HttpClient { Timeout = TimeSpan.FromSeconds(100) };
        using HttpResponseMessage response = await client.GetAsync(url).ConfigureAwait(false);

        ThrowIfUniProtRejected(response.StatusCode, response.ReasonPhrase, accession, url);

        await using Stream body = await response.Content.ReadAsStreamAsync().ConfigureAwait(false);
        await using var file = new FileStream(temp, FileMode.Create, FileAccess.Write, FileShare.None);
        await body.CopyToAsync(file).ConfigureAwait(false);
        return temp;
    }

    /// <summary>
    /// Counts the modification-like features UniProt annotates, by type, straight from the XML.
    /// </summary>
    /// <remarks>
    /// mzLib deliberately loads only modifications with a defined target and a defined mass —
    /// a glycosylation site annotated as "N-linked (GlcNAc...)" has no mass to add, and a peptide
    /// carrying an ambiguous-mass PTM is not identifiable by mass spectrometry anyway. That is the
    /// right call, and this does not second-guess it.
    /// <para>
    /// What it adds is the part that was missing: telling the caller. Without this, albumin
    /// reports 14 modification sites and there is no way to learn that UniProt annotates 38, that
    /// 24 were excluded, or why. The number was always correct under a rule the user could not
    /// see. Reporting the census makes the rule visible without changing it.
    /// </para>
    /// </remarks>
    private static Dictionary<string, int> CensusAnnotatedFeatures(string xmlPath)
    {
        var census = new Dictionary<string, int>(StringComparer.Ordinal);
        var document = new System.Xml.XmlDocument();
        document.Load(xmlPath);

        foreach (System.Xml.XmlNode node in document.GetElementsByTagName("feature"))
        {
            string? type = node.Attributes?["type"]?.Value;
            if (string.IsNullOrEmpty(type))
                continue;
            if (type is not ("modified residue" or "glycosylation site" or "lipid moiety-binding region" or "cross-link"))
                continue;
            census[type] = census.TryGetValue(type, out int n) ? n + 1 : 1;
        }

        return census;
    }

    /// <summary>The mass of one electron, in daltons.</summary>
    private const double ElectronMass = 0.00054857990;

    /// <summary>
    /// The formal charge a modification carries, recovered from its own recorded mass.
    /// </summary>
    /// <remarks>
    /// Some modifications produce a permanently charged residue. Trimethylation of a lysine
    /// ε-amine gives a quaternary ammonium — <c>NH2 → N(CH3)3⁺</c> — which adds C₃H₇ and removes
    /// an electron, and UniProt records the delta as 43.054227 rather than the neutral-formula
    /// 43.054775. (Unimod's 42.04695 is C₃H₆, the neutral convention search engines use because
    /// they assume every charge comes from an added proton.)
    /// <para>
    /// mzLib is right to follow UniProt here, but a consumer computing m/z must know: a mass that
    /// already carries a charge needs fewer protons added, not the same number. The difference
    /// between the chemical formula's mass and the recorded mass is exactly one electron per
    /// formal charge, so the charge is recoverable without a second lookup table.
    /// </para>
    /// </remarks>
    private static int FormalChargeOf(Modification modification)
    {
        if (modification.ChemicalFormula is null || modification.MonoisotopicMass is null)
            return 0;

        double deficit = modification.ChemicalFormula.MonoisotopicMass - modification.MonoisotopicMass.Value;
        int charge = (int)Math.Round(deficit / ElectronMass);

        // Only trust a clean whole number of electrons; anything else is rounding in the source
        // data rather than a formal charge, and guessing would be worse than reporting none.
        return Math.Abs(deficit - charge * ElectronMass) < 1e-5 ? charge : 0;
    }

    /// <summary>Flattens a digested peptide and its fragments into the wire shape.</summary>
    private static object ToWirePeptide(PeptideWithSetModifications peptide, List<Product> products) => new
    {
        base_sequence = peptide.BaseSequence,
        full_sequence = peptide.FullSequence,
        monoisotopic_mass = peptide.MonoisotopicMass,
        length = peptide.Length,
        one_based_start = peptide.OneBasedStartResidue,
        one_based_end = peptide.OneBasedEndResidue,
        missed_cleavages = peptide.MissedCleavages,
        modification_count = peptide.AllModsOneIsNterminus.Count,
        // Charges the INTACT peptide already carries before protonation, used by Peptide.mz to add
        // (z - fixed_charges) protons rather than z. This is a whole-peptide sum and deliberately
        // does NOT apply to fragments: a c or z ion carries only the fixed charges of the residues
        // within its own span, so fragment m/z would need per-fragment accounting. Fragments
        // therefore expose neutral_mass only, and the limitation is documented on the Python side.
        fixed_charges = peptide.AllModsOneIsNterminus.Values.Sum(FormalChargeOf),
        // AllModsOneIsNterminus is keyed with slot 1 as the N-TERMINUS, so residue i lives at
        // key i+1. Exposing that key as a "one-based position" was a lie about what the number
        // is: 474 of 498 modifications pointed one residue past their own target, and a peptide
        // MAR reported position 4 for its arginine — position 4 of a 3-mer.
        //
        // Terminal modifications are reported as such rather than squeezed into a residue index
        // they do not have.
        modifications = peptide.AllModsOneIsNterminus
            .Select(kv => new
            {
                one_based_residue = kv.Key is 1 ? (int?)null
                    : kv.Key >= peptide.Length + 2 ? (int?)null
                    : kv.Key - 1,
                terminus = kv.Key is 1 ? "N"
                    : kv.Key >= peptide.Length + 2 ? "C"
                    : null,
                id = kv.Value.IdWithMotif,
                mass = kv.Value.MonoisotopicMass,
                formal_charge = FormalChargeOf(kv.Value),
            })
            .ToList(),
        fragments = products.Select(p => new
        {
            product_type = p.ProductType.ToString(),
            fragment_number = p.FragmentNumber,
            neutral_mass = p.NeutralMass,
            neutral_loss = p.NeutralLoss,
            residue_position = p.ResiduePosition,
        }).ToList(),
    };
}
