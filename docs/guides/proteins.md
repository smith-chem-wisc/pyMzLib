# Protein databases

A search result tells you *which accession* a peptide matched. It does not tell you what organism
that accession is, what it does, which gene encodes it, or whether the peptide could have come from
somewhere else. **The protein database you searched already knows all four**, and
`pymzlib.proteins` reads them out of it with mzLib:

| Question | Function | mzLib |
|---|---|---|
| What organism, taxon, gene and mass is this accession? What are its GO terms? | [`read()`](#look-up-organism-and-taxon-for-a-list-of-accessions) | `ProteinDbLoader`, `Protein.GoTerms` ([#1336][1336]), `Protein.EnsemblGeneReferences` |
| Which Ensembl gene is this protein, in a way I can reproduce next year? | [`resolve_genes()`](#resolve-proteins-to-ensembl-genes-reproducibly) | `EnsemblGeneResolver` ([#1338][1338]) |
| Does this peptide identify one protein, one gene, or neither? | [`classify_peptides()`](#decide-whether-a-peptide-is-unique-with-i-l) | `PeptideUniquenessClassifier` ([#1348][1348]) |

```python
import pymzlib

db = pymzlib.proteins.read("uniprotkb_human_proteome.xml.gz")
print(db.record_count)                         # one row per protein
print(db.taxonomy()["P04406"])                 # '9606'
```

Every example on this page runs against the small databases committed in pyMzLib's test suite
(`pkg/python/tests/fixtures/proteins/`): two real UniProt XML entries (GAPDH, PSCA), human albumin
and three DYNC1I2 isoforms as FASTA, mouse AIFM1 as FASTA, and bovine albumin as a contaminant.
The outputs shown are what those files produce.

## Before you start: XML or FASTA?

**Use the UniProt XML if you have the choice.** Both formats give you accession, organism, gene
name and sequence. Only the XML carries cross-references, so only the XML has GO terms and Ensembl
genes:

| | UniProt XML (`.xml`, `.xml.gz`) | FASTA (`.fasta`, `.fa`, `.faa`, `.fas`, `.gz`) |
|---|---|---|
| organism | `<organism>` | `OS=` |
| NCBI taxonomy id | `<dbReference type="NCBI Taxonomy">` | `OX=` (UniProt headers; absent in older and non-UniProt FASTA) |
| gene names | primary, synonyms, ORF and locus names | `GN=` only |
| GO terms | ✅ | ❌ |
| Ensembl gene links | ✅ | ❌ |

A FASTA's empty GO table is **the format being silent**, not a protein with no annotation. pyMzLib
will not let those two look alike: every result lists each input in `files`, and a FASTA's entry
says so explicitly.

```python
db = pymzlib.proteins.read(["human_subset.xml", "human_extra.fasta"], tables=("proteins", "go_terms"))
for f in db.files:
    print(f.file_type, f.absent_fields)
# UniProtXml []
# Fasta ['go_terms', 'ensembl_genes', 'ensembl_gene_ids']
```

Download the XML from UniProt with `format=xml` on the same query you would use for FASTA, e.g.
`https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640&format=xml&compressed=true`.
It is larger (the human reference proteome is a few hundred MB compressed) and takes longer to read,
which is the price of the annotations.

## Look up organism and taxon for a list of accessions

**The problem.** MetaMorpheus writes the organism *name* in its protein-group table, not the
NCBI taxonomy id, and a mixed-species search (a host and a pathogen, a xenograft, a spike-in) needs
the id to join against anything taxonomic. You have the accessions; the database you searched has
the answers.

```python
import pymzlib

accessions = ["P04406", "Q9Z0X1", "P02769", "P04406-1"]
db = pymzlib.proteins.read(
    ["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"],
    contaminants=["contaminants.fasta"],
    accessions=accessions,
)

print(db.taxonomy())
# {'P04406': '9606', 'Q9Z0X1': '10090', 'P02769': '9913'}
print(db.organisms()["P02769"])
# 'Bos taurus'
print(db.accessions_not_found)
# ['P04406-1']
```

Three things in that output are worth reading twice:

- **`accessions_not_found` is not empty.** Matching is *exact*: `P04406-1` is not `P04406`. The
  filter never guesses which isoform you meant, and a miss is reported rather than dropped.
  Always check it.
- **The taxon is a string** (`'9606'`), because it is an identifier, not a quantity. `None` means
  the entry records none - a FASTA header without `OX=` - never a guessed id.
- **The bovine albumin came from `contaminants=`**, so its row has `is_contaminant == True`. Marking
  contaminants matters more for the other two functions than for this one.

The whole protein table is column-major, so it goes straight into pandas:

```python
import pandas as pd

proteins = pd.DataFrame(db.columns)
proteins[["accession", "organism", "ncbi_taxonomy_id", "primary_gene_name", "length", "monoisotopic_mass"]]
```

| accession | organism | ncbi_taxonomy_id | primary_gene_name | length | monoisotopic_mass |
|---|---|---|---|---|---|
| P04406 | Homo sapiens | 9606 | GAPDH | 335 | 36030.40 |
| Q9Z0X1 | Mus musculus | 10090 | Aifm1 | 612 | 66723.84 |
| P02769 | Bos taurus | 9913 | ALB | 607 | 69248.44 |

`monoisotopic_mass` is in daltons, for the **unmodified sequence as written** - initiator
methionine, signal peptide and propeptide included - so it is the precursor's mass, not the mature
protein's. It is `None` when the sequence holds a letter with no defined mass (`X`, `B`, `Z`, `J`).

**Many accessions, many databases.** The accessions travel on stdin, not the command line, so a list
of 20,000 is fine. Databases are read in the order given; `threads=4` reads four at once, and the
result is identical at any thread count - only memory changes, since each database is read whole.

## Get GO terms for a protein group

**The problem.** A protein group from your search is significant, and you want to know what its
members do - but only from annotations that were not assigned automatically.

GO terms come from the UniProt XML you already have, with no extra download: UniProt ships every
entry's GO annotations as `<dbReference type="GO">`, with the aspect as a prefix on the term
(`C:cytoplasm`) and one reference per line of evidence. mzLib (#1336) turns that into one term per
GO id with its evidence unioned, and pyMzLib returns it as a long table:

```python
group = ["O43653"]      # PSCA, from a protein group's accession list
db = pymzlib.proteins.read("human_subset.xml", accessions=group, tables=("go_terms",))
go = pd.DataFrame(db.go_terms.columns)
go[["go_id", "aspect", "term_name", "evidence_codes", "projects"]]
```

| go_id | aspect | term_name | evidence_codes | projects |
|---|---|---|---|---|
| GO:0031225 | CellularComponent | anchored component of membrane | [ECO:0000501] | [UniProtKB-KW] |
| GO:0070062 | CellularComponent | extracellular exosome | [ECO:0007005] | [UniProtKB] |
| GO:0005576 | CellularComponent | extracellular region | [ECO:0000304] | [Reactome] |
| GO:0016020 | CellularComponent | membrane | [ECO:0000318] | [GO_Central] |
| GO:0005886 | CellularComponent | plasma membrane | [ECO:0000314] | [HPA] |
| GO:0033130 | MolecularFunction | acetylcholine receptor binding | [ECO:0000314] | [UniProtKB] |
| GO:0070373 | BiologicalProcess | negative regulation of ERK1 and ERK2 cascade | [ECO:0000314] | [UniProtKB] |
| GO:0099601 | BiologicalProcess | regulation of neurotransmitter receptor activity | [ECO:0000314] | [UniProtKB] |

Evidence codes are [ECO](https://www.evidenceontology.org/) ids, sorted. To drop annotations that
rest only on automatic assertion (`ECO:0000501`), filter on them - it is a column, not a hidden mode:

```python
automatic = {"ECO:0000501"}
curated = go[go["evidence_codes"].map(lambda codes: bool(set(codes) - automatic))]
```

What the table does and does not say:

- **`aspect` is `"Unknown"` when UniProt gave no `C:`/`F:`/`P:` prefix.** It is never inferred:
  an organelle claim built on a guessed aspect is worse than one that admits it does not know.
- **One row per GO id per protein.** A term UniProt repeats with three lines of evidence is one row
  with three evidence codes, so counting rows counts terms.
- **These are the terms as annotated, not propagated.** A protein annotated to `plasma membrane` is
  not also listed under `membrane` unless UniProt says so; walking the GO graph is not done here.
  Enrichment tools expect to do that themselves.

For a whole group, pass all its accessions and group by `accession` - or by `go_id` to see which
terms the members share.

## Resolve proteins to Ensembl genes, reproducibly

**The problem.** You want gene-level results - to join with RNA-seq, to count genes, to compare
with last year's analysis - and "the gene" has to mean the same thing each time. A gene *symbol*
cannot do that: symbols are renamed, one symbol can name several loci, and the `Gene` column of a
protein-group table is only the first name per protein. A stable Ensembl gene id counted against
one Ensembl release can.

`resolve_genes()` reads the Ensembl links the UniProt XML carries for each protein and counts them
against **a gene set you supply**: an Ensembl GTF for one release.

```python
genes = pymzlib.proteins.resolve_genes(
    ["human_subset.xml", "human_extra.fasta"],
    contaminants=["contaminants.fasta"],
    gtf="Homo_sapiens.GRCh38.116.gtf",                 # yours: the release's .gtf.gz
    xref="Homo_sapiens.GRCh38.116.uniprot.tsv",        # optional second opinion
)
print(genes.gene_set.release, genes.gene_set.genome_build, genes.gene_set.sha256[:12])
# 116 GRCh38.p14 cfe0494115dd
print(genes.outcome_counts)
# {'resolved': 1, 'multi_gene': 0, 'off_primary_only': 1, 'not_in_source': 4,
#  'unrecognized_accession': 0, 'contaminant_not_mapped': 1}
```

| accession | outcome | n_genes | gene_id | versioned_gene_id | gene_symbol | off_primary_genes | ensembl_xref_agrees |
|---|---|---|---|---|---|---|---|
| P04406 | resolved | 1 | ENSG00000111640 | ENSG00000111640.15 | GAPDH | 0 | True |
| O43653 | off_primary_only | 0 | | | | 1 | |
| P02768 | not_in_source | 0 | | | | 0 | |
| Q13409 | not_in_source | 0 | | | | 0 | |
| P02769 | contaminant_not_mapped | 0 | | | | 0 | |

(The example's gene set is a three-gene test GTF, which is why PSCA's gene is "off the assembly".)

### Every protein gets exactly one outcome

`None` in `gene_id` can mean five different things, so the table never leaves it at `None`:

| outcome | what happened | what to do |
|---|---|---|
| `resolved` | exactly one gene on the gene set's assembly | use `gene_id` |
| `multi_gene` | several genes on the assembly: **one row per gene**, never a pick | decide per question; `n_genes` says how many |
| `off_primary_only` | the database links only genes outside the gene set (ALT haplotypes, patches) | the source has an answer, just not on the assembly everything else is counted against |
| `not_in_source` | a well-formed accession the database links to no gene | for a FASTA, always - see below |
| `unrecognized_accession` | not a UniProt or RefSeq accession (a decoy or entrapment prefix, a mangled id) | a data problem, not a biology result |
| `contaminant_not_mapped` | from a `contaminants=` database: never mapped, never dropped | nothing - that is the point |

### Why you supply the GTF, and why you pin it

Nothing is downloaded or defaulted. **Which release you resolve against is part of the answer**:
between releases, genes are merged, split and retired, versions increment, symbols change, and the
set of genes on the primary assembly moves. A resolution with no release attached cannot be
reproduced or compared.

- **Use the primary-assembly GTF**, `Homo_sapiens.GRCh38.<release>.gtf.gz` from
  `https://ftp.ensembl.org/pub/release-<release>/gtf/homo_sapiens/` - **not** the
  `chr_patch_hapl_scaff` file. UniProt links a protein to every Ensembl transcript that encodes it,
  alternative haplotypes included, where one locus is described many times. Counting those as
  separate genes made human release 116 look 6.99% multi-gene when it is 0.36%.
- **Keep Ensembl's file name.** The release is read from it (`.116.gtf.gz` → `"116"`); a renamed file
  gives `gene_set_release == None` and a caveat.
- **Record the hashes.** Every row carries `search_database_sha256` (of the *decompressed* database,
  so `.xml` and `.xml.gz` agree), `gene_set_sha256`, and `ensembl_xref_sha256`. Together they name
  exactly the inputs a row was computed from.
- **A GTF is large** (human 116: 141 MB compressed, 4.66 GB unzipped) and only its gene rows are
  read. If you resolve often, write the gene rows once with mzLib's `EnsemblGeneSetWriter` (about
  0.5 MB) and pass that as `gene_set=` instead. It carries the GTF's provenance, so rows are keyed
  exactly as against the GTF.

### The optional Ensembl cross-reference

UniProt and Ensembl each map accessions to genes, and they are expected to disagree for real
reasons (readthrough genes, identical paralogs). Passing Ensembl's own
`Homo_sapiens.GRCh38.<release>.uniprot.tsv.gz` as `xref=` keeps both answers:
`ensembl_xref_agrees` says per row whether Ensembl links the same pair, `ensembl_xref_info_type`
gives Ensembl's evidence (`DIRECT`, `SEQUENCE_MATCH`, `INFERRED_PAIR`), and a gene only Ensembl
links gets a row of its own with `source == "ensembl_xref"`. Those rows never change a protein's
outcome. Without `xref=`, agreement is `None` - **unknown, not false**.

### A FASTA cannot resolve

A FASTA header has no Ensembl links, so **every FASTA protein is `not_in_source`**, whatever Ensembl
knows about it. The result's `caveats` counts them. Resolve against the UniProt XML of the same
proteome - the accessions are the same, so you can resolve a FASTA search's results that way.

## Decide whether a peptide is unique, with I = L

**The problem.** Protein inference, isoform claims and "is this contamination?" all start from the
same question: *could this peptide have come from anything else in the search space?*

```python
calls = pymzlib.proteins.classify_peptides(
    ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
    ["human_subset.xml", "human_extra.fasta"],
    contaminants=["contaminants.fasta"],
)
pd.DataFrame(calls.columns)
```

| peptide | sharing | accessions | shared_gene_keys |
|---|---|---|---|
| VGVNGFGR | Unique | [P04406] | [ensembl:ENSG00000111640, entry:P04406, gene:Homo sapiens:GAPDH] |
| LVLNGNPLTLFQER | Unique | [P04406] | [ensembl:ENSG00000111640, entry:P04406, gene:Homo sapiens:GAPDH] |
| ALSEQINIFFDYSGR | SharedWithinGene | [Q13409, Q13409-2, Q13409-3] | [entry:Q13409, gene:Homo sapiens:DYNC1I2] |
| YLYEIAR | SharedAcrossGenes | [P02768, P02769] | [] |
| AEFVEVTK | Unique | [P02769] | [entry:P02769, gene:Bos taurus:ALB] |
| PEPTIDEK | NotInDatabase | [] | [] |

Reading it row by row:

- **`LVLNGNPLTLFQER` is GAPDH's `LVINGNPITIFQER`.** I and L have the same mass, so a mass
  spectrometer cannot tell them apart and neither does the classifier: every `I` is matched as `L`.
  The `peptide` column echoes what you passed, unfolded. `calls.i_and_l_equivalent` is always
  `True`, so the rule travels with the result.
- **`ALSEQINIFFDYSGR` is in all three DYNC1I2 isoforms.** It supports the gene, not any one isoform.
  An isoform-level claim needs a `Unique` peptide.
- **`YLYEIAR` is in human *and* bovine albumin.** It is evidence for albumin, but not for *human*
  albumin: on its own it cannot distinguish a sample protein from BSA contamination.
  `AEFVEVTK`, `Unique` to bovine albumin, points at the contaminant. This is why `contaminants=`
  matters: leave the contaminant database out and `YLYEIAR` would be called `Unique` to human
  albumin.

### The rules, exactly

- **Containment, not digestion.** A protein contains a peptide if its sequence contains it
  *anywhere*, whatever the protease. That is deliberately conservative: a peptide is `Unique` only
  if no other sequence contains it at all, so an isoform claim cannot be undone by an entry at a
  site your search's cleavage rules happened to skip.
- **Identical sequences are one sequence.** Two accessions with the same sequence cannot be told
  apart by any peptide, so a peptide in both is `Unique` and lists both accessions.
- **"Same gene" means sharing a gene key**: an Ensembl gene id (`ensembl:`, from the XML),
  organism plus primary gene name (`gene:Homo sapiens:GAPDH` - so human and bovine ALB stay apart),
  or the UniProt entry (`entry:`, so `P12345-2` meets `P12345`). A FASTA isoform and its XML entry
  meet on the entry key.
- **Decoys are ignored** (`decoy_proteins_ignored` counts them); **contaminants are included**.
- **Base sequences only, upper case.** Strip modifications first: mzLib refuses `PEPT[Phospho]IDEK`
  and `peptidek` with a `UsageError` naming the peptide, rather than guess what you meant.

Every database is part of one search space, so `classify_peptides` has no `on_error="skip"`: if a
database fails to load, the call fails, because carrying on would report the peptides only that
database contains as unique.

## Many databases, threads and failures

All three functions take one database or a list, and `contaminants=` as a separate list. They are
read in one bridge process, in order: `databases` first, then `contaminants`, and every table's
`source_index` column points into `files` in that order.

| Argument | Default | Meaning |
|---|---|---|
| `threads` | `1` | Databases read at once; `-1` means every core. The result is **byte-identical** at any value - this is a memory choice, since each database is read whole. |
| `on_error` | `"fail"` | `"fail"` raises on the first unreadable database in input order; `"skip"` records it in `files[i].error` and reads the rest (`read` and `resolve_genes` only). |

mzLib's loader reads each database exactly as a search would, except that **no decoys are
generated**. A decoy already in the file (an accession starting `DECOY`) is kept, with
`is_decoy == True`. A UniProt XML that records a **genotype** (a Spritz-style variant database) still
has those variants applied, which renames the accession - `P38936_C117Y` - and changes its length
and mass; that file's `caveats` say how many.

## Errors

| Raised | When |
|---|---|
| `UsageError` | no database; a blank or repeated path; a file that is not `.xml` or a FASTA extension; a missing file (under `on_error="fail"`); an unknown table; `threads` of `0` or below `-1`; an empty `accessions` list; neither or both of `gtf` / `gene_set`; a missing GTF or xref; a peptide that is not upper-case A-Z |
| `BridgeError` | mzLib could not parse a database, a GTF gene row has no `gene_id`, or an xref table's columns are not Ensembl's layout |

## References

- The Gene Ontology Consortium: Ashburner, M. *et al.* Gene Ontology: tool for the unification of
  biology. *Nature Genetics* **25**, 25-29 (2000).
  [doi:10.1038/75556](https://doi.org/10.1038/75556)
- The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids
  Research* **53**, D609-D617 (2025). [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
- Yates, A. D. *et al.* Ensembl 2026. *Nucleic Acids Research* **54**, D1053-D1060 (2026).
  [doi:10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239)

The same three verbs are specified once, language-neutrally, for pyMzLib, mzLibRust and mzLibR:
`proteins read`, `genes resolve` and `proteins classify-peptides`. See the
[API reference](../reference.md#pymzlibproteins) for every field.

[1336]: https://github.com/smith-chem-wisc/mzLib/pull/1336
[1338]: https://github.com/smith-chem-wisc/mzLib/pull/1338
[1348]: https://github.com/smith-chem-wisc/mzLib/pull/1348
