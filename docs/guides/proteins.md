# Protein databases

A search result tells you *which accession* a peptide matched. It does not tell you what organism
that accession is, what it does, which gene encodes it, or whether the peptide could have come from
somewhere else. **The protein database you searched already knows all four**, and
`pymzlib.proteins` reads them out of it with mzLib:

| You want to know | Call | mzLib |
|---|---|---|
| What organism, taxon, gene and mass is this accession? What are its GO terms? | [`read()`](#look-up-organism-and-taxon-for-a-list-of-accessions) | `ProteinDbLoader`, `Protein.GoTerms` ([#1336][1336]), `Protein.EnsemblGeneReferences` |
| Which Ensembl gene is this protein, in a way I can reproduce next year? | [`resolve_genes()`](#resolve-proteins-to-ensembl-genes-reproducibly) | `EnsemblGeneResolver` ([#1338][1338]) |
| Does this peptide identify one protein, one gene, or neither? | [`classify_peptides()`](#decide-whether-a-peptide-is-unique-with-i-l) | `PeptideUniquenessClassifier` ([#1348][1348]) |
| Which GO terms does each protein group carry, with every member kept and the ontology's ancestors filled in? | [`annotate_go()`](#annotate-protein-groups-with-go-keeping-every-member) | `GoGroupAnnotator` ([#1353][1353]), `ToGoAnnotationGroups` ([#1366][1366]) |
| Where do I get a go.obo, and how do I keep the release pinned? | [`update_go()`](#first-get-a-goobo-once-and-keep-it) | `Loaders.UpdateGeneOntology` ([#1353][1353]) |

```pycon
>>> import pymzlib
>>> databases = ["human_subset.xml", "human_extra.fasta", "mouse_aifm1.fasta"]
>>> db = pymzlib.proteins.read(databases, contaminants=["contaminants.fasta"],
...                            tables=pymzlib.proteins.TABLES)
>>> db.record_count                            # one row per protein
8
>>> db.taxonomy()["P04406"]
'9606'

```

Every `>>>` example on this page runs in CI against output recorded from the real bridge, on the
small databases committed in pyMzLib's test suite (`pkg/python/tests/fixtures/proteins/`): two real
UniProt XML entries (GAPDH, PSCA), human albumin and three DYNC1I2 isoforms as FASTA, mouse AIFM1
as FASTA, and bovine albumin as a contaminant. The Gene Ontology section adds a real MetaMorpheus
protein-group table and the UniProt entries it names. Blocks titled **Not run** say why.

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

```pycon
>>> for f in db.files:
...     print(f.file_type, f.contaminant, f.absent_fields)
UniProtXml False []
Fasta False ['go_terms', 'ensembl_genes', 'ensembl_gene_ids']
Fasta False ['go_terms', 'ensembl_genes', 'ensembl_gene_ids']
Fasta True ['go_terms', 'ensembl_genes', 'ensembl_gene_ids']

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

```pycon
>>> found = pymzlib.proteins.read(
...     databases,
...     contaminants=["contaminants.fasta"],
...     accessions=["P04406", "Q9Z0X1", "P02769", "P04406-1"],
... )
>>> found.taxonomy()
{'P04406': '9606', 'Q9Z0X1': '10090', 'P02769': '9913'}
>>> found.organisms()["P02769"]
'Bos taurus'
>>> found.accessions_not_found
['P04406-1']

```

Three things in that output are worth reading twice:

- **`accessions_not_found` is not empty.** Matching is *exact*: `P04406-1` is not `P04406`. The
  filter never guesses which isoform you meant, and a miss is reported rather than dropped.
  Always check it.
- **The taxon is a string** (`'9606'`), because it is an identifier, not a quantity. `None` means
  the entry records none - a FASTA header without `OX=` - never a guessed id.
- **The bovine albumin came from `contaminants=`**, so its row has `is_contaminant == True`. Marking
  contaminants matters more for the other two functions than for this one.

The whole protein table is column-major, one list per column:

```pycon
>>> for row in found.records:
...     print(row["accession"], row["organism"], row["primary_gene_name"], row["length"],
...           round(row["monoisotopic_mass"], 2), row["is_contaminant"])
P04406 Homo sapiens GAPDH 335 36030.4 False
Q9Z0X1 Mus musculus Aifm1 612 66723.84 False
P02769 Bos taurus ALB 607 69248.44 True

```

so it goes straight into pandas:

```python title="Not run: needs pandas, which pyMzLib does not depend on"
import pandas as pd
proteins = pd.DataFrame(found.columns)
```

`monoisotopic_mass` is in daltons, for the **unmodified sequence as written** - initiator
methionine, signal peptide and propeptide included - so it is the precursor's mass, not the mature
protein's. It is `None` when the sequence holds a letter with no defined mass (`X`, `B`, `Z`, `J`).

**Many accessions, many databases.** The accessions travel on stdin, not the command line, so a list
of 20,000 is fine. Databases are read in the order given; `threads=4` reads four at once, and the
result is identical at any thread count - only memory changes, since each database is read whole.

## Look up the GO terms a protein was annotated with

**The problem.** You want to know what a protein does, from the annotations UniProt itself records
for it - but only those that were not assigned automatically.

GO terms come from the UniProt XML you already have, with no extra download: UniProt ships every
entry's GO annotations as `<dbReference type="GO">`, with the aspect as a prefix on the term
(`C:cytoplasm`) and one reference per line of evidence. mzLib (#1336) turns that into one term per
GO id with its evidence unioned, and pyMzLib returns it as a long table:

Here are PSCA's (O43653), from the read above:

```pycon
>>> psca = [r for r in db.go_terms.records if r["accession"] == "O43653"]
>>> for r in psca:
...     print(r["go_id"], r["aspect"], r["term_name"], r["evidence_codes"], r["projects"])
GO:0031225 CellularComponent anchored component of membrane ['ECO:0000501'] ['UniProtKB-KW']
GO:0070062 CellularComponent extracellular exosome ['ECO:0007005'] ['UniProtKB']
GO:0005576 CellularComponent extracellular region ['ECO:0000304'] ['Reactome']
GO:0016020 CellularComponent membrane ['ECO:0000318'] ['GO_Central']
GO:0005886 CellularComponent plasma membrane ['ECO:0000314'] ['HPA']
GO:0033130 MolecularFunction acetylcholine receptor binding ['ECO:0000314'] ['UniProtKB']
GO:0070373 BiologicalProcess negative regulation of ERK1 and ERK2 cascade ['ECO:0000314'] ['UniProtKB']
GO:0099601 BiologicalProcess regulation of neurotransmitter receptor activity ['ECO:0000314'] ['UniProtKB']

```

Evidence codes are [ECO](https://www.evidenceontology.org/) ids, sorted. To drop annotations that
rest only on automatic assertion (`ECO:0000501`), filter on them; it is a column, not a hidden mode:

```pycon
>>> automatic = {"ECO:0000501"}
>>> curated = [r for r in psca if set(r["evidence_codes"]) - automatic]
>>> len(psca), len(curated)
(8, 7)

```

What the table does and does not say:

- **`aspect` is `"Unknown"` when UniProt gave no `C:`/`F:`/`P:` prefix.** It is never inferred:
  an organelle claim built on a guessed aspect is worse than one that admits it does not know.
- **One row per GO id per protein.** A term UniProt repeats with three lines of evidence is one row
  with three evidence codes, so counting rows counts terms.
- **These are the terms as annotated, not propagated.** A protein annotated to `plasma membrane` is
  not also listed under `membrane` unless UniProt says so. To put a term's ancestors in as well -
  and to do it for whole protein groups - use [`annotate_go()`](#annotate-protein-groups-with-go-keeping-every-member),
  next.

## Annotate protein groups with GO, keeping every member

**The problem.** Your search reports *protein groups*, and a group can have several members: two
histone H2A.Z variants that no peptide could tell apart are one group, `P0C0S5|Q71UI9`. A GO
enrichment tool wants one set of terms per group. The usual shortcut - take the first accession -
throws away half the evidence, and the first accession is not even a choice MetaMorpheus made: it
sorts members alphabetically and never picks a leading protein. And the terms UniProt records are
the most specific ones, so a group annotated to *mitochondrial inner membrane* is not, as written,
"in the mitochondrion" at all.

`annotate_go()` answers both, with mzLib's
`GoGroupAnnotator` ([#1353][1353]) reading the protein-group table MetaMorpheus already wrote
([#1366][1366]):

- **one row per (group, term)** that *any* member holds - directly, or through an ancestor in the
  ontology (`is_a` and `part_of`, nothing weaker);
- **every row says who carries the term**, how directly, and on what evidence, so the union, the
  consensus and the direct-only views are all filters you apply, not choices made for you;
- **every non-decoy group gets at least one row**, and a group with no term gets one row that says
  why.

| Question | Code | mzLib |
|---|---|---|
| Which GO terms does each group carry, members kept apart? | `annotate_go(groups, database, go_obo=...)` | `GoGroupAnnotator` ([#1353][1353]), `ToGoAnnotationGroups` ([#1366][1366]) |
| Which terms do *all* members share? | rows with `n_with == n_members` | the same rows |
| Which organelles, complexes, pathways - in *my* vocabulary? | `category_map=` | `GoCategoryResolver` |
| Which GO release did I use, exactly? | `go.go.sha256`, `go.header` | `GeneOntologyGraph` |
| How do I get a go.obo, once, on purpose? | `update_go("go.obo")` | `Loaders.UpdateGeneOntology` |

Every example in this section is executed when the docs are built, against output recorded from the
real bridge, so what you see is what pyMzLib returns. The data is real: mzLib's own copy of a
MetaMorpheus 1.1.11 protein-group table from PXD036557 (six groups: tubulin alpha-1B, ADP/ATP
translocase 2, the H2A.Z pair, 14-3-3 zeta, bovine serum albumin as a contaminant, and one decoy),
the five human UniProt entries it names, cut whole from UniProt's reviewed human proteome, and GO
release 2026-07-26. The go.obo committed for the tests is trimmed to the 412 terms these proteins
reach; against the full 48,340-term release the same call gives the same 563 rows.

### First, get a go.obo - once, and keep it

Terms are added, renamed and moved between GO releases, so a GO result means something only
relative to one release. `annotate_go()` therefore **never downloads anything**: it reads the
go.obo you name. Fetching one is its own function, which you call on purpose:

```pycon
>>> update = pymzlib.proteins.update_go("go.obo")
>>> update.go.release, update.go.term_count, update.changed
('releases/2026-07-26', 48340, True)

```

`update_go()` streams the whole file (about 37 MB) from GO's PURL, which always serves the current
release. Call it again later and a newer release replaces the file, while the old one is kept beside
it as `go.obo.<timestamp>`, so a result you already published can still be reproduced. Keep the
go.obo with your project, the way you keep the FASTA you searched.

### Annotate the groups

Give it the protein-group table, the UniProt XML the search used, and the go.obo:

```pycon
>>> go = pymzlib.proteins.annotate_go(
...     "PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
...     go_obo="go-pxd036557.obo", category_map="organelle_map.tsv")
>>> go.table_row_count, go.decoy_group_count, go.group_count, go.row_count
(6, 1, 5, 563)

```

Six rows in the table, one of them a decoy, which is skipped because decoys carry no GO: five groups,
563 rows. Any MetaMorpheus protein-group table works - `AllQuantifiedProteinGroups.tsv`,
`AllProteinGroups.tsv` from a search without quantification, or one file's `_ProteinGroups.tsv`.
Hand `go.columns` to `pandas.DataFrame` for a data frame; `go.records` gives one dict per row.

### Read one row

Both histones are annotated to the nucleosome, each on its own evidence:

```pycon
>>> nucleosome = next(r for r in go.records
...                   if r["protein_group"] == "P0C0S5|Q71UI9" and r["go_id"] == "GO:0000786")
>>> nucleosome["go_name"], nucleosome["aspect"], nucleosome["n_with"], nucleosome["n_members"]
('nucleosome', 'cellular_component', 2, 2)
>>> nucleosome["accession_used"], nucleosome["evidence_by_member"]
(['P0C0S5', 'Q71UI9'], {'P0C0S5': ['ECO:0000353'], 'Q71UI9': ['ECO:0000353']})
>>> nucleosome["propagated"], nucleosome["inherited"]
(False, False)

```

| column | what it tells you |
|---|---|
| `accession_used` | the members carrying the term - directly or through a more specific one |
| `accession_direct` | of those, the members annotated to this exact term |
| `n_with` of `n_members` | how many members carry it; `n_with == n_members` is the consensus |
| `evidence`, `evidence_by_member` | [ECO](https://www.evidenceontology.org/) codes, pooled and per member |
| `propagated` | `True`: no member is annotated to this term; it is implied by a more specific one |
| `inherited` | `True`: the members' terms were borrowed - an isoform (`P04406-2`) or a sequence variant (`P04406_A20T`) the database does not annotate takes its entry's terms |
| `entrapment_members` | the group's entrapment members, on every row of the group, so you decide whether they count |
| `q_value` | the group's q-value, from the table - **never filtered on** |

The column names are mzLib's own `GoAnnotationTsv` schema, so this table, the file `out=` writes
and the file mzLib writes inside a search all mean the same thing.

### Union, consensus, direct only: filters, not modes

The table is the **union** over a group's members. Everything narrower is one comparison:

```pycon
>>> histones = [r for r in go.records if r["protein_group"] == "P0C0S5|Q71UI9"]
>>> union = {r["go_id"] for r in histones}
>>> consensus = {r["go_id"] for r in histones if r["n_with"] == r["n_members"]}
>>> len(union), len(consensus)
(106, 57)
>>> euchromatin = next(r for r in histones if r["go_name"] == "euchromatin")
>>> euchromatin["accession_used"]
['P0C0S5']

```

Of the 106 terms either H2A.Z variant carries, 57 are carried by both. *Euchromatin* is not one of
them: only `P0C0S5` is annotated to it. Whether that should count for the group is a scientific
question, which is why the table answers it per member rather than deciding.

Propagation is what makes broad questions answerable. 14-3-3 zeta has 205 rows, and only 42 of them
are terms UniProt annotates directly; the other 163 are their ancestors:

```pycon
>>> zeta = [r for r in go.records if r["protein_group"] == "P63104"]
>>> len(zeta), sum(r["propagated"] is False for r in zeta)
(205, 42)

```

Keep `propagated` in mind before an enrichment test: most tools propagate terms themselves, and
feeding them already-propagated rows counts every ancestor twice. Give such a tool the rows with
`propagated is False`.

### A group with no term still gets a row

Bovine serum albumin is in the table as a contaminant. It gets exactly one row, with no term, and
a status that says why:

```pycon
>>> [(r["protein_group"], r["annotation_status"]) for r in go.records if r["go_id"] is None]
[('P02769', 'contaminant')]

```

| `annotation_status` | means |
|---|---|
| `annotated` | the row has a term (every row with a `go_id` is `annotated`) |
| `no_go_terms` | the group's members are in the database, but none has a GO annotation |
| `no_entry` | a member is not in the database at all - annotate against the database the search used |
| `contaminant` | the group is a contaminant and none of its members has a term |

So a group never silently disappears, and *no annotation* is never confused with *not looked up*.

### Count groups, not rows

`go.header` is the provenance header mzLib's writer puts at the top of the file, read back. Its
counters count **groups**, and only groups at q ≤ 0.01:

```pycon
>>> go.header["counter_q_value_max"]
'0.01'
>>> {k: go.header[k] for k in ("status_annotated", "status_no_go_terms", "status_no_entry", "status_contaminant")}
{'status_annotated': '4', 'status_no_go_terms': '0', 'status_no_entry': '0', 'status_contaminant': '1'}

```

The rows themselves are never filtered by q-value: every group is in the table, and you choose the
cut-off.

### Put terms into your own categories

GO has tens of thousands of terms; a figure usually wants a handful of categories. mzLib ships **no
vocabulary** - the categories are your science - but it applies yours with one rule: a term belongs
to a category when one of the category's anchors is the term itself or one of its ancestors, and
within a category the most specific anchor wins. A map is a small tab-separated file:

```text
#!category_map_format 1
#!map_name organelle
#!map_version 1
category	subcategory	anchor_go_id
nucleus		GO:0005634
nucleus	chromatin	GO:0000785
mitochondrion		GO:0005739
mitochondrion	inner_membrane	GO:0005743
cytoskeleton		GO:0005856
cytoskeleton	microtubule	GO:0005874
cytosol		GO:0005829
plasma_membrane		GO:0005886
extracellular		GO:0005576
```

The categories come back as their own table, joined to the annotations on `go_id`:

```pycon
>>> go.categories.map_name, go.categories.anchor_count, go.categories.row_count
('organelle', 9, 30)
>>> inner = next(c for c in go.categories.records if c["go_id"] == "GO:0005743")
>>> inner["category"], inner["subcategory"]
('mitochondrion', 'mitochondrion:inner_membrane')
>>> place = {}
>>> for c in go.categories.records:
...     place.setdefault(c["go_id"], set()).add(c["category"])
>>> where = {}
>>> for r in go.records:
...     for category in place.get(r["go_id"], ()):
...         where.setdefault(r["protein_group"], set()).add(category)
>>> sorted(where["P05141"])
['cytoskeleton', 'mitochondrion', 'nucleus', 'plasma_membrane']

```

A term under none of your anchors has no category row, so absence means "outside your map", not
"unknown". The map's name, version and sha256 travel with the result (`go.categories.sha256`): a
version is a claim, the hash is proof two runs used the same map.

### A large run: write the file, return the summary

A real search annotates hundreds of groups into tens of thousands of rows: the full MetaMorpheus
table from the same PXD036557 search gave 715 target groups and 87,513 rows, in about 20 s, against
the whole human proteome ([#1353][1353]). Write them to a file with mzLib's own writer and bring back
only the summary:

```pycon
>>> big = pymzlib.proteins.annotate_go(
...     "PXD036557_AllQuantifiedProteinGroups.tsv", "pxd036557_proteins.xml",
...     go_obo="go-pxd036557.obo", out="go_annotations.tsv", limit=0)
>>> big.written.path, big.written.row_count, big.returned_count
('go_annotations.tsv', 563, 0)

```

`limit` and `offset` shape only what comes back; the file always gets every row. It is a
tab-separated table whose header records what it was computed from:

```text
#!go_annotation_format 1
#!mzlib_version 1.0.0+0a808fec346e6e8f334e455490faab463ea65457
#!mzlib_release none
#!go_release releases/2026-07-26
#!go_obo_sha256 c1cdfd098c5cf395ce4bbefdea017fd698fdf22d135b28623e74f6f714e1c5fa
#!annotation_db_sha256 484e16a9a2fa5f4754f90d904f9030621c6ce3257488fd39fdff9b4fb8a0b5db
#!source_file_sha256 751a00c5109bdf350db75d2d4053e0fd78fea8326c5a38bce95b1ce0d8445ae6
#!counter_q_value_max 0.01
#!n_multi_member_groups 1
#!status_annotated 4
#!status_no_go_terms 0
#!status_no_entry 0
#!status_contaminant 1
protein_group	accession_used	accession_direct	...	go_id	go_name	aspect	...
```

`mzlib_release` reads `none` because the bridge builds mzLib from source; `mzlib_version` names the
exact commit. `out=` must end in `.tsv`, and is checked before anything is read, so a typo costs
milliseconds, not a proteome load.

### When the database is newer than the ontology

UniProt is released more often than you will move your go.obo, so sooner or later the database
cites a GO id your release does not have. By default that **fails**, naming every missing id,
because a silently dropped term is a silently wrong answer. If you would rather keep the run,
pass `skip_unknown_go_ids=True`: each such id is dropped, listed in `go.unresolved_go_ids`, counted
in the file header as `unresolved_go_ids`, and described in `go.caveats`. A protein whose only terms
were dropped then reads `no_go_terms`. The better fix is usually a newer go.obo.

### What it will not do, and says so

- **One annotation database.** Every row carries the sha256 of the database its terms came from, so
  two databases have no honest single value. Use the UniProt XML of the proteome you searched. A
  FASTA is refused outright: it carries no GO, and every group would read `no_go_terms` whatever the
  proteins are.
- **No hidden download.** A missing go.obo is an error that tells you to run `update_go()`.
- **No filtering.** Not by q-value, not by evidence, not by member. Every filter is a column.
- **Obsolete terms are kept** if your database still cites them; they have no ancestors in the
  release, so they never propagate.

To reproduce a GO result, record three hashes with it: `go.go.sha256` (the ontology),
`go.annotation_database.sha256` (the database) and `go.groups_file_sha256` (the table). The file
`out=` writes already carries all three.

## Resolve proteins to Ensembl genes, reproducibly

**The problem.** You want gene-level results - to join with RNA-seq, to count genes, to compare
with last year's analysis - and "the gene" has to mean the same thing each time. A gene *symbol*
cannot do that: symbols are renamed, one symbol can name several loci, and the `Gene` column of a
protein-group table is only the first name per protein. A stable Ensembl gene id counted against
one Ensembl release can.

`resolve_genes()` reads the Ensembl links the UniProt XML carries for each protein and counts them
against **a gene set you supply**: an Ensembl GTF for one release.

```pycon
>>> genes = pymzlib.proteins.resolve_genes(
...     ["human_subset.xml", "human_extra.fasta"],
...     contaminants=["contaminants.fasta"],
...     gtf="Homo_sapiens.GRCh38.116.gtf",                 # yours: the release's .gtf.gz
...     xref="Homo_sapiens.GRCh38.116.uniprot.tsv",        # optional second opinion
... )
>>> genes.gene_set.release, genes.gene_set.genome_build
('116', 'GRCh38.p14')
>>> genes.outcome_counts                                   # doctest: +NORMALIZE_WHITESPACE
{'resolved': 1, 'multi_gene': 0, 'off_primary_only': 1, 'not_in_source': 4,
 'unrecognized_accession': 0, 'contaminant_not_mapped': 1}
>>> for r in genes.records:
...     print(r["accession"], r["outcome"], r["gene_id"], r["gene_symbol"], r["ensembl_xref_agrees"])
P04406 resolved ENSG00000111640 GAPDH True
O43653 off_primary_only None None None
P02768 not_in_source None None None
Q13409 not_in_source None None None
Q13409-2 not_in_source None None None
Q13409-3 not_in_source None None None
P02769 contaminant_not_mapped None None None

```

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

```pycon
>>> calls = pymzlib.proteins.classify_peptides(
...     ["VGVNGFGR", "LVLNGNPLTLFQER", "ALSEQINIFFDYSGR", "YLYEIAR", "AEFVEVTK", "PEPTIDEK"],
...     ["human_subset.xml", "human_extra.fasta"],
...     contaminants=["contaminants.fasta"],
... )
>>> for r in calls.records:
...     print(r["peptide"], r["sharing"], r["accessions"])
VGVNGFGR Unique ['P04406']
LVLNGNPLTLFQER Unique ['P04406']
ALSEQINIFFDYSGR SharedWithinGene ['Q13409', 'Q13409-2', 'Q13409-3']
YLYEIAR SharedAcrossGenes ['P02768', 'P02769']
AEFVEVTK Unique ['P02769']
PEPTIDEK NotInDatabase []
>>> calls.records[2]["shared_gene_keys"]
['entry:Q13409', 'gene:Homo sapiens:DYNC1I2']

```

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

Results come back as JSON tables only: `read` has no `out=` option, because its three tables
would need three files.

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

For `annotate_go()` and `update_go()`:

| Raised | When |
|---|---|
| `UsageError` | a missing groups table, database, go.obo or category map; a FASTA database; an `out` or `categories_out` that is not `.tsv`, names an input, or sits in a folder that does not exist; `categories_out` without `category_map`; a negative `limit` or `offset` - all before any file is read. For `update_go()`, a folder that does not exist |
| `BridgeError` | the database cites a GO id the release lacks (type `InvalidDataException`, every id named) unless `skip_unknown_go_ids=True`; a go.obo, category map or table mzLib cannot parse; a group written twice with different members. For `update_go()`, GO's server is unreachable (type `ServiceUnavailable`: retry later) |

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
resources the answers come from:

--8<-- "docs/reference/_generated/cite.proteins.md"

The same verbs are specified once, language-neutrally, for pyMzLib, mzLibRust and mzLibR:
`proteins read`, `genes resolve`, `proteins classify-peptides`, `proteins annotate-go` and
`proteins update-go`. See the
[wire-verb reference](../reference/proteins.md) for every parameter and field with its unit, and the
[API reference](../reference.md#pymzlibproteins) for the Python classes.

[1336]: https://github.com/smith-chem-wisc/mzLib/pull/1336
[1338]: https://github.com/smith-chem-wisc/mzLib/pull/1338
[1348]: https://github.com/smith-chem-wisc/mzLib/pull/1348
[1353]: https://github.com/smith-chem-wisc/mzLib/pull/1353
[1366]: https://github.com/smith-chem-wisc/mzLib/pull/1366
