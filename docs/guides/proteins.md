# Protein databases

A search result tells you *which accession* a peptide matched. It does not tell you what organism
that accession is, what it does, which gene encodes it, or whether the peptide could have come from
somewhere else. **The protein database you searched already knows all four**, and
`pymzlib.proteins` reads them out of it with mzLib:

| You want to know | Call | mzLib |
|---|---|---|
| What organism, taxon, gene and mass is this accession? What are its GO terms? | [`read()`](#look-up-organism-and-taxon-for-a-list-of-accessions) | `ProteinDbLoader`, `Protein.GoTerms` ([#1336][1336]), `Protein.EnsemblGeneReferences` |
| Which Ensembl gene is this protein, in a way I can reproduce next year? | [`resolve_genes()`](#resolve-proteins-to-ensembl-genes-reproducibly) | `EnsemblGeneResolver` ([#1338][1338]) |
| Which mouse genes are orthologous to this human gene, and if none, why not? | [`build_orthology()`](#orthologs-across-species-an-ensembl-compara-snapshot) | `OrthologySnapshot`, `OrthologySnapshotWriter` ([#1381][1381]) |
| Does this peptide identify one protein, one gene, or neither? | [`classify_peptides()`](#decide-whether-a-peptide-is-unique-with-i-l) | `PeptideUniquenessClassifier` ([#1348][1348]) |
| Which GO terms does each protein group carry, with every member kept and the ontology's ancestors filled in? | [`annotate_go()`](#annotate-protein-groups-with-go-keeping-every-member) | `GoGroupAnnotator` ([#1353][1353]), `ToGoAnnotationGroups` ([#1366][1366]) |
| Where do I get a go.obo, and how do I keep the release pinned? | [`update_go()`](#first-get-a-goobo-once-and-keep-it) | `Loaders.UpdateGeneOntology` ([#1353][1353]) |
| Where is residue 120 of this protein in its isoform, or in the bovine protein? | [`align()`, `align_many()`](#map-a-residue-from-one-sequence-onto-another) | `PairwiseAligner` ([#1428][1428]) |
| Which entry, isoform and variants does `P04637-2_S70N` or `NP_000537.3_R72P` name? | [`parse_accessions()`](#read-an-accession-into-entry-isoform-and-variants) | `VariantApplication.ParseAccession` ([#1382][1382]) |

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
protein-group table and the UniProt entries it names. The orthology section uses mzLib's own
three-species test world (`pkg/python/tests/fixtures/orthology/`). Blocks titled **Not run** say why.

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

## Map a residue from one sequence onto another

**The problem.** A position means something only in the sequence it was counted in. A phosphosite
reported on the canonical isoform, a variant called against last year's UniProt release, a cysteine
you know in the human protein: before you look for it in an isoform, an ortholog or a newer entry,
you need the residue that *corresponds* to it there. Counting the same number along the other
sequence is wrong as soon as there is one insertion or deletion before it.

`align()` aligns two sequences with mzLib's `PairwiseAligner` ([#1428][1428]) and maps every residue
to its counterpart in the other sequence, or to `None` when it has none. Here the canonical DYNC1I2
sequence (Q13409) is aligned with its isoforms 2B and 2C, read from the database with
`sequences=True`:

```pycon
>>> db = pymzlib.proteins.read(
...     ["human_extra.fasta"], contaminants=["contaminants.fasta"],
...     accessions=["Q13409", "Q13409-2", "Q13409-3", "P02768", "P02769"], sequences=True)
>>> seq = dict(zip(db.columns["accession"], db.columns["sequence"]))
>>> isoforms = pymzlib.proteins.align_many([(seq["Q13409"], seq["Q13409-2"]),
...                                         (seq["Q13409"], seq["Q13409-3"])])
>>> two_b, two_c = isoforms.alignments
>>> for p in (75, 76, 82, 100, 120, 133):
...     print(p, two_b.map_source_to_target(p), two_c.map_source_to_target(p))
75 75 75
76 None None
82 76 76
100 94 94
120 114 None
133 127 107

```

Reading it row by row:

- **Up to residue 75 the numbers agree.** Both isoforms begin with the same 75 residues.
- **Residues 76 to 81 have no counterpart in either isoform**: `None`, never the nearest residue.
  Both isoforms lack those six residues of the canonical sequence.
- **After the gap the numbers are shifted**: canonical residue 100 is residue 94 of each isoform.
- **Isoform 2C lacks 20 more** (canonical 113 to 132), so canonical 120 has no counterpart in 2C,
  and canonical 133 is 2C's residue 107.

`map_target_to_source()` goes the other way. To see where the gaps fall, print the aligned rows,
which hold `-` wherever a residue faces a gap:

```pycon
>>> print(two_c.aligned_source[60:140], two_c.aligned_target[60:140], sep="\n")
EALLQSMGLTPESPIVFSEYWVPPPMSPSSKSVSTPSEAGSQDSGDGAVGSRTLHWDTDPSVLQLHSDSDLGRGPIKLGM
EALLQSMGLTPESPI------VPPPMSPSSKSVSTPSEAGSQDSGDGAVGSR--------------------RGPIKLGM
>>> two_b.identical, two_b.substituted, two_b.source_residues_facing_gap
(632, 0, 6)

```

**What each number means.** `identical` and `substituted` count the paired residues that are the
same letter and that differ; `source_residues_facing_gap` and `target_residues_facing_gap` count
the residues with no counterpart. `score` is the alignment's BLOSUM62 score in half-bits (see the
traps). `aligner_id` names the method and every parameter, so store it beside any position you map:

```pycon
>>> isoforms.aligner_id
'global-affine;BLOSUM62;open=11;extend=1;end-gaps=free'

```

### Orthologs: the counterpart can be a different amino acid

Across species most residues still pair, but many pair with a *different* residue. Human serum
albumin's free thiol, Cys34 of the mature protein, is residue 58 of the UniProt sequence, which
counts the 24-residue signal peptide and propeptide. In bovine albumin it is a cysteine at 58 too:

```pycon
>>> albumin = pymzlib.proteins.align(seq["P02768"], seq["P02769"])
>>> albumin.map_source_to_target(58), seq["P02769"][58 - 1]
(58, 'C')
>>> albumin.identical, albumin.substituted, albumin.source_residues_facing_gap
(465, 142, 2)
>>> albumin.map_source_to_target(26), seq["P02768"][26 - 1], seq["P02769"][26 - 1]
(26, 'A', 'T')

```

Residue 26 pairs, but it is A in human and T in bovine albumin. **Always check the residue at the
mapped position** before you carry a modification across: a phosphosite on S is not a site on A.
`align_many()` has the whole residue map as a table, one row per residue of both sequences, with
both letters side by side. Hand `isoforms.residue_map.columns` to `pandas.DataFrame` and compare
the `residue` and `mapped_residue` columns.

### Many pairs: one call

`align_many()` takes a list of `(source, target)` pairs and aligns them in one bridge process, in
order, so a few thousand isoform pairs cost one process start, not a few thousand. Its
`residue_map` table has a `pair_index` column, and `alignments` gives one `Alignment` per pair. For a
large batch where you only want the scores and counts, pass `residue_map=False`: the map has one row
per residue of both sequences, and leaving it out keeps the answer small.

### The options, and whose defaults they are

Every option you leave out is **mzLib's to set**. pyMzLib sends nothing for it, and the result
says what the aligner used. The isoform alignment above passed no options, so these are mzLib's
own defaults in this build:

```pycon
>>> isoforms.matrix, isoforms.gap_open, isoforms.gap_extend, isoforms.end_gaps, isoforms.max_cells
('BLOSUM62', 11, 1, 'free', 200000000)

```

| Argument | Meaning |
|---|---|
| `gap_open` | cost of opening a gap, in BLOSUM62's half-bits; a gap of *k* residues costs `gap_open + k * gap_extend` (BLAST's convention; mzLib's defaults are BLASTP's) |
| `gap_extend` | cost of each residue in a gap, the first included |
| `end_gaps` | `"free"`: an overhang at either end costs nothing (EMBOSS needle's default), so a signal peptide or a longer tail is not penalised. `"charged"`: every gap costs |
| `max_cells` | the largest table mzLib will build, `(len(source) + 1) * (len(target) + 1)` cells of one byte; a larger pair is refused rather than allowed to exhaust memory |

The matrix is always BLOSUM62, mzLib's default and the only one offered.

```pycon
>>> charged = pymzlib.proteins.align("MLLLPEPTIDE", "PEPTIDE", end_gaps="charged")
>>> charged.aligned_target, charged.score, charged.aligner_id
('----PEPTIDE', 24, 'global-affine;BLOSUM62;open=11;extend=1;end-gaps=charged')

```

With free end gaps the same pair scores 39: the four-residue overhang costs `11 + 4 * 1` only when
it is charged.

### The traps

- **`score` is not an identity, a bit score or an E-value.** It is a raw sum of BLOSUM62 scores in
  half-bits, less the gaps, and it grows with length. Compare scores only between alignments with
  the same `aligner_id`. pyMzLib computes no percent identity, because there are several
  definitions: divide `identical` by the length you mean.
- **It is a global alignment.** Every residue of both sequences is placed, paired or facing a gap.
  With free end gaps a short sequence is placed on the best-scoring stretch of a long one, even when
  it has no real match there.
- **One optimal alignment, chosen the same way every time.** Where two alignments score the same,
  mzLib's traceback picks one, always the same one. Isoform 2B above is a real case: canonical V76
  and V82 are both V, either could pair with the isoform's V76, and mzLib pairs V82, leaving V76 to
  W81 facing the gap. Another aligner may give the equally good other answer.
- **Pass bare one-letter sequences.** mzLib refuses only `-` and `.`, as a `UsageError`. Any other
  character BLOSUM62 does not name - selenocysteine `U`, pyrrolysine `O`, a digit, a space, a
  modification bracket - is scored as an unknown residue `X`, silently. Strip modifications first.
- **Long pairs need a larger `max_cells`.** Time and memory grow with the product of the two
  lengths. Two titin isoforms need about 1.2 billion cells, past the default; mzLib accepts up to
  2,147,483,591. Two identical sequences are aligned residue for residue without the table, so the
  limit does not apply to them.

## Read an accession into entry, isoform and variants

**The problem.** The accessions a search reports are not all plain entries. UniProt isoforms carry
`-2`; RefSeq accessions carry a `.3` version; and mzLib names a protein with applied sequence
variants after its entry, `P04637_S70N_A80T` or `NP_000537.3_R72P`. To join results to anything you
keep per entry - genes, orthologs, GO terms - you need the entry back, without losing the isoform or
the variants. Cutting at the first `_` turns `NP_000537` into `NP`.

`parse_accessions()` reads each accession with mzLib's own grammars
(`VariantApplication.ParseAccession`, [#1382][1382]):

```pycon
>>> parsed = pymzlib.proteins.parse_accessions(["P04637", "P04637-2", "Q13409-2",
...     "NP_000537.3_R72P", "P04637_S70N_A80T", "P04637_2", "DECOY_P02768",
...     "sp|P04637|P53_HUMAN", "P04637|P02768"])
>>> for r in parsed.records:
...     print(r["input"], r["namespace"], r["entry_accession"], r["isoform"], r["version"],
...           r["applied_variants"])
P04637 uniprot P04637 None None None
P04637-2 uniprot P04637 2 None None
Q13409-2 uniprot Q13409 2 None None
NP_000537.3_R72P refseq NP_000537 None 3 R72P
P04637_S70N_A80T uniprot P04637 None None S70N_A80T
P04637_2 unrecognized P04637_2 None None None
DECOY_P02768 unrecognized DECOY_P02768 None None None
sp|P04637|P53_HUMAN unrecognized sp|P04637|P53_HUMAN None None None
P04637|P02768 unrecognized P04637|P02768 None None None
>>> parsed.namespace_counts
{'unrecognized': 4, 'uniprot': 4, 'refseq': 1}

```

What each column means:

- **`entry_accession`** is the UniProt entry or the unversioned RefSeq accession. It answers
  "same entry?", **not** "same sequence?": `P04637` and `P04637-2` share it, but they are different
  sequences. Join per-entry data on it; join sequence-level data on `input`.
- **`isoform`** and **`version`** are the `-n` and `.n` numbers, `None` when there is none.
- **`parent_accession`** is the entry the variants were applied to, isoform or version kept:
  `NP_000537.3` for `NP_000537.3_R72P`. Look the parent sequence up under it.
- **`applied_variants`** is the variant suffix exactly as written, `None` when there is none. It is
  not checked against any sequence.
- **`namespace`** is `uniprot`, `refseq` or `unrecognized`, spelled as `resolve_genes()` spells it,
  so the two tables join.

`parsed.entry_of()` is the `input` to `entry_accession` dictionary, for grouping.

**It parses and never repairs.** Everything outside the two grammars comes back `unrecognized`,
with `entry_accession` equal to the input, so a data problem stays visible instead of being fixed
silently:

- **`P04637_2` is not a variant of P04637.** It is mzLib's loader counter for a *different* entry
  whose accession collided with P04637's.
- **Decoy, contaminant and entrapment prefixes** (`DECOY_P02768`) and **FASTA headers**
  (`sp|P04637|P53_HUMAN`) are left as they are.
- **A protein group** (`P04637|P02768`) is not split. Split groups on `|` first, if that is what you
  mean.
- **Lower case and spaces count.** `" P04637"` is unrecognized: pyMzLib sends every accession
  exactly as you give it. Only RefSeq *protein* accessions (`NP_`, `XP_`, `YP_`, `WP_`, `AP_`) are
  recognized, so a transcript (`NM_`) or an Ensembl id (`ENSP...`) is unrecognized too.
- **An isoform or version past 2,147,483,647** reads as `None`, while the namespace still says
  `uniprot` or `refseq`: mzLib keeps the entry and drops the number.

## Orthologs across species: an Ensembl Compara snapshot

**The problem.** You have a mouse result and a human question, or identifications from three species
to line up: which human gene is this mouse gene's ortholog? Ensembl Compara answers that, per
species pair, in bulk files with three traps:

1. **The dumps partition.** Each species has its own homology dump, and Compara puts each of a
   pair's relationships in *one* of the pair's two dumps, arbitrarily. In release 116 every
   human-mouse orthology is in the **mouse** dump and none is in the human one. Reading the dump
   named after your species gives a complete-looking, badly incomplete answer.
2. **Orthology is not transitive.** Human A ~ mouse B and mouse B ~ rat C does not make A ~ C.
   Chaining pairs into groups invents tuples Compara never asserted.
3. **"No ortholog" means three things.** The gene's tree may hold a gene of the other species with no
   ortholog called, hold none of that species, or the gene may be in no tree at all.

`build_orthology()` builds mzLib's orthology snapshot (`OrthologySnapshot`, [#1381][1381]), which
makes the first trap structural - **mzLib refuses a call without every species' dump** - and keeps
the other two visible. This example is mzLib's own test world: four human genes, three mouse and two
rat, with every human-mouse row in the mouse dump:

```pycon
>>> gtfs = {"homo_sapiens": "Homo_sapiens.GRCh38.116.gtf",
...         "mus_musculus": "Mus_musculus.GRCm39.116.gtf",
...         "rattus_norvegicus": "Rattus_norvegicus.GRCr8.116.gtf"}
>>> dumps = [f"homologies/{s}/Compara.116.protein_default.homologies.tsv.gz" for s in gtfs]
>>> snap = pymzlib.proteins.build_orthology(
...     "116", gtfs=gtfs, dumps=dumps,
...     gene_trees="vertebrates.GeneTree_content.default.e116.txt.gz",
...     pair_status=("homo_sapiens", "mus_musculus"),     # every human gene, in mouse
... )
>>> for r in snap.pair_status.records:
...     print(r["gene_id"], r["gene_biotype"], r["status"])
ENSG00000000001 protein_coding has_ortholog
ENSG00000000002 protein_coding has_ortholog
ENSG00000000003 protein_coding tree_lacks_target_species
ENSG00000000004 protein_coding not_in_any_tree
>>> snap.snapshot_id[:12]
'5487ebcabe48'

```

Record `snap.snapshot_id` with any result you draw from the snapshot. It is mzLib's sha256 over the
format, the release, the species and every input's sha256, so the same inputs always give the same
id, and a different id means different inputs.

### What to download, once per release

Nothing is downloaded for you. For release 116, from `https://ftp.ensembl.org/pub/release-116/`:

| Input | Path under the release | Size (human, 116) | Argument |
|---|---|---|---|
| Each species' primary-assembly GTF | `gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz` | 141 MB | `gtfs={species: path}` |
| Each species' homology dump | `tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz` | 109 MB | `dumps=[...]` |
| The gene trees, one file for all species | `compara/vertebrates.GeneTree_content.default.e116.txt.gz` | 79 MB | `gene_trees=` |

- **The species names are Ensembl production names**: `homo_sapiens`, `mus_musculus`,
  `rattus_norvegicus`, exactly as the dumps' `species` column spells them. They are the keys of
  `gtfs=`, and they are the snapshot's species.
- **Take the dump in each species' own folder.** Every species' dump has the same file name, so keep
  them in their folders. One folder up, `homologies/` also holds a 10 GB file of the same name; it
  is not one species' dump.
- **Keep Ensembl's file names.** mzLib holds the gene trees and every dump to your `release` only when
  their names carry it (`.e116.txt.gz`, `Compara.116.`). A renamed file is not checked, and a caveat
  says which.
- **A GTF's own number is not compared.** Ensembl publishes some species with their own annotation
  release (yeast is `Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz` in release 116). mzLib holds a gene
  set to the dumps by content instead: a homology that names a gene the gene set lacks is refused.
- **A compact gene table works too.** If you already made one for [`resolve_genes()`](#resolve-proteins-to-ensembl-genes-reproducibly)
  with mzLib's `EnsemblGeneSetWriter`, pass it in `gene_sets={species: path}` instead of `gtfs=`. It
  carries the GTF's provenance, so the snapshot and its id are the same.

### Every species' dump, always

The result shows where the rows came from. The human dump holds one row (human-rat); all three
human-mouse rows came from the mouse dump:

```pycon
>>> {species: dump.rows_kept for species, dump in snap.dumps.items()}
{'homo_sapiens': 1, 'mus_musculus': 5, 'rattus_norvegicus': 1}
>>> for p in snap.pairs:
...     print(p.species_a, p.species_b, p.ortholog_count, p.paralog_count)
homo_sapiens homo_sapiens 0 0
homo_sapiens mus_musculus 3 0
homo_sapiens rattus_norvegicus 1 0
mus_musculus mus_musculus 0 1
mus_musculus rattus_norvegicus 2 0
rattus_norvegicus rattus_norvegicus 0 0
>>> snap.pair("mus_musculus", "homo_sapiens").row_count     # either order
3

```

Leave the mouse dump out and mzLib refuses: *"No homology dump for mus_musculus. Each pair's
relationships are split between its two species' dumps, so every species' dump is required."* That
is a `UsageError`. A pair with nothing in it, such as the human paralogs here, still has its entry.

### What each status means

`pair_status=(from_species, to_species)` gives every gene of `from_species`'s gene set exactly one
status in `to_species`. It says how far Compara's inference got, not what biology holds:

| status | what it means | read it as |
|---|---|---|
| `has_ortholog` | at least one ortholog in `to_species` | an answer; the Parquet has the rows |
| `no_edge_in_shared_tree` | the gene's tree holds a `to_species` gene, and Compara called no ortholog to it | the nearest the data comes to "no ortholog exists" |
| `tree_lacks_target_species` | the gene's tree holds no `to_species` gene | nothing was compared |
| `not_in_any_tree` | the gene is in the gene set but in no gene tree | no call was attempted |

`snap.pair_status.status_counts` counts genes per status, every status present, and
`snap.pair_status.genes("no_edge_in_shared_tree")` lists one status's genes. **A status is not
symmetric**: ask for `("mus_musculus", "homo_sapiens")` separately. **Every biotype is a row**, so
choose your denominator by filtering on `gene_biotype`. mzLib's own human-to-mouse check on release
116 counted protein-coding plus IG and TR genes, which is not the same set as `protein_coding`
alone.

### Write the snapshot, then query it

The relationships themselves go to disk. Give `out=` a directory that does not exist yet:

```pycon
>>> written = pymzlib.proteins.build_orthology(
...     "116", gtfs=gtfs, dumps=dumps,
...     gene_trees="vertebrates.GeneTree_content.default.e116.txt.gz",
...     restrict_gene_trees=True, out="compara-116",
... )
>>> written.snapshot_id == snap.snapshot_id
True
>>> for f in written.written.files[:7]:
...     print(f.path, f.rows)
genes/homo_sapiens.parquet 4
genes/mus_musculus.parquet 3
genes/rattus_norvegicus.parquet 2
members/homo_sapiens.parquet 3
members/mus_musculus.parquet 3
members/rattus_norvegicus.parquet 2
pairs/homo_sapiens__homo_sapiens.parquet 0
>>> written.written.files[-1].path, written.written.files[-1].rows
('views.sql', None)

```

The directory holds `genes/SPECIES.parquet` (the gene set), `members/SPECIES.parquet` (gene to
Compara tree), `pairs/A__B.parquet` (every relationship between A and B, A sorting first; A == B
holds the paralogs), `views.sql` and `manifest.json` (every input and file with its sha256). The
same inputs always give the same bytes. **A directory that is not empty is refused**: a snapshot is
written once, and next release's goes in a new directory to compare with this one.

`restrict_gene_trees=True` keeps only your species' genes from the gene-tree file. Without it,
mzLib's default keeps every gene in it - 4.2 million for the vertebrates of release 116 - and a
caveat says how many. The snapshot is the same either way, as the equal ids show.

**mzLib writes Parquet and ships no reader, and pyMzLib reads none** (it has no third-party
dependency). Query the files with DuckDB, which runs the macros mzLib writes in `views.sql`:

```python title="Not run: needs DuckDB (pip install duckdb), which pyMzLib does not depend on"
import duckdb

con = duckdb.connect()
con.execute(open("compara-116/views.sql", encoding="utf-8").read())

# Every human protein-coding gene's status in mouse: the same as pair_status=, from the files.
con.sql("""
    SELECT status, count(*) AS genes
    FROM pair_status('compara-116', 'homo_sapiens', 'mus_musculus')
    WHERE gene_biotype = 'protein_coding'
    GROUP BY status ORDER BY status
""").fetchall()
# [('has_ortholog', 2), ('not_in_any_tree', 1), ('tree_lacks_target_species', 1)]

# The orthologs, oriented from human: one row per relationship, one-to-many kept as several rows.
con.sql("SELECT gene, partner_gene, relationship_type, identity, partner_identity "
        "FROM orthologs('compara-116', 'homo_sapiens', 'mus_musculus') ORDER BY gene, partner_gene").fetchall()
# [('ENSG00000000001', 'ENSMUSG00000000001', 'ortholog_one2one', 91.0, 90.0),
#  ('ENSG00000000002', 'ENSMUSG00000000002', 'ortholog_one2many', 71.0, 70.0),
#  ('ENSG00000000002', 'ENSMUSG00000000003', 'ortholog_one2many', 61.0, 60.0)]

# One gene per species where all three pairs are orthologs; never chained.
con.sql("SELECT * FROM species_set3('compara-116', 'homo_sapiens', 'mus_musculus', "
        "'rattus_norvegicus')").fetchall()
# [('ENSG00000000001', 'ENSMUSG00000000001', 'ENSRNOG00000000001', True)]
```

H2 is orthologous to mouse M2, and M2 to rat R2, but H2 ~ R2 was never called, so `species_set3`
does not return H2. That is trap 2, handled. Or read a table with pyarrow:

```python title="Not run: needs pyarrow (pip install pyarrow), which pyMzLib does not depend on"
import pandas as pd
import pyarrow.parquet as pq

pairs = pq.read_table("compara-116/pairs/homo_sapiens__mus_musculus.parquet",
                      filters=[("relationship_class", "=", "ortholog")])
df = pairs.to_pandas(types_mapper=pd.ArrowDtype)   # goc_score stays an integer; a null stays <NA>
```

**Keep the nulls.** A plain `pd.read_parquet()` turns an integer column with a null, such as
`goc_score`, into floats with `NaN`, and `is_high_confidence`'s null into `None` beside `True` and
`False`. `types_mapper=pd.ArrowDtype` keeps both as they are. In R, the `arrow` and `duckdb`
packages read the same files.

### What the Parquet holds

`pairs/A__B.parquet` has one row per relationship, ordered by Compara's `homology_id`. Values are
Compara's, verbatim, and **a NULL stays null, never 0 or false**:

| column | unit | meaning; null means |
|---|---|---|
| `relationship_type` | | Compara's `homology_type`: `ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`, `within_species_paralog`, ... |
| `relationship_class` | | `ortholog`, `paralog` or `homoeolog` |
| `gene_a`, `protein_a`, `gene_b`, `protein_b` | | side A is the species that sorts first; each value moves with its gene |
| `identity_a`, `identity_b` | percent | the share of **that side's** protein identical to the other's in the alignment, so the two differ; null: Compara gave none |
| `goc_score` | 0 to 100 | gene order conservation; null: not computed for this pair |
| `wga_coverage` | percent | whole-genome alignment coverage; null: no alignment for this pair |
| `dn`, `ds` | substitutions per site | Compara's dN and dS; null for most comparisons |
| `is_high_confidence` | | Compara's flag; null: Compara declines to say, which is not `False` |
| `source_dump` | | the **file name** of the dump the row came from: provenance, not meaning |

`members/SPECIES.parquet`'s `group_id` is `compara-116:` plus the tree id, so trees are never
compared across releases. The full column lists, with types, are in the
[wire-verb reference](../reference/proteins.md#build_orthology).

### The traps

- **Never chain pairs.** Use `species_set3` (or mzLib's `OrthologySnapshot.SpeciesSet` for more
  species), which requires every two members of a tuple to be orthologs.
- **One-to-many is several rows.** Nothing is picked: a human gene with two mouse orthologs has two
  rows, both `ortholog_one2many`. Decide per question.
- **The gene set and the dumps must be one release.** A homology naming a gene its gene set lacks is
  refused as a `BridgeError` (`InvalidDataException`): that is what a GTF and dumps from different
  releases look like. Use the primary-assembly GTF, as for `resolve_genes()`.
- **`source_dump` cannot tell species apart.** Ensembl gives every species' dump the same file name,
  so the column holds that one name everywhere, and a caveat says so. `manifest.json`'s `inputs`
  list each dump with its species and sha256.
- **Join with `resolve_genes()` on the gene set.** Its `gene_set_sha256` column and
  `snap.gene_sets[species].sha256` are the same value when both were counted against the same GTF.
  Check that before joining a resolution to a snapshot.
- **It is slow and large.** Every input is read whole. For human, mouse and rat at release 116, mzLib
  reports about 20 s to build and 40 s more to write 13 files (9.0 MB). Pass `timeout=None` (the
  default) and use `restrict_gene_trees=True` when memory is short.

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

For `build_orthology()`:

| Raised | When |
|---|---|
| `UsageError` | no species, no dumps, or a species in both `gtfs` and `gene_sets`; a missing GTF, gene table, dump or gene-tree file; `out` naming a file or a directory that is not empty; `pair_status` not two different species of the snapshot (all before any input is read). From mzLib: a species without a dump, two dumps of one species, a dump of a species not given, or a dump or gene-tree file named for another release |
| `BridgeError` | a file mzLib's reader refuses (type `InvalidDataException`, naming the file and line): a dump whose columns are not Compara's, an unknown `homology_type`, a gene-tree row without four cells. Or data that breaks the snapshot's rules: a homology naming a gene its gene set lacks, one gene in two species' gene sets |

For `align()`, `align_many()` and `parse_accessions()`, every error is a `UsageError`:

| Function | When |
|---|---|
| `align()`, `align_many()` | no pairs; a pair that is not two strings; a sequence with a tab or line break; two empty sequences; an option of the wrong type, or `end_gaps` not `"free"` or `"charged"`. From mzLib: a gap cost or `max_cells` out of range, a sequence holding `-` or `.`, or a pair whose table exceeds `max_cells` (the message names the pair, and nothing is returned). A `position` outside the sequence in `map_source_to_target()` |
| `parse_accessions()` | no accessions; one that is not a string, is blank, or holds a line break; one string where a list was expected |

## Cite

If this guide's results go into a paper, cite mzLib (see [Citing](../index.md#citing)) and the
resources the answers come from:

--8<-- "docs/reference/_generated/cite.proteins.md"

The same verbs are specified once, language-neutrally, for pyMzLib, mzLibRust and mzLibR:
`proteins read`, `genes resolve`, `genes orthology`, `proteins classify-peptides`,
`proteins annotate-go`, `proteins update-go`, `proteins align` and `proteins parse-accession`. See the
[wire-verb reference](../reference/proteins.md) for every parameter and field with its unit, and the
[API reference](../reference.md#pymzlibproteins) for the Python classes.

[1336]: https://github.com/smith-chem-wisc/mzLib/pull/1336
[1338]: https://github.com/smith-chem-wisc/mzLib/pull/1338
[1348]: https://github.com/smith-chem-wisc/mzLib/pull/1348
[1353]: https://github.com/smith-chem-wisc/mzLib/pull/1353
[1366]: https://github.com/smith-chem-wisc/mzLib/pull/1366
[1381]: https://github.com/smith-chem-wisc/mzLib/pull/1381
[1382]: https://github.com/smith-chem-wisc/mzLib/pull/1382
[1428]: https://github.com/smith-chem-wisc/mzLib/pull/1428
